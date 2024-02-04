# 
# gauge.py                                                               
# 
# D. Clarke
# 
# Tools for reading and writing gauge configurations in Python. To interact with binary configurations, we can use
# Python's built-in struct module. More info on that: https://docs.python.org/3.7/library/struct.html. Inspired by
# QCDUtils https://github.com/mdipierro/qcdutils.
#

import struct
import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import numbaON
from latqcdtools.math.math import rel_check
from latqcdtools.base.check import checkType
numbaON()
from latqcdtools.math.SU3 import SU3
from latqcdtools.physics.gauge import gaugeField


class confReader:

    """ Base class for reading configurations. """

    def __init__(self, Ns, Nt, nproc=None):
        """ Initialize a confReader object.

        Args:
            Ns (int): Spatial extension of lattice. 
            Nt (int): Euclidian time extension of lattice. 
            nrows (int): Number of saved rows for compression. 
        """        
        self.Ns = Ns
        self.Nt = Nt
        checkType(Ns,int)
        checkType(Nt,int)
        if nproc==None:
            self.nproc=self.Nt
        else:
            checkType(nproc,int)
            self.nproc=nproc
        self.offset = None     # How many bytes is the header?
        self.endianness = '>'  # Big endian by default
        self.precision = 'f'   # Single precision by default
        self.file = None       # File handle
        self.linkTrace = None  # <tr U>
        self.plaquette = None  # <tr U^[]>
        self.nrows = None      # number of rows
        self.gauge = gaugeField(self.Ns,self.Nt,self.nproc)


    def unpack(self, data) -> SU3():
        """ Unpack a string of bytes from file into an SU3 object. """
        numData = int(len(data)/self.getByteSize())
        linkTuple = struct.unpack( self.endianness + str(numData) + self.precision, data )
        link = SU3()
        i = 0
        for j in range(self.rows):
            for k in range(3):
                link[j,k] = complex( linkTuple[i], linkTuple[i+1] )
                i += 2
        if self.rows==2:
            link.su3unitarize()
        return link


    def pack(self, items):
        """ Packs an SU3 object into a string of bytes. """
        numData = len(items)
        return struct.pack(self.endianness + str(numData) + self.precision, *items)


    def getByteSize(self):
        if self.precision == 'f':
            return 4
        elif self.precision == 'd':
            return 8
        else:
            logger.TBError('Unknown precision',self.precision,'(expected f or d)')



class NERSCReader(confReader):

    def readHeader(self,fileName):

        """ Extract metadata from the header. Do a few checks to make sure the read-in went smoothly. This method
        assumes all NERSC headers have the minimal entries
            BEGIN_HEADER
            DATATYPE
            DIMENSION_1
            DIMENSION_2
            DIMENSION_3
            DIMENSION_4
            FLOATING_POINT
            END_HEADER
        """

        self.file = open(fileName,'rb')

        # Extract offset.
        header = self.file.read()
        minOffset = 78 # Minimal number of characters for the entries above, not including what comes on RHS.
        self.offset = header.find(b'END_HEADER\n')+11
        if self.offset < minOffset:
            logger.TBError(fileName,'not in NERSC format.')

        # Convert header to metaData dictionary.
        entries = header[:self.offset-1].split(b'\n')[1:-1]
        metaData = {}
        for entry in entries:
            LHS = entry.split(b'=')[0].strip()
            RHS = entry.split(b'=')[1].strip()
            metaData[LHS] = RHS

        # Check that the header was read correctly.
        Nx = int(metaData[b'DIMENSION_1'])
        Ny = int(metaData[b'DIMENSION_2'])
        Nz = int(metaData[b'DIMENSION_3'])
        Nt = int(metaData[b'DIMENSION_4'])
        if Nx != self.Ns:
            logger.TBError('Read Nx =',Nx,'. Expected ',self.Ns)
        if Ny != self.Ns:
            logger.TBError('Read Ny =',Ny,'. Expected ',self.Ns)
        if Nz != self.Ns:
            logger.TBError('Read Nz =',Nz,'. Expected ',self.Ns)
        if Nt != self.Nt:
            logger.TBError('Read Nt =',Nt,'. Expected ',self.Nt)

        # Extract endianness
        if metaData[b'FLOATING_POINT'].endswith(b'SMALL'):
            self.endianness = '<'
        elif metaData[b'FLOATING_POINT'].endswith(b'BIG'):
            self.endianness = '>'
        else:
            logger.TBError('Unrecognized endianness.')

        # Extract precision
        if metaData[b'FLOATING_POINT'].startswith(b'IEEE32'):
            self.precision = 'f'
        elif metaData[b'FLOATING_POINT'].startswith(b'IEEE64'):
            self.precision = 'd'
        else:
            logger.TBError('Unrecognized precision.')

        # Extract number of rows
        if metaData[b'DATATYPE'].endswith(b'3x3'):
            self.rows = 3
        else:
            self.rows = 2

        # Extract trace of average link
        try:
            self.linkTrace = float( metaData[b'LINK_TRACE'] )
        except:
            pass

        # Extract average plaquette
        try:
            self.plaquette = float( metaData[b'PLAQUETTE'] )
        except:
            pass


    def readConf(self,fileName):
        """ Read in the configuration. Check what you find in the binary against its metadata. """

        self.readHeader(fileName)

        bytesPerLink = 3*self.rows*2*self.getByteSize()

        for t in range(self.Nt):
            for z in range(self.Ns):
                for y in range(self.Ns):
                    for x in range(self.Ns):
                        for mu in range(4):
                            byteIndex = self.offset + bytesPerLink*( mu + 4*( x + self.Ns*( y + self.Ns*( z + self.Ns*t ) ) ) )
                            self.file.seek(byteIndex)
                            data = self.file.read(bytesPerLink)
                            link = self.unpack(data)
                            link.su3unitarize()
                            self.gauge.setLink(link,x,y,z,t,mu)
        logger.details('Loaded conf into gaugeField.')

        # Check <tr U>
        if self.linkTrace is not None:
            linkTrace = self.gauge.getLinkTrace()
            if not rel_check(linkTrace, self.linkTrace):
                logger.TBError('<tr U> is wrong. Compare:',linkTrace,'with',self.linkTrace)
            logger.details('Configuration',fileName,'has correct <tr U>.')

        # Check <tr U^[]>
        if self.plaquette is not None:
            plaq = self.gauge.getPlaquette()
            if not rel_check(plaq, self.plaquette):
                logger.TBError('<tr U^[]> is wrong. Compare:',plaq,'with',self.plaquette)
            logger.details('Configuration', fileName, 'has correct <tr U_plaq>.')

        return self.gauge
