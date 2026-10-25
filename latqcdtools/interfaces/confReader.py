# 
# confReader.py                                                               
# 
# D. Clarke, K. Ebira
# 
# Tools for reading gauge configurations in Python. To interact with binary configurations, we can use
# Python's built-in struct module. More info on that: https://docs.python.org/3.7/library/struct.html. NERSC reader inspired
# by https://github.com/mdipierro/qcdutils. ILDG reader from https://github.com/nftqcd/nthmc/blob/master/lib/fieldio.py
#

import numpy as np
import struct, os
import latqcdtools.base.logger as logger
from latqcdtools.base.speedify import numbaON,DEFAULTTHREADS
from latqcdtools.math.math import rel_check
from latqcdtools.base.check import checkType
from latqcdtools.interfaces.lime import LIMEMAGIC,trimNull,xmlFind,scidacChecksum
numbaON()
from latqcdtools.math.SU3 import SU3, fastUnitarize, lam
from latqcdtools.physics.gauge import gaugeField


DYNISOMILCMAGIC  = 20103 # Magic number at the start of dyniso's MILC-format files
DYNISOASCIIMAGIC = 59354 # Version number at the start of dyniso's ascii files


def nerscChecksum(data, precision) -> int:
    """
    NERSC checksum: the sum, modulo 2^32, of the link data viewed as unsigned 32-bit words. The words are
    formed from the values of the floats, not from their byte order on disk, so files with different
    endianness but the same links have the same checksum. This matches SIMULATeQCD.

    Args:
        data (array-like): complex link entries in the order they appear in the file 
        precision (str): 'f' for IEEE32 or 'd' for IEEE64 

    Returns:
        int: checksum 
    """
    if precision == 'f':
        dtype = '>c8'
    elif precision == 'd':
        dtype = '>c16'
    else:
        logger.TBRaise('Unrecognized precision',precision)
    words = np.frombuffer(np.asarray(data).astype(dtype).tobytes(), dtype='>u4')
    return int(np.sum(words, dtype=np.uint64) & 0xffffffff)


def detectConf(fileName):
    """
    Detect the format and lattice size of a configuration from its contents. NERSC files start with BEGIN_HEADER,
    ILDG files with the LIME magic number, and dyniso files are one of comp, MILC, or ascii, as described in
    DynisoReader. All of these store the lattice dimensions near the start of the file.

    Args:
        fileName (str): path to the configuration

    Returns:
        tuple: (format, Ns, Nt), where format is 'NERSC', 'ILDG', or 'dyniso'
    """
    checkType(str,fileName=fileName)
    fileSize = os.path.getsize(fileName)
    with open(fileName,'rb') as f:
        start = f.read(512)

    if start.startswith(b'BEGIN_HEADER'):
        fmt  = 'NERSC'
        dims = []
        for line in start.split(b'\n'):
            if b'=' in line:
                key = line.split(b'=')[0].strip().upper()
                if key in (b'DIMENSION_1', b'DIMENSION_2', b'DIMENSION_3', b'DIMENSION_4'):
                    dims.append(int(line.split(b'=')[1]))
        return fmt, *_cubicDims(dims, fileName)

    if len(start) >= 4 and struct.unpack('>i',start[:4])[0] == LIMEMAGIC:
        return 'ILDG', *_cubicDims(_ildgDims(fileName), fileName)

    # dyniso binary formats are in the byte order of the machine that wrote them, so try both.
    for endianness in ['<','>']:
        ints = struct.unpack(endianness+'5i',start[:20])
        if ints[0] == DYNISOMILCMAGIC:
            return 'dyniso', *_cubicDims(list(ints[1:5]), fileName)
        dims = list(ints[0:4])
        if min(dims) > 0 and fileSize == 16 + 32*4*dims[0]*dims[1]*dims[2]*dims[3]:
            return 'dyniso', *_cubicDims(dims, fileName)

    try:
        tokens = start.decode('ascii').split()
        if int(tokens[0]) == DYNISOASCIIMAGIC:
            dims = tokens[1:5]
        else:
            dims = tokens[0:4]
        intDims = []
        for dim in dims:
            intDims.append(int(dim))
        return 'dyniso', *_cubicDims(intDims, fileName)
    except (UnicodeDecodeError, ValueError, IndexError):
        logger.TBRaise('Could not detect the format of',fileName)


def _cubicDims(dims, fileName):
    """
    Turn [Nx, Ny, Nz, Nt] into (Ns, Nt), checking that the spatial volume is cubic.
    """
    if len(dims) != 4:
        logger.TBRaise('Could not find 4 lattice dimensions in',fileName)
    if not (dims[0] == dims[1] == dims[2]):
        logger.TBRaise('Only Ns^3 x Nt lattices are supported, but',fileName,'has dimensions',dims)
    return dims[0], dims[3]


def _ildgDims(fileName):
    """
    Get [Nx, Ny, Nz, Nt] from the scidac-private-file-xml or ildg-format record of a LIME file.
    """
    with open(fileName,'rb') as f:
        while True:
            header = f.read(144)
            if len(header) < 144:
                break
            magi, _, _, size, rtype = struct.unpack('>ihHq128s',header)
            if magi != LIMEMAGIC:
                logger.TBRaise(f'lime magic number does not match, got {magi}')
            rtype = trimNull(rtype)
            if rtype == b'scidac-private-file-xml':
                data = trimNull(f.read(size))
                dims = []
                for dim in xmlFind(data,'dims').strip().split():
                    dims.append(int(dim))
                return dims
            elif rtype == b'ildg-format':
                data = trimNull(f.read(size))
                dims = []
                for tag in ['lx','ly','lz','lt']:
                    dims.append(int(xmlFind(data,tag)))
                return dims
            else:
                f.seek(size,os.SEEK_CUR)
            f.seek(7-(size+7)%8,os.SEEK_CUR)
    logger.TBRaise('Could not find lattice dimensions in',fileName)


class confReader:

    """ 
    Base class for reading configurations. 
    """

    def __init__(self, Ns, Nt, nproc=None):
        """ 
        Initialize a confReader object.

        Args:
            Ns (int): Spatial extension of lattice. 
            Nt (int): Euclidian time extension of lattice. 
            nrows (int): Number of saved rows for compression. 
        """        
        self.Ns = Ns
        self.Nt = Nt
        self.Nd = 4
        self.Nc = 3
        checkType('int',Ns=Ns)
        checkType('int',Nt=Nt)
        if nproc==None:
            self.nproc=min(self.Nt,DEFAULTTHREADS)
        else:
            checkType('int',nproc=nproc)
            self.nproc=nproc
        self.offset = None     # How many bytes is the header?
        self.endianness = '>'  # Big endian by default
        self.precision = 'f'   # Single precision by default
        self.file = None       # File handle
        self.linkTrace = None  # <tr U>
        self.plaquette = None  # <tr U^[]>
        self.nrows = None      # number of rows
        self.checksum = None   # checksum
        self.gauge = gaugeField(self.Ns,self.Nt,self.nproc)


    def __repr__(self) -> str:
        return f"confReader(Ns={self.Ns}, Nt={self.Nt})"


    def unpack(self, data) -> SU3:
        """ 
        Unpack a string of bytes from file into an SU3 object. 
        """
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


    def getByteSize(self):
        if self.precision == 'f':
            return 4
        elif self.precision == 'd':
            return 8
        else:
            logger.TBRaise('Unknown precision',self.precision,'(expected f or d)')


    def checkLatDims(self, Nx, Ny, Nz, Nt):
        if Nx != self.Ns:
            logger.TBRaise('Read Nx =',Nx,'. Expected ',self.Ns)
        if Ny != self.Ns:
            logger.TBRaise('Read Ny =',Ny,'. Expected ',self.Ns)
        if Nz != self.Ns:
            logger.TBRaise('Read Nz =',Nz,'. Expected ',self.Ns)
        if Nt != self.Nt:
            logger.TBRaise('Read Nt =',Nt,'. Expected ',self.Nt)



class NERSCReader(confReader):

    def __repr__(self) -> str:
        return f"NERSCReader(Ns={self.Ns}, Nt={self.Nt})"

    def readHeader(self,fileName):

        """ 
        Extract metadata from the header. Do a few checks to make sure the read-in went smoothly. This method
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
            logger.TBRaise(fileName,'not in NERSC format.')

        # Convert header to metaData dictionary.
        entries = header[:self.offset-1].split(b'\n')[1:-1]
        metaData = {}
        for entry in entries:
            if b'=' in entry:
                LHS = entry.split(b'=')[0].strip().upper()
                RHS = entry.split(b'=')[1].strip()
                metaData[LHS] = RHS

        # Check that the header was read correctly.
        Nx = int(metaData[b'DIMENSION_1'])
        Ny = int(metaData[b'DIMENSION_2'])
        Nz = int(metaData[b'DIMENSION_3'])
        Nt = int(metaData[b'DIMENSION_4'])
        self.checkLatDims(Nx, Ny, Nz, Nt)

        # Extract endianness
        fp = metaData[b'FLOATING_POINT'].upper()
        if fp.endswith(b'SMALL') or fp.endswith(b'LITTLE'):
            self.endianness = '<'
        elif fp.endswith(b'BIG') or fp in (b'IEEE32', b'IEEE64'):
            # Like SIMULATeQCD, treat a missing endianness suffix as big endian.
            self.endianness = '>'
        else:
            logger.TBRaise('Unrecognized endianness.')

        # Extract precision
        if fp.startswith(b'IEEE32'):
            self.precision = 'f'
        elif fp.startswith(b'IEEE64'):
            self.precision = 'd'
        else:
            logger.TBRaise('Unrecognized precision.')

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

        # Extract checksum
        try:
            self.checksum = metaData[b'CHECKSUM'].decode('ascii').strip()
        except:
            self.checksum = None


    def readConf(self, fileName):
        """ 
        Read in the configuration. Check what you find in the binary against its metadata.

        Parameters
        ----------
        fileName : str
            Path to the NERSC configuration file.

        Returns
        -------
        gaugeField
            The loaded gauge field configuration.

        Raises
        ------
        ToolboxException
            If the header format is invalid, checksum fails, or observables do not match.
        """
        self.readHeader(fileName)

        self.file.seek(self.offset)
        payload = self.file.read()
        self.file.close()

        latprec = 8 if self.precision == 'f' else 16
        dtype_str = f'{self.endianness}c{latprec}'
        vol = (self.Ns**3) * self.Nt
        total_elements = vol * self.Nd * self.rows * 3
        raw = np.frombuffer(payload, dtype=dtype_str, count=total_elements)

        # Validate checksum against header if present
        if self.checksum is not None:
            expected_checksum = int(self.checksum, 16)
            computed_checksum = nerscChecksum(raw, self.precision)
            if expected_checksum != computed_checksum:
                logger.TBRaise(f'Checksum error: expected {expected_checksum:08x}, computed {computed_checksum:08x}')
            logger.details('Configuration', fileName, 'has correct checksum.')

        if self.rows == 3:
            self.gauge.field = raw.reshape((self.Nt, self.Ns, self.Ns, self.Ns, self.Nd, 3, 3)).astype(complex)
        elif self.rows == 2:
            subfield = raw.reshape((self.Nt, self.Ns, self.Ns, self.Ns, self.Nd, 2, 3)).astype(complex)
            r0 = subfield[..., 0, :]
            r1 = subfield[..., 1, :]
            norm0 = np.sqrt(np.sum(np.abs(r0)**2, axis=-1, keepdims=True))
            r0 = r0 / norm0
            proj = np.sum(r1 * np.conj(r0), axis=-1, keepdims=True)
            r1 = r1 - proj * r0
            norm1 = np.sqrt(np.sum(np.abs(r1)**2, axis=-1, keepdims=True))
            r1 = r1 / norm1
            r2 = np.conj(np.cross(r0, r1))
            self.gauge.field = np.stack([r0, r1, r2], axis=-2)
        else:
            logger.TBRaise(f'Unsupported number of rows: {self.rows}')

        logger.details('Loaded conf into gaugeField.')

        # Check <tr U>
        if self.linkTrace is not None:
            linkTrace = self.gauge.getLinkTrace()
            if not rel_check(linkTrace, self.linkTrace):
                logger.TBRaise('<tr U> is wrong. Compare:',linkTrace,'with',self.linkTrace)
            logger.details('Configuration',fileName,'has correct <tr U>.')

        # Check <tr U^[]>
        if self.plaquette is not None:
            plaq = self.gauge.getPlaquette()
            if not rel_check(plaq, self.plaquette):
                logger.TBRaise('<tr U^[]> is wrong. Compare:',plaq,'with',self.plaquette)
            logger.details('Configuration', fileName, 'has correct <tr U_plaq>.')

        return self.gauge



class ILDGReader(confReader):

    def __repr__(self) -> str:
        return f"ILDGReader(Ns={self.Ns}, Nt={self.Nt})"

    def readConf(self,fileName):
        """ 
        Read in the configuration. Check the checksums if they exist. Ignore spin. 
        """

        self.file       = open(fileName,'rb')
        self.nrows      = 3
        self.offset     = 144
        self.endianness = '>'

        latsize = -1
        bytesPerLink = -1
        latdatacount = -1
        latdims = []
        latnc = -1
        latprec = -1
        lattype = b'Unknown'

        while True:

            header = self.file.read(self.offset)
            if not header:
                break

            magi, vers, mbeg_end_res, size, type = struct.unpack('>ihHq128s',header)
            if magi!=LIMEMAGIC:
                logger.TBRaise(f'lime magic number does not match, got {magi}')
            type = trimNull(type)

            if type==b'scidac-binary-data' or type==b'ildg-binary-data':
                lattype = type
                latsize = size
                latdata = self.file.read(size)
                next = 7-(size+7)%8
            else:
                data = self.file.read(size)
                data = trimNull(data)
                if type==b'scidac-private-file-xml':
                    spacetime = int(xmlFind(data,'spacetime'))
                    latdims = [int(x) for x in xmlFind(data,'dims').strip().split()]
                    if len(latdims)!=spacetime:
                        logger.TBRaise(f'got spacetime {spacetime} but dims {latdims}')
                elif type==b'scidac-file-xml':
                    logger.info(f'file metadata: {data}')
                elif type==b'scidac-private-record-xml':
                    logger.info('date: ',xmlFind(data,'date'))
                    logger.info('recordtype: ',xmlFind(data,'recordtype'))
                    logger.info('datatype: ',xmlFind(data,'datatype'))
                    precision = xmlFind(data,'precision').lower()
                    if precision==b'f':
                        latprec = 8    # complex float
                    elif precision==b'd':
                        latprec = 16   # complex double
                    else:
                        logger.TBRaise(f'unknown precision {precision}')
                    latnc = int(xmlFind(data,'colors'))
                    bytesPerLink = int(xmlFind(data,'typesize'))
                    latdatacount = int(xmlFind(data,'datacount'))
                elif type==b'scidac-record-xml':
                    logger.info(f'record metadata: {data}')
                elif type==b'scidac-checksum':
                    latsuma = int(xmlFind(data,'suma'),16)
                    latsumb = int(xmlFind(data,'sumb'),16)
                    pass
                elif type==b'ildg-format':
                    field = xmlFind(data,'field')
                    if field!=b'su3gauge':
                        logger.TBRaise(f'unsupported ildg field type {field}')
                    precision = int(xmlFind(data,'precision'))
                    if precision==32:
                        latprec = 8    # complex float
                    elif precision==64:
                        latprec = 16   # complex double
                    else:
                        logger.TBRaise(f'unknown precision {precision}')
                    lx = int(xmlFind(data,'lx'))
                    ly = int(xmlFind(data,'ly'))
                    lz = int(xmlFind(data,'lz'))
                    lt = int(xmlFind(data,'lt'))
                    latdims = [lx,ly,lz,lt]
                    latnc = self.Nc
                    bytesPerLink = latnc*latnc*latprec
                    latdatacount = self.Nd
                else:
                    logger.info(f'unused type: {type}  data: {data}')
                next = 7-(size+7)%8
            self.file.seek(next,os.SEEK_CUR)
        self.file.close()
        ndim = len(latdims)

        # Carry out some checks on the read-in data.
        if latsize<=0 or bytesPerLink<=0 or latdatacount<=0 or ndim==0 or latnc<0 or latprec<0 or lattype==b'Unknown':
            logger.TBRaise(f'unsupported file: {fileName}')
        if latprec==8:
            self.precision = 'f'
        elif latprec==16:
            self.precision = 'd'
        else:
            logger.TBRaise(f'unknown complex precision {latprec}')
        if latnc != self.Nc:
            logger.TBRaise(f'unsupported number of colors: {latnc}')
        if latdatacount != ndim:
            logger.TBRaise(f'There should be one link per direction, but header suggests {latdatacount} links.')
        self.checkLatDims(*latdims)
        if ndim!=self.Nd:
            logger.TBRaise(f'unsupported number of dimensions: {ndim}.')
        vol = self.Ns**3 * self.Nt 
        if latsize != vol*self.Nd*bytesPerLink:
            logger.TBRaise(f'incorrect lattice size, expect {vol*ndim*bytesPerLink}, but got {latsize}')
        if bytesPerLink != latprec*latnc*latnc:
            logger.TBRaise(f'incorrect bytes per link, expected {latprec*latnc*latnc}, but got {bytesPerLink}')
        if not (lattype==b'scidac-binary-data' or lattype==b'ildg-binary-data'):
            logger.TBRaise(f'unknown lattice format: {lattype}')

        # Check checksums
        if lattype==b'scidac-binary-data':
            if latsuma is None or latsumb is None:
                logger.TBRaise(f'No SciDAC checksum.')
            suma,sumb = scidacChecksum(latdata, vol, ndim*bytesPerLink)
            if suma!=latsuma or sumb!=latsumb:
                logger.TBRaise(f'Checksum error: expected {latsuma:x} {latsumb:x}, computed {suma:x} {sumb:x}')
            logger.info('SciDAC checksum OK.')

        # Convert from bytes to numpy array.
        lat = np.frombuffer(latdata,dtype='>c'+str(latprec),count=vol*ndim*latnc*latnc)
        lat = np.reshape(lat,latdims[::-1]+[ndim,latnc,latnc])
        
        # Load into gaugeField.
        for t in range(self.Nt):
            for z in range(self.Ns):
                for y in range(self.Ns):
                    for x in range(self.Ns):
                        for mu in range(self.Nd):
                            link = lat[t,z,y,x,mu]
                            self.gauge.setLink(link,x,y,z,t,mu)
        logger.details('Loaded conf into gaugeField.')

        return self.gauge 



class DynisoReader(confReader):

    """
    Read the configurations written by dyniso's save_lattice. The format is detected from the file:
        comp:  4 ints nx, ny, nz, nt. Then each link U is stored as the 8 floats a_k = Im tr(log(U) lambda_k), with
               lambda_k the Gell-Mann matrices, so that U = exp(i/2 sum_k a_k lambda_k). dyniso's checkpoints
               are in this format.
        milc:  MILC header (magic number 20103, nx, ny, nz, nt, 64-character date, order, two checksums,
               which dyniso leaves at 0), then each link as a 3x3 complex float matrix.
        ascii: version number 59354, nx ny nz nt, two couplings, then each element of each link as a line
               "re im". The files read by dyniso's reload_ascii start directly with nx ny nz nt.
    The links run over t, z, y, x, then direction, in all three formats. The binary formats are in the byte order
    of the machine that wrote them, which is detected from the header. Like dyniso, we reunitarize the links
    after reading them.
    """

    def readConf(self,fileName):
        """
        Read in the configuration.
        """
        nlinks = self.Nd*self.Ns**3*self.Nt
        fileSize = os.path.getsize(fileName)
        if fileSize == 16 + 32*nlinks:
            links = self._readComp(fileName)
        elif fileSize == 96 + 72*nlinks:
            links = self._readMILC(fileName)
        else:
            links = self._readASCII(fileName,nlinks)
        self.gauge.field[:] = np.reshape(links,(self.Nt,self.Ns,self.Ns,self.Ns,self.Nd,self.Nc,self.Nc))
        for t in range(self.Nt):
            for z in range(self.Ns):
                for y in range(self.Ns):
                    for x in range(self.Ns):
                        for mu in range(self.Nd):
                            fastUnitarize(self.gauge.field[t,z,y,x,mu])
        logger.details('Loaded conf into gaugeField.')
        return self.gauge


    def _endianness(self, header, nskip):
        """
        Find the byte order for which the 4 ints after the first nskip ints of header are the lattice dimensions.
        """
        for endianness in ['<','>']:
            dims = struct.unpack(endianness+str(nskip+4)+'i',header[:4*(nskip+4)])[nskip:]
            if dims == (self.Ns,self.Ns,self.Ns,self.Nt):
                return endianness
        logger.TBRaise('Could not find lattice dimensions',self.Ns,self.Ns,self.Ns,self.Nt,'in the header.')


    def _readComp(self, fileName):
        with open(fileName,'rb') as f:
            data = f.read()
        endianness = self._endianness(data,0)
        alpha = np.frombuffer(data,dtype=endianness+'f4',offset=16).astype(float)
        alpha = np.reshape(alpha,(-1,8))
        # U = exp(iM) with Hermitian M = sum_k a_k lambda_k / 2. Exponentiate in the eigenbasis of M, as dyniso does.
        M = np.zeros((len(alpha),self.Nc,self.Nc),dtype=complex)
        for k in range(1,9):
            M += alpha[:,k-1,None,None]*lam[k]/2
        w, V = np.linalg.eigh(M)
        return (V*np.exp(1j*w)[:,None,:]) @ np.conj(np.swapaxes(V,1,2))


    def _readMILC(self, fileName):
        with open(fileName,'rb') as f:
            data = f.read()
        endianness = self._endianness(data,1)
        magic = struct.unpack(endianness+'i',data[:4])[0]
        if magic != DYNISOMILCMAGIC:
            logger.TBRaise('Expected MILC magic number',DYNISOMILCMAGIC,'but got',magic)
        return np.frombuffer(data,dtype=endianness+'c8',offset=96).astype(complex)


    def _readASCII(self, fileName, nlinks):
        try:
            with open(fileName,'r') as f:
                tokens = f.read().split()
        except UnicodeDecodeError:
            logger.TBRaise(fileName,'is not a dyniso configuration of size',self.Ns,self.Ns,self.Ns,self.Nt)
        nvalues = 2*self.Nc*self.Nc*nlinks
        if len(tokens) == 7 + nvalues:
            if int(tokens[0]) != DYNISOASCIIMAGIC:
                logger.TBRaise('Expected version number',DYNISOASCIIMAGIC,'but got',tokens[0])
            dims = tokens[1:5]
        elif len(tokens) == 4 + nvalues:
            dims = tokens[0:4]
        else:
            logger.TBRaise(fileName,'is not a dyniso configuration of size',self.Ns,self.Ns,self.Ns,self.Nt)
        self.checkLatDims(int(dims[0]),int(dims[1]),int(dims[2]),int(dims[3]))
        values = np.array(tokens[len(tokens)-nvalues:],dtype=float)
        return values[0::2] + 1j*values[1::2]
