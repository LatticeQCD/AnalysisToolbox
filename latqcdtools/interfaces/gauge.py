# 
# gauge.py                                                               
# 
# D. Clarke
# 
# Tools for reading and writing gauge configurations in Python. To interact with binary configurations, we can use
# Python's built-in struct module. More info on that: https://docs.python.org/3.7/library/struct.html.
#
# Read/write code was created by closely following QCDUtils https://github.com/mdipierro/qcdutils as an example.
# Unfortunately that was written for Python 2. Anyway it is useful to have such methods incorporated into the toolbox.
#


import struct
import latqcdtools.base.logger as logger


def reconstruct(link):
    """ Reconstruct a unitary 3x3 matrix from its first two rows. """
    (a1re, a1im, a2re, a2im, a3re, a3im, b1re, b1im, b2re, b2im, b3re, b3im) = link
    c1re =   a2re * b3re - a2im * b3im - a3re * b2re + a3im * b2im
    c1im = -(a2re * b3im + a2im * b3re - a3re * b2im - a3im * b2re)
    c2re =   a3re * b1re - a3im * b1im - a1re * b3re + a1im * b3im
    c2im = -(a3re * b1im + a3im * b1re - a1re * b3im - a1im * b3re)
    c3re =   a1re * b2re - a1im * b2im - a2re * b1re + a2im * b1im
    c3im = -(a1re * b2im + a1im * b2re - a2re * b1im - a2im * b1re)
    return (a1re, a1im, a2re, a2im, a3re, a3im,
            b1re, b1im, b2re, b2im, b3re, b3im,
            c1re, c1im, c2re, c2im, c3re, c3im)


class gaugeField:

    """ Base class for gauge fields. """

    def __init__(self, Ns=None, Nt=None):
        self.Ns = Ns           # Spatial extension
        self.Nt = Nt           # Euclidean time extension
        self.offset = None     # How many bytes is the header?
        self.endianness = '>'  # Big endian by default
        self.precision = 'f'   # Single precision by default
        self.file = None       # File handle
        if (Ns is None) or (Nt is None):
           logger.TBError("gaugeField objects must be initialized with a size!")
        logger.info("Initialized "+str(self.Ns)+"^3x"+str(self.Nt)+" gaugeField object.")

    def unpack(self, data):
        """ Unpack a string of bytes from file into a list of float/double. """
        numData = int(len(data)/self.getByteSize())
        items = struct.unpack(self.endianness + str(numData) + self.precision, data)
        return items

    def pack(self, items):
        """ Packs a list of float/double into a string of bytes. """
        numData = len(items)
        return struct.pack(self.endianness + str(numData) + self.precision, *items)

    def getByteSize(self):
        if self.precision == 'f':
            return 4
        elif self.precision == 'd':
            return 8
        else:
            logger.TBError('Unknown precision',self.precision,'(expected f or d)')


class fieldNERSC(gaugeField):

    def readHeader(self,fileName):

        """ Extract offset, endianness, and precision from the header. Do a few checks to make sure the read-in went
        smoothly. This method assumes all NERSC headers have the minimal entries
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


#    template<class f1, class f2>
#    void to_buf(f1 *buf, const GSU3<f2> &U) const {
#        int i = 0;
#        for (int j = 0; j < rows; j++)
#            for (int k = 0; k < 3; k++) {
#                buf[i++] = U(j, k).cREAL;
#                buf[i++] = U(j, k).cIMAG;
#            }
#    }

#    template<class floatT>
#    void put(GSU3<floatT> U) {
#        char *start = &buf[index];
#        if (float_size == 4)
#            to_buf((float *) start, U);
#        else if (float_size == 8)
#            to_buf((double *) start, U);
#        index += su3_size;
#    }

    def readConf(self,fileName):

        self.readHeader(fileName)

        bytesPerLink = 6*2*self.getByteSize()

        for t in range(self.Nt):
            for z in range(self.Ns):
                for y in range(self.Ns):
                    for x in range(self.Ns):
                        for mu in range(4):
                            byteIndex = self.offset + bytesPerLink*( mu + 4*( x + self.Ns*( y + self.Ns*( z + self.Ns*t ) ) ) )
                            self.file.seek(byteIndex)
                            data = self.file.read(bytesPerLink)
                            link = self.unpack(data)

                        exit()
#        if self.reunitarize:
#            new_items = []
#            for i in range(4):
#                new_items += reunitarize(items[i*12:(i+1)*12])
#            items = new_items
#        return items