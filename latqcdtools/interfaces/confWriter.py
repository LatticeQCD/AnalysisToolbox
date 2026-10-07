# 
# confWriter.py                                                               
# 
# K. Ebira, D. Clarke
# 
# Tools for writing gauge configurations in Python.
#

import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType
from latqcdtools.base.fileSystem import createFilePath
from latqcdtools.physics.gauge import gaugeField


class confWriter:
    """
    Base class for writing gauge field configurations.

    Parameters
    ----------
    gauge : gaugeField
        SU(3) gauge field to be written.
    """

    def __init__(self, gauge):
        checkType(gaugeField, gauge=gauge)
        self.gauge = gauge
        self.Ns = gauge.Ns
        self.Nt = gauge.Nt
        self.Nd = 4
        self.Nc = 3

    def __repr__(self) -> str:
        return f"confWriter(Ns={self.Ns}, Nt={self.Nt})"


class NERSCWriter(confWriter):
    """
    Writer for NERSC format gauge configurations.

    Parameters
    ----------
    gauge : gaugeField
        SU(3) gauge field to be written.
    """

    def __repr__(self) -> str:
        return f"NERSCWriter(Ns={self.Ns}, Nt={self.Nt})"

    def writeConf(self, fileName, format='3x3', endianness='>', precision='d'):
        """
        Write gauge configuration to file in NERSC format.

        Parameters
        ----------
        fileName : str
            Path to the output file.
        format : str, optional
            '3x3' to write all 3 rows per link, or '2x3' to write only the first 2 rows.
            Default is '3x3'.
        endianness : str, optional
            Byte ordering: '>' for big-endian (default) or '<' for little-endian.
        precision : str, optional
            Floating point precision: 'd' for IEEE64 (default) or 'f' for IEEE32.

        Raises
        ------
        ToolboxException
            If format, endianness, or precision is invalid.
        """
        checkType(str, fileName=fileName)
        createFilePath(fileName)

        if str(format) in ('3x3', '3'):
            datatype = '4D_SU3_GAUGE_3x3'
            subfield = self.gauge.field
        elif str(format) in ('2x3', '2'):
            datatype = '4D_SU3_GAUGE'
            subfield = self.gauge.field[..., :2, :]
        else:
            logger.TBRaise(f"Unsupported format '{format}'. Expected '3x3' or '2x3'.")

        end_u = str(endianness).upper()
        if end_u in ('>', 'BIG'):
            end_code = '>'
            end_str = 'BIG'
        elif end_u in ('<', 'LITTLE', 'SMALL'):
            end_code = '<'
            end_str = 'LITTLE'
        else:
            logger.TBRaise(f"Unsupported endianness '{endianness}'.")

        prec_u = str(precision).upper()
        if prec_u in ('D', '64', 'IEEE64', 'DOUBLE'):
            c_dtype = 'c16'
            prec_str = 'IEEE64'
        elif prec_u in ('F', '32', 'IEEE32', 'SINGLE'):
            c_dtype = 'c8'
            prec_str = 'IEEE32'
        else:
            logger.TBRaise(f"Unsupported precision '{precision}'.")

        payload = subfield.astype(f'{end_code}{c_dtype}').tobytes()

        # NERSC 32-bit checksum: sum of big-endian unsigned 32-bit words
        chk_arr = np.frombuffer(payload, dtype='>u4')
        checksum = int(np.sum(chk_arr, dtype=np.uint64) & 0xffffffff)

        link_trace = self.gauge.getLinkTrace()
        plaquette = self.gauge.getPlaquette()

        header = (
            f"BEGIN_HEADER\n"
            f"DATATYPE = {datatype}\n"
            f"DIMENSION_1 = {self.Ns}\n"
            f"DIMENSION_2 = {self.Ns}\n"
            f"DIMENSION_3 = {self.Ns}\n"
            f"DIMENSION_4 = {self.Nt}\n"
            f"CHECKSUM = {checksum:08x}\n"
            f"LINK_TRACE = {link_trace:.10f}\n"
            f"PLAQUETTE = {plaquette:.10f}\n"
            f"FLOATING_POINT = {prec_str}{end_str}\n"
            f"END_HEADER\n"
        ).encode('ascii')

        with open(fileName, 'wb') as f:
            f.write(header)
            f.write(payload)

        logger.details(f"Wrote NERSC configuration to {fileName}.")


def writeConf(fileName, gauge, format='3x3', endianness='>', precision='d'):
    """
    Convenience function to write a gaugeField configuration to file in NERSC format.

    Parameters
    ----------
    fileName : str
        Path to the output file.
    gauge : gaugeField
        The SU(3) gauge field configuration to write.
    format : str, optional
        '3x3' to write full 3 rows per link, or '2x3' to write first 2 rows.
        Default is '3x3'.
    endianness : str, optional
        Byte order: '>' for big-endian (default) or '<' for little-endian.
    precision : str, optional
        'd' for IEEE64 double precision (default) or 'f' for IEEE32 single precision.

    Raises
    ------
    ToolboxException
        If format, endianness, or precision is invalid, or if fileName/gauge types are invalid.
    """
    if not isinstance(fileName, str) and isinstance(gauge, str):
        fileName, gauge = gauge, fileName
    writer = NERSCWriter(gauge)
    writer.writeConf(fileName, format=format, endianness=endianness, precision=precision)
