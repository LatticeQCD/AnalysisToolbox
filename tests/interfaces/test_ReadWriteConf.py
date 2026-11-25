# 
# test_ReadWriteConf.py                                                               
# 
# D. Clarke, K. Ebira, Google Antigravity
# 
# Testing reading and writing gauge configurations.
# 

import os
import re
import tempfile
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.logger import ToolboxException
from latqcdtools.interfaces.confReader import NERSCReader, ILDGReader, DynisoReader, detectConf
from latqcdtools.interfaces.confWriter import NERSCWriter, writeConf
from latqcdtools.base.utilities import timer
from latqcdtools.math.math import rel_check
from latqcdtools.testing import concludeTest


def _datasetPath(filename: str) -> str:
    """Resolve a file in the repository datasets/ directory regardless of working directory."""
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
    return os.path.join(repo_root, 'datasets', filename)


def testReadWriteConf():

    timing = timer()

    path_nersc = _datasetPath('nersc.l8t4b3360')
    path_ildg = _datasetPath('ildg.l8t4b3360')

    reader = NERSCReader(Ns=8, Nt=4)
    gauge1 = reader.readConf(path_nersc)
    timing.printTiming()
    gauge1.checkSU3()
    timing.printTiming()

    lpass = True

    if reader.checksum is None or reader.checksum.lower() != 'b379560a':
        logger.TBFail(f'Initial NERSC checksum mismatch: expected b379560a, got {reader.checksum}')
        lpass = False

    reader = ILDGReader(Ns=8, Nt=4)
    gauge2 = reader.readConf(path_ildg)
    timing.printTiming()
    gauge2.checkSU3()
    timing.printTiming()

    for mu in range(4):
        for t in range(4):
            for z in range(8):
                for y in range(8):
                    for x in range(8):
                        link1 = gauge1.getLink(x,y,z,t,mu)
                        link2 = gauge2.getLink(x,y,z,t,mu)
                        lpass *= rel_check(link1, link2)

    # The dyniso files were written by dyniso's save_lattice code from the NERSC conf, in single precision. The
    # ascii format keeps 6 significant digits.
    for fmt, prec in [('comp',1e-6), ('milc',1e-6), ('ascii',1e-5)]:
        reader = DynisoReader(Ns=8, Nt=4)
        gauge3 = reader.readConf(_datasetPath('dyniso_' + fmt + '.l8t4b3360'))
        timing.printTiming()
        gauge3.checkSU3()
        diff = np.max(np.abs(gauge3.field-gauge1.field))
        if diff > prec:
            logger.TBFail('dyniso',fmt,'differs from NERSC conf by',diff)
            lpass = False
        else:
            logger.TBPass('dyniso',fmt,'agrees with NERSC conf.')

    # The format and lattice size should be detected from each file.
    for fileName, fmt in [(path_nersc,'NERSC'), (path_ildg,'ILDG'),
                          (_datasetPath('dyniso_comp.l8t4b3360'),'dyniso'),
                          (_datasetPath('dyniso_milc.l8t4b3360'),'dyniso'),
                          (_datasetPath('dyniso_ascii.l8t4b3360'),'dyniso')]:
        if detectConf(fileName) != (fmt, 8, 4):
            logger.TBFail('detectConf gave',detectConf(fileName),'for',fileName)
            lpass = False

    # Temporary test files for writeConf testing isolated inside a TemporaryDirectory
    with tempfile.TemporaryDirectory() as tmpdir:
        file_3x3 = os.path.join(tmpdir, 'temp_test_3x3.nersc')
        file_2x3 = os.path.join(tmpdir, 'temp_test_2x3.nersc')
        file_2x3_le = os.path.join(tmpdir, 'temp_test_2x3_le.nersc')
        file_le = os.path.join(tmpdir, 'temp_test_le.nersc')
        file_f = os.path.join(tmpdir, 'temp_test_f.nersc')
        file_f_le = os.path.join(tmpdir, 'temp_test_f_le.nersc')
        file_swap = os.path.join(tmpdir, 'temp_test_swap.nersc')
        file_corrupt = os.path.join(tmpdir, 'temp_test_corrupt.nersc')
        file_corrupt_le = os.path.join(tmpdir, 'temp_test_corrupt_le.nersc')
        file_trunc = os.path.join(tmpdir, 'temp_test_trunc.nersc')
        file_trunc_le = os.path.join(tmpdir, 'temp_test_trunc_le.nersc')

        # Test 1: Round-trip 3x3 write and read
        writer = NERSCWriter(gauge1)
        writer_repr = repr(writer)
        if 'NERSCWriter' not in writer_repr:
            logger.TBFail('NERSCWriter __repr__ failed')
            lpass = False

        writer.writeConf(file_3x3, format='3x3')
        reader_3x3 = NERSCReader(Ns=8, Nt=4)
        gauge_3x3 = reader_3x3.readConf(file_3x3)
        diff_3x3 = np.max(np.abs(gauge1.field - gauge_3x3.field))
        if diff_3x3 > 1e-12:
            logger.TBFail(f'3x3 roundtrip diff too large: {diff_3x3}')
            lpass = False
        else:
            logger.details('3x3 roundtrip test passed.')

        # Test 2: Round-trip 2x3 write and read using convenience function writeConf
        writeConf(file_2x3, gauge1, format='2x3')
        reader_2x3 = NERSCReader(Ns=8, Nt=4)
        gauge_2x3 = reader_2x3.readConf(file_2x3)
        diff_2x3 = np.max(np.abs(gauge1.field - gauge_2x3.field))
        if diff_2x3 > 1e-6:
            logger.TBFail(f'2x3 roundtrip diff too large: {diff_2x3}')
            lpass = False
        else:
            logger.details('2x3 roundtrip test passed.')

        # Also verify 2x3 format with little-endian and checksum invariance
        writeConf(file_2x3_le, gauge1, format='2x3', endianness='<')
        reader_2x3_le = NERSCReader(Ns=8, Nt=4)
        gauge_2x3_le = reader_2x3_le.readConf(file_2x3_le)
        diff_2x3_le = np.max(np.abs(gauge1.field - gauge_2x3_le.field))
        if diff_2x3_le > 1e-6:
            logger.TBFail(f'2x3 little-endian roundtrip diff too large: {diff_2x3_le}')
            lpass = False
        elif reader_2x3_le.checksum.lower() != reader_2x3.checksum.lower():
            logger.TBFail(f'2x3 checksum endian mismatch: BE={reader_2x3.checksum}, LE={reader_2x3_le.checksum}')
            lpass = False
        else:
            logger.details('2x3 little-endian roundtrip and checksum test passed.')

        # Test 3: Round-trip little-endian configuration
        writeConf(file_le, gauge1, endianness='<')
        reader_le = NERSCReader(Ns=8, Nt=4)
        gauge_le = reader_le.readConf(file_le)
        diff_le = np.max(np.abs(gauge1.field - gauge_le.field))
        if diff_le > 1e-12:
            logger.TBFail(f'Little-endian roundtrip diff too large: {diff_le}')
            lpass = False
        else:
            logger.details('Little-endian roundtrip test passed.')

        if reader_le.checksum.lower() != 'b379560a':
            logger.TBFail(f'Little-endian checksum mismatch: expected b379560a, got {reader_le.checksum}')
            lpass = False
        else:
            logger.details('Little-endian checksum preservation test passed.')

        # Test 4: Round-trip single precision configuration
        writeConf(file_f, gauge1, precision='f')
        reader_f = NERSCReader(Ns=8, Nt=4)
        gauge_f = reader_f.readConf(file_f)
        diff_f = np.max(np.abs(gauge1.field - gauge_f.field))
        if diff_f > 1e-6:
            logger.TBFail(f'Single precision roundtrip diff too large: {diff_f}')
            lpass = False
        else:
            logger.details('Single precision roundtrip test passed.')

        # Also verify single precision with little-endian and checksum invariance
        writeConf(file_f_le, gauge1, precision='f', endianness='<')
        reader_f_le = NERSCReader(Ns=8, Nt=4)
        gauge_f_le = reader_f_le.readConf(file_f_le)
        diff_f_le = np.max(np.abs(gauge1.field - gauge_f_le.field))
        if diff_f_le > 1e-6:
            logger.TBFail(f'Single precision little-endian roundtrip diff too large: {diff_f_le}')
            lpass = False
        elif reader_f_le.checksum.lower() != reader_f.checksum.lower():
            logger.TBFail(f'Single precision checksum endian mismatch: BE={reader_f.checksum}, LE={reader_f_le.checksum}')
            lpass = False
        else:
            logger.details('Single precision little-endian roundtrip and checksum test passed.')

        # The checksum is computed from the link values, so it must not depend on endianness. For the
        # original conf it must reproduce the checksum in that file's header.
        reader_orig = NERSCReader(Ns=8, Nt=4)
        reader_orig.readHeader(path_nersc)
        reader_orig.file.close()
        if reader_le.checksum != reader_3x3.checksum or int(reader_3x3.checksum,16) != int(reader_orig.checksum,16):
            logger.TBFail('Checksum depends on endianness:',reader_3x3.checksum,reader_le.checksum,reader_orig.checksum)
            lpass = False
        else:
            logger.details('Endianness-independent checksum test passed.')

        # Test 5: Swapped argument order for writeConf
        writeConf(gauge1, file_swap)
        reader_swap = NERSCReader(Ns=8, Nt=4)
        gauge_swap = reader_swap.readConf(file_swap)
        diff_swap = np.max(np.abs(gauge1.field - gauge_swap.field))
        if diff_swap > 1e-12:
            logger.TBFail(f'Swapped argument order roundtrip diff too large: {diff_swap}')
            lpass = False
        else:
            logger.details('Swapped argument order test passed.')

        # Test 6: Checksum failure detection without hardcoded magic checksum
        with open(file_3x3, 'rb') as f:
            content = f.read()
        corrupted_content = re.sub(rb'CHECKSUM = [0-9a-fA-F]+', b'CHECKSUM = 00000000', content)
        with open(file_corrupt, 'wb') as f:
            f.write(corrupted_content)

        corrupt_raised = False
        try:
            reader_bad = NERSCReader(Ns=8, Nt=4)
            reader_bad.readConf(file_corrupt)
        except ToolboxException:
            corrupt_raised = True

        if not corrupt_raised:
            logger.TBFail('NERSCReader failed to detect checksum mismatch')
            lpass = False
        else:
            logger.details('Checksum failure test passed.')

        # Also verify checksum failure detection on little-endian file
        with open(file_le, 'rb') as f:
            content_le = f.read()
        corrupted_content_le = re.sub(rb'CHECKSUM = [0-9a-fA-F]+', b'CHECKSUM = 00000000', content_le)
        with open(file_corrupt_le, 'wb') as f:
            f.write(corrupted_content_le)

        corrupt_le_raised = False
        try:
            reader_bad_le = NERSCReader(Ns=8, Nt=4)
            reader_bad_le.readConf(file_corrupt_le)
        except ToolboxException:
            corrupt_le_raised = True

        if not corrupt_le_raised:
            logger.TBFail('NERSCReader failed to detect little-endian checksum mismatch')
            lpass = False
        else:
            logger.details('Little-endian checksum failure test passed.')

        # Test 7: Truncated payload detection
        with open(file_trunc, 'wb') as f:
            f.write(content[:-3])
        trunc_raised = False
        try:
            reader_trunc = NERSCReader(Ns=8, Nt=4)
            reader_trunc.readConf(file_trunc)
        except ToolboxException:
            trunc_raised = True

        if not trunc_raised:
            logger.TBFail('NERSCReader failed to detect truncated file')
            lpass = False
        else:
            logger.details('Payload truncation test passed.')

        # Also verify truncation detection on little-endian file
        with open(file_trunc_le, 'wb') as f:
            f.write(content_le[:-3])
        trunc_le_raised = False
        try:
            reader_trunc_le = NERSCReader(Ns=8, Nt=4)
            reader_trunc_le.readConf(file_trunc_le)
        except ToolboxException:
            trunc_le_raised = True

        if not trunc_le_raised:
            logger.TBFail('NERSCReader failed to detect truncated little-endian file')
            lpass = False
        else:
            logger.details('Little-endian payload truncation test passed.')

        # Verify 4-byte-aligned truncation raises, while trailing padding bytes remain tolerated for backwards compatibility
        with open(file_trunc, 'wb') as f:
            f.write(content[:-16])
        aligned_raised = False
        try:
            reader_aligned = NERSCReader(Ns=8, Nt=4)
            reader_aligned.readConf(file_trunc)
        except ToolboxException:
            aligned_raised = True
        if not aligned_raised:
            logger.TBFail('NERSCReader failed to detect aligned truncation')
            lpass = False
        else:
            logger.details('Payload aligned truncation test passed.')

        with open(file_trunc, 'wb') as f:
            f.write(content + b'\x00\x00\x00\x00')
        try:
            reader_padded = NERSCReader(Ns=8, Nt=4)
            g_padded = reader_padded.readConf(file_trunc)
            diff_padded = np.max(np.abs(g_padded.field - gauge1.field))
            if diff_padded > 1e-14:
                logger.TBFail(f'NERSC trailing padding link mismatch: max diff {diff_padded}')
                lpass = False
            else:
                logger.details('NERSC trailing padding compat test passed.')
        except Exception as e:
            logger.TBFail(f'NERSCReader rejected trailing padding bytes: {e}')
            lpass = False

        # Verify uppercase DATATYPE 3X3 is handled and non-hex CHECKSUM raises ToolboxException
        upper_dt = re.sub(rb'4D_SU3_GAUGE_3x3', b'4D_SU3_GAUGE_3X3', content, count=1)
        with open(file_trunc, 'wb') as f:
            f.write(upper_dt)
        try:
            reader_upper = NERSCReader(Ns=8, Nt=4)
            g_upper = reader_upper.readConf(file_trunc)
            if np.max(np.abs(g_upper.field - gauge1.field)) > 1e-14:
                logger.TBFail('NERSCReader failed to read uppercase 4D_SU3_GAUGE_3X3')
                lpass = False
        except Exception as e:
            logger.TBFail(f'NERSCReader crashed on uppercase 4D_SU3_GAUGE_3X3: {e}')
            lpass = False

        bad_hex = re.sub(rb'CHECKSUM = [0-9a-fA-F]+', b'CHECKSUM = ZZZZZZZZ', content)
        with open(file_corrupt, 'wb') as f:
            f.write(bad_hex)
        bad_hex_raised = False
        try:
            reader_bad_hex = NERSCReader(Ns=8, Nt=4)
            reader_bad_hex.readConf(file_corrupt)
        except ToolboxException:
            bad_hex_raised = True
        if not bad_hex_raised:
            logger.TBFail('NERSCReader failed to raise ToolboxException on non-hex CHECKSUM')
            lpass = False

    concludeTest(lpass) 

if __name__ == '__main__':
    testReadWriteConf()