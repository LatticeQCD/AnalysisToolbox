# 
# test_ReadWriteConf.py                                                               
# 
# D. Clarke, K. Ebira
# 
# Testing reading and writing gauge configurations.
# 

import os
import numpy as np
import latqcdtools.base.logger as logger
from latqcdtools.base.logger import ToolboxException
from latqcdtools.interfaces.confReader import NERSCReader, ILDGReader
from latqcdtools.interfaces.confWriter import NERSCWriter, writeConf
from latqcdtools.base.utilities import timer
from latqcdtools.math.math import rel_check
from latqcdtools.testing import concludeTest


def testReadWriteConf():

    timing = timer()

    path_nersc = '../../datasets/nersc.l8t4b3360'
    if not os.path.exists(path_nersc) and os.path.exists('datasets/nersc.l8t4b3360'):
        path_nersc = 'datasets/nersc.l8t4b3360'

    path_ildg = '../../datasets/ildg.l8t4b3360'
    if not os.path.exists(path_ildg) and os.path.exists('datasets/ildg.l8t4b3360'):
        path_ildg = 'datasets/ildg.l8t4b3360'

    reader = NERSCReader(Ns=8, Nt=4)
    gauge1 = reader.readConf(path_nersc)
    timing.printTiming()
    gauge1.checkSU3()
    timing.printTiming()

    reader = ILDGReader(Ns=8, Nt=4)
    gauge2 = reader.readConf(path_ildg)
    timing.printTiming()
    gauge2.checkSU3()
    timing.printTiming()

    lpass = True
    for mu in range(4):
        for t in range(4):
            for z in range(8):
                for y in range(8):
                    for x in range(8):
                        link1 = gauge1.getLink(x,y,z,t,mu)
                        link2 = gauge2.getLink(x,y,z,t,mu)
                        lpass *= rel_check(link1, link2)

    # Temporary test files for writeConf testing
    file_3x3 = 'temp_test_3x3.nersc'
    file_2x3 = 'temp_test_2x3.nersc'
    file_le = 'temp_test_le.nersc'
    file_f = 'temp_test_f.nersc'
    file_swap = 'temp_test_swap.nersc'
    file_corrupt = 'temp_test_corrupt.nersc'
    file_trunc = 'temp_test_trunc.nersc'
    all_temp_files = [file_3x3, file_2x3, file_le, file_f, file_swap, file_corrupt, file_trunc]

    try:
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
        import re
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

        # Test 7: Truncated payload detection
        with open(file_trunc, 'wb') as f:
            f.write(content[:-3])
        trunc_raised = False
        try:
            reader_trunc = NERSCReader(Ns=8, Nt=4)
            reader_trunc.readConf(file_trunc)
        except (ToolboxException, ValueError):
            trunc_raised = True

        if not trunc_raised:
            logger.TBFail('NERSCReader failed to detect truncated file')
            lpass = False
        else:
            logger.details('Payload truncation test passed.')

    finally:
        for f in all_temp_files:
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass

    concludeTest(lpass) 

if __name__ == '__main__':
    testReadWriteConf()