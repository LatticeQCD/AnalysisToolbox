# 
# test_Lime.py                                                               
# 
# D. Clarke, K. Ebira
# 
# Testing the Lime interface. 
# 
from latqcdtools.interfaces.lime import limeHeader, trimNull
from latqcdtools.testing import print_results,concludeTest
import latqcdtools.base.logger as logger


def testLime():

    lpass = True

    testByteString = limeHeader(1,0,123,b'test') 

    lpass *= print_results( len(testByteString), 144, text="header length")

    if trimNull(testByteString) !=  b'Eg\x89\xab':
        lpass = False
        logger.TBFail('null trim')

    import struct
    # Issue #18: Verify mbeg and mend flags are preserved when both are True
    hdr_both = limeHeader(1, 1, 100, b'both')
    m_both = struct.unpack('>ihHq128s', hdr_both)[2]
    expected_both = (1 << 15) | (1 << 14)
    if m_both != expected_both:
        logger.TBFail(f'limeHeader(1, 1) m-flag: got {m_both}, expected {expected_both}')
        lpass = False

    hdr_end_only = limeHeader(0, 1, 100, b'end')
    m_end_only = struct.unpack('>ihHq128s', hdr_end_only)[2]
    expected_end_only = 1 << 14
    if m_end_only != expected_end_only:
        logger.TBFail(f'limeHeader(0, 1) m-flag: got {m_end_only}, expected {expected_end_only}')
        lpass = False

    hdr_none = limeHeader(0, 0, 100, b'none')
    m_none = struct.unpack('>ihHq128s', hdr_none)[2]
    if m_none != 0:
        logger.TBFail(f'limeHeader(0, 0) m-flag: got {m_none}, expected 0')
        lpass = False

    concludeTest(lpass)


if __name__ == '__main__':
    testLime()
