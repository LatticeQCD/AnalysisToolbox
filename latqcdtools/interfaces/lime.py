# 
# lime.py                                                               
# 
# D. Clarke, X.-Y. Jin
# 
# Some tools for working with the LIME format. This is a port of X.-Y. Jin's original Python code,
# which can be found here: https://github.com/nftqcd/nthmc/blob/master/lib/fieldio.py. Information
# about the LIME format can be found here:
#
# To help with understanding LIME format, note the following organizational hierarchy:
#     message > record > word
# Message begin (MB) and message end (ME) flags have to be set according to the following rules:
#
# (M1)  MB=1 for the first record of the LIME file.
# (M2)  For any two consecutive records with ME flag me1 and MB flag mb2, respectively, the relation me1=mb2 must hold.
# (M3)  ME=1 for the last record of the LIME file.
#
# Special thanks to H. Simma for the explanation. Since we are not using the message organizational level, we
# will always set these flags to 1.
# 

import os, re, struct, zlib
import latqcdtools.base.logger as logger
from latqcdtools.base.check import checkType

LIMEMAGIC = 0x456789ab


def limeHeader(mbeg, mend, size, type) -> bytes:
    """
    Get the header for a LIME record.

    Args:
        mbeg (bool): start flag 
        mend (bool): end flag 
        size (int): data length in bytes 
        type (bytes): description of data entry 

    Returns:
        bytes: LIME header byte sequence 
    """
    m = 0
    if mbeg:
        m = 1<<15 # starting flag in binary
    elif mend:
        m = 1<<14 # ending flag in binary

    # '> ensures big-endian byte order
    # i packs LIMEMAGIC as a 4-byte integer
    # h packs the integer 1 as a 2-byte short (stands in for version number)
    # H packs m as a 2-byte unsigned short
    # q packs size as an 8-byte long long
    # 128s packs type as a 128-byte string (padded with null bytes if necessary)
    return struct.pack('>ihHq128s', LIMEMAGIC, 1, m, size, type)


def scidacChecksum(latdata, vol, sitesize):
    """
    Compute SciDAC checksum for lattice data
    
    Args:
        latdata (bytes)
        vol (int): number of sites on lattice 
        sitesize (int): size in bytes of data at each site 

    Returns:
        checksum values (suma, sumb) 
    """
    suma = 0
    sumb = 0
    for site in range(vol):
        # Compute the CRC32 checksum for this site
        c = zlib.crc32(latdata[site*sitesize:(site+1)*sitesize]) & 0xffffffff
        s29 = site%29
        s31 = site%31
        # Updates sum by XORing it with c rotated left by s* bits and right by (32-s*) bits
        suma ^= (c<<s29 | c>>(32-s29)) & 0xffffffff
        sumb ^= (c<<s31 | c>>(32-s31)) & 0xffffffff
    return suma, sumb


def xmlFind(s,xmlTag) -> bytes:
    """
    Find substring in byte string s between xmlTags.
    
    Args:
        s (bytes)
        xmlTag (bytes)

    Returns:
        bytes: byte string between tags 
    """
    rs = '<'+xmlTag+r'>([^<]*)</'+xmlTag+'>'
    m = re.findall(rs.encode('utf-8'),s)
    if len(m)!=1:
        logger.TBRaise(f'finding record {xmlTag} in unsupported xml: {s}')
    return m[0]


def trimNull(byteString) -> bytes:
    """
    Remove null byte and everything after it from a byte string.

    Args:
        byteString (bytes)

    Returns:
        bytes: trimmed byte string 
    """
    n = byteString.find(b'\0')
    if n>=0:
        return byteString[:n]
    else:
        return byteString
    

def printLimeHeaders(filename):
    """
    Print the headers of a file in LIME format.

    Args:
        filename (str)
    """
    checkType(str,filename=filename)
    f = open(filename,'rb')  
    while True:
        header = f.read(144)
        if not header:
            break
        magi, vers, mbeg_end_res, size, type = struct.unpack('>ihHq128s',header)
        if magi!=LIMEMAGIC:
            logger.TBRaise(f'lime magic number does not match, got {magi}')
        logger.info('* HEADER')
        logger.info(f'   lime version: {vers}')
        mbeg = 1 & mbeg_end_res>>15
        mend = 1 & mbeg_end_res>>14
        mres = ~(3<<14) & mbeg_end_res
        logger.info(f'  mbeg mend res: {mbeg}, {mend}, {mres}')
        type = trimNull(type)
        logger.info(f'           size: {size} bytes')
        logger.info(f'           type: {type}')
        logger.info(f'       position: currently {f.tell()} bytes from start of file')
        if type==b'scidac-binary-data' or type==b'ildg-binary-data':
            next = (size+7)//8*8
        else:
            data = f.read(size)
            data = trimNull(data)
            logger.info(f'           data: {data}')
            next = 7-(size+7)%8
        f.seek(next,os.SEEK_CUR)
    f.close()
