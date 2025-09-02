latqcdtools.interfaces.lime
=============

```Python
limeHeader(mbeg, mend, size, type) -> bytes:
'''
Get the header for a LIME record.

Args:
    mbeg (bool): start flag 
    mend (bool): end flag 
    size (int): data length in bytes 
    type (bytes): description of data entry 

Returns:
    bytes: LIME header byte sequence 
'''
```
```Python
printLimeHeaders(filename):
'''
Print the headers of a file in LIME format.

Args:
    filename (str)
'''
```
```Python
scidacChecksum(latdata, vol, sitesize):
'''
Compute SciDAC checksum for lattice data

Args:
    latdata (bytes)
    vol (int): number of sites on lattice 
    sitesize (int): size in bytes of data at each site 

Returns:
    checksum values (suma, sumb) 
'''
```
```Python
trimNull(byteString) -> bytes:
'''
Remove null byte and everything after it from a byte string.

Args:
    byteString (bytes)

Returns:
    bytes: trimmed byte string 
'''
```
```Python
xmlFind(s, xmlTag) -> bytes:
'''
Find substring in byte string s between xmlTags.

Args:
    s (bytes)
    xmlTag (bytes)

Returns:
    bytes: byte string between tags 
'''
```
