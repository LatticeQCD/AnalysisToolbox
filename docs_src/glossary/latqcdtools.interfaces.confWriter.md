latqcdtools.interfaces.confWriter
=============

```python
writeConf(fileName, gauge, format='3x3', endianness='>', precision='d'):
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
```
```python
class NERSCWriter(gauge):
"""
Writer for NERSC format gauge configurations.

Parameters
----------
gauge : gaugeField
    SU(3) gauge field to be written.
"""
```
```python
class confWriter(gauge):
"""
Base class for writing gauge field configurations.

Parameters
----------
gauge : gaugeField
    SU(3) gauge field to be written.
"""
```
