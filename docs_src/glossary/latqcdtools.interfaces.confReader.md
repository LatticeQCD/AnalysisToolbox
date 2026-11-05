latqcdtools.interfaces.confReader
=============

```python
_cubicDims(dims, fileName):
"""
Turn [Nx, Ny, Nz, Nt] into (Ns, Nt), checking that the spatial volume is cubic.
"""
```
```python
_ildgDims(fileName):
"""
Get [Nx, Ny, Nz, Nt] from the scidac-private-file-xml or ildg-format record of a LIME file.
"""
```
```python
detectConf(fileName):
"""
Detect the format and lattice size of a configuration from its contents. NERSC files start with BEGIN_HEADER,
ILDG files with the LIME magic number, and dyniso files are one of comp, MILC, or ascii, as described in
DynisoReader. All of these store the lattice dimensions near the start of the file.

Args:
    fileName (str): path to the configuration

Returns:
    tuple: (format, Ns, Nt), where format is 'NERSC', 'ILDG', or 'dyniso'
"""
```
```python
nerscChecksum(data, precision) -> int:
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
```
```python
class DynisoReader(Ns, Nt, nproc=None):
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
```
```python
class ILDGReader(Ns, Nt, nproc=None):
```
```python
class NERSCReader(Ns, Nt, nproc=None):
```
```python
class confReader(Ns, Nt, nproc=None):
"""
Base class for reading configurations. 
"""
```
