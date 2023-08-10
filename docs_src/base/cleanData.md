# Cleaning, splicing, and organizing arrays

You can often think of all the data in a project you are dealing with as a giant table.
For instance the table will be indexed by things like configuration number, ensemble stream/series,
quark mass, and so on. Sometimes it is important to find data that lies at the intersection
of some sets, or to look only at data whose bare coupling $\beta$ fall within some prescribed
range, and so on. The methods in
```Python
latqcdtools.base.cleanData
```
try to help with this.

These manipulations are already possible within numpy, but they can be difficult to read. (At least
they are difficult for David to read.) Hence this module contains wrappers for such manipulations.
They include
- `deleteRow`: Delete a row of a 2-d array.
- `deleteCol`: Delete a column of a 2-d array.
- `clipRange`: Pick a column `col` of a data table `data`. Throw out all data whose corresponding
entries in `col` fall outside the range [`minVal`,`maxVal`].
- `intersectAtCol`: Return only those rows of `table1` and `table2` that have 
identical elements in column `col`.
- `spliceAtCol`:  Assuming two tables `table1` and `table2` have common values in column `col`, 
create a new table, where `table1` has corresponding entries less than `atVal` in `col`, 
and `table2` has corresponding entries greater than `atVal.
- `restrictAtCol`: Return only those rows of a table where `col` has exactly the value `atVal`.
- `excludeAtCol`: Return everything except those rows of a table where `col` has exactly the 
value `atVal`.
