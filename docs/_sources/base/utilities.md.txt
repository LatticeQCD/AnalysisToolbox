# Utilities

In this module we collect a bunch of useful, general utilities.

## Arguments and formatted output

- `getArg`: Python's argument parser is not very careful. This wrapper complains if any invalid arguments were passed.
- `printArg`: Will print an argument to screen, only if it isn't `None`.
- `printDict`: Prints a dictionary out with each key-value pair getting a line. This is more readable than the standard print.
- `cleanOutput`: Often one wants to print multiple different numbers on one line, for example a row of a table. It helps
readability if all these use the same indentation and have the same width on the screen. You can create formatted strings
that have the same convention by passing values separated by commas to `cleanOutput`.
- `printClean`: This is a wrapper for `cleanOutput` that puts the result to screen instead of giving you a string. (You may
want the string if, e.g., you are writing a custom output to file.)

## Interaction with shell

- `shell`: Calls a Bash command from within Python and captures the output of the command.
- `shellVerbose`: Does the same thing, except prints the output to screen rather than capturing it.
- `deleteFile`: Deletes a file, if it exists.

## Playing with arrays

- `find_nearest_idx`: Finds the index of an array whose element is closest to some specified number.
- `convertToNumpy`: Converts all arguments to numpy arrays.
- `isArrayLike`: Determine whether argument is array-like.
- `isHigherDimensional`: Determine whether argument is array-like with at least two indices.
- `unvector`: Remove outermost brackets of an array-like object, if possile.
- `envector`: Convert a scalar to a numpy array. `envector` should undo an `unvector` and vice-versa.
- `naturalSort`: Sort list of strings so that, e.g. '10' comes after '9' rather than before it.

In addition, this module contains the `timer` class.