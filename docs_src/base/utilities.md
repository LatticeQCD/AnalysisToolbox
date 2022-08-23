# Utilities

In this module we collect a bunch of useful, general utilities. These include:
- `getArg`: Python's argument parser is not very careful. This wrapper complains if any invalid arguments were passed.
- `printArg`: Will print an argument to screen, only if it isn't `None`.
- `printDict`: Prints a dictionary out with each key-value pair getting a line. This is more readable than the standard print.
- `shell`: Calls a Bash command from within Python and captures the output of the command.
- `shellVerbose`: Does the same thing, except prints the output to screen rather than capturing it.
- `find_nearest_idx`: Finds the index of an array whose element is closest to some specified number.

In addition, this module contains the `timer` class.