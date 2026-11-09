# Security

Quanty is a programming language. A Quanty script is a Lua program and, like a Python or shell script, it can read, write and delete files and run other programs, with the rights of the user who runs it. Quanty therefore does not try to protect the computer from the scripts it runs. A bug in Quanty, for example an out-of-bounds read or write in its C functions, cannot do anything that a script could not already do on purpose.

This means:

- Only run scripts that you trust, as you would with any other program.
- Run Quanty as an ordinary user, never as root or administrator.

## Reporting bugs

All bugs, including crashes and wrong memory access, are welcome on the forum at [www.quanty.org](https://www.quanty.org). Please include the Quanty version, the operating system, and the smallest script that shows the problem. Fixes go into the `develop` branch and appear in the next release.
