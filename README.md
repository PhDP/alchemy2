Alchemy 2
=========
[![Build status](https://travis-ci.org/PhDP/alchemy2.svg?branch=master)](https://travis-ci.org/PhDP/alchemy2)

A copy of the [official SVN repository](http://code.google.com/p/alchemy-2/)
with very very few modifications (see *Changes*).

Linux x86\_64 binaries can be found [here](https://github.com/PhDP/alchemy2/releases).

Installation
------------
Alchemy was developed/tested on Linux (Fedora Core 5) and... it fails to
compile with modern versions of gcc, clang, and bison.

I've tested this version of the code with with Bison 2 and g++ 4.4, see
*install.sh* for a simple install script used on Ubuntu to get the
right versions of bison and gcc.

There is no target for tests, so, yeah...

Changes
-------
Since alchemy is the reference implementation for markov logic, I tried to make
as few changes as possible. I only made a few modifications to simplify
compilation:

* I added files for a full example in *examples*.
* I added install.sh, README.md, .travis.cl, and .gitignore.
* I added the manual and tutorial PDF files in the *doc* folder.

License
-------
[MIT](http://opensource.org/licenses/MIT).

