#!/bin/bash

# Install dependencies
sudo apt-get install g++-4.4 flex

# Make the folder to put the bin
mkdir -p bin/obj

# Download & install Bison 2
wget http://ftp.gnu.org/gnu/bison/bison-2.0.tar.gz
tar -xvf bison-2.0.tar.gz
cd bison-2.0
./configure
make
sudo make install
cd ..

# Compile Alchemy
cd src
make depend
make

