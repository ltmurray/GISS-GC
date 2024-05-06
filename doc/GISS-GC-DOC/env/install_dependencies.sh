#!/bin/bash

# Go to home directory
cd $HOME/opt/gcc
rm -rf bin include lib share src

# Make Installation directories
mkdir -p src

# Set Installation Path 
export INSTALL_DIR=$HOME/opt/gcc

export LD_LIBRARY_PATH=$INSTALL_DIR/lib:$LD_LIBRARY_PATH     # Linux
export DYLD_LIBRARY_PATH=$INSTALL_DIR/lib:$DYLD_LIBRARY_PATH # macOS

export CPPFLAGS="-fPIC -I${INSTALL_DIR}/include"
export LDFLAGS="-L${INSTALL_DIR}/lib"

# Downloading zlib
cd $INSTALL_DIR/src
wget -c -4 https://github.com/madler/zlib/archive/refs/tags/v1.2.12.tar.gz
# Unpack the tar ﬁle
tar -xvzf v1.2.12.tar.gz
cd zlib-1.2.12/
# Conﬁgure and Install Zlib
CC= CXX= ./configure --prefix=$INSTALL_DIR
make -j
make install
make -j check

export FC=mpif90
export F90=mpif90
export F77=mpif90
export CC=mpicc
export CXX=mpicxx

# Downloading HDF5
cd $INSTALL_DIR/src
wget -c -4 https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-1_14_3.tar.gz
# Unpack the tar ﬁle
tar -xvzf hdf5-1_14_3.tar.gz
cd hdf5-hdf5-1_14_3
# Conﬁgure and Install HDF5
./configure --prefix=$INSTALL_DIR --with-zlib=$INSTALL_DIR --enable-hl --enable-fortran --enable-parallel
make -j
make install
#make -j check

# Downloading netCDF
cd $INSTALL_DIR/src
wget -c -4 https://github.com/Unidata/netcdf-c/archive/refs/tags/v4.9.2.tar.gz
# Unpack the tar ﬁle
tar -xzvf v4.9.2.tar.gz
cd netcdf-c-4.9.2/
# Conﬁgure and Install netCDF-C
./configure --prefix=$INSTALL_DIR --disable-dap
make -j 
make install
#make check
# Check netCDF-C installation parameters
nc-conﬁg --version
#nc-conﬁg --all

# Downloading netCDF-fortran
cd $INSTALL_DIR/src
wget -c -4 https://github.com/Unidata/netcdf-fortran/archive/refs/tags/v4.6.0.tar.gz
# Unpack the tar ﬁle
tar -xvzf v4.6.0.tar.gz
cd netcdf-fortran-4.6.0/
# Conﬁgure and Install netCDF-fortran
./configure --prefix=$INSTALL_DIR --disable-shared
make -j
make install
#make check
# Check netCDF fortran installation parameters
#nf-config --all

# Downloading PnetCDF
cd $INSTALL_DIR/src
wget -c -4  https://parallel-netcdf.github.io/Release/pnetcdf-1.12.3.tar.gz
# Unpack the tar ﬁle
tar -xvzf pnetcdf-1.12.3.tar.gz
cd pnetcdf-1.12.3
# Setting the environment 
export MPICC=mpicc
export MPICXX=mpicxx
export MPIFC=mpifort
CC=$MPICC 
CXX=$MPICXX 
FC=$MPIFC 
F77=$MPIFC
# Conﬁgure and Install PnetCDF
./configure --prefix=$INSTALL_DIR --enable-shared --enable-fortran --enable-profiling --enable-large-file-test
make -j
#make tests
#make check
#make ptest
#make ptests
make install

exit 0;
