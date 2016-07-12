INSTALL_DIR=$PWD
CPP_FL="--std=c++11 -march=native -O3"
C_FL="--std=gnu99 -march=native -O3"

# Install QMT
echo "Installing QMT..."
git clone --recursive https://github.com/usqcd-software/qmt.git
cd qmt
mkdir build
mkdir install
cd build
../configure --prefix=${INSTALL_DIR}/qmt/install/ CC=gcc-4.8 CXX=g++-4.8 CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
make
make install
cd ${INSTALL_DIR}

# Install QDP++
echo "Installing QDP++..."
git clone --recursive https://github.com/usqcd-software/qdpxx.git
cd qdpxx
autoreconf -f
mkdir scalar-build
mkdir scalar-install
cd scalar-build
../configure --prefix=${INSTALL_DIR}/qdpxx/scalar-install --with-qmt=${INSTALL_DIR}/qmt/install/ --enable-parallel-arch=scalar --enable-sse3 CC=gcc-4.8 CXX=g++-4.8 CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
make
make install
cd ${INSTALL_DIR}

# Install QUDA
echo "Installing QUDA"
git clone --recursive https://github.com/lattice/quda.git
cd quda
mkdir build
mkdir install
sed -i '401 a dev = 0;' ./lib/interface_quda.cpp
cd build
cmake -j4 .. -DQUDA_GPU_ARCH=sm_50 -DCMAKE_C_COMPILER=gcc-4.8 -DCMAKE_CXX_COMPILER=g++-4.8
make
cd ${INSTALL_DIR}

# Install Chroma
echo "Installing Chroma..."
git clone --recursive https://github.com/JeffersonLab/chroma.git
cd chroma
mkdir scalar-build
mkdir scalar-install
./autogen.sh
cd scalar-build
../configure --prefix=${INSTALL_DIR}/chroma/scalar-install/ --with-qdp=${INSTALL_DIR}/qdpxx/scalar-install/ --with-qmt=${INSTALL_DIR}/qmt/install/ --with-quda=${INSTALL_DIR}/quda/build/ --enable-sse3 CC=gcc-4.8 CXX=g++-4.8 CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
make -j8
make install -j4
cd ${INSTALL_DIR}

# Install Chroma utils
echo "Installing Chroma utils"
git clone --recursive https://github.com/JeffersonLab/chroma_utils.git
cd chroma_utils
mkdir build
mkdir install
autoreconf -f
automake --add-missing
cd build
../configure --prefix=${INSTALL_DIR}/chroma_utils/install/ CC=gcc-4.8 CXX=g++-4.8 CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
make
make install
cd ${INSTALL_DIR}

# Create symbolic links
sudo ln -s ${INSTALL_DIR}/chroma/scalar-install/bin/chroma /usr/local/bin/chroma
sudo ln -s ${INSTALL_DIR}/chroma/scalar-install/bin/purgaug /usr/local/bin/purgaug

echo "Done"

# Alternative flags for CXXFLAGS:  -O3 -finline-limit=50000  -funroll-all-loops
# Alternative chroma options: --enable-sse-wilson-dslash (not sure why it doesn't work)
