INSTALL_DIR=$PWD
CPP_COMP="g++-4.8"
C_COMP="gcc-4.8"
CPP_FL="--std=c++11 -march=native -O3"
C_FL="--std=gnu99 -march=native -O3"

# Install QMT
echo "Installing QMT..."
git clone --recursive https://github.com/usqcd-software/qmt.git
cd qmt
mkdir build
mkdir install
cd build
../configure --prefix=${INSTALL_DIR}/qmt/install/ CC=${C_COMP} CXX=${CPP_COMP} CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
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
../configure --prefix=${INSTALL_DIR}/qdpxx/scalar-install --with-qmt=${INSTALL_DIR}/qmt/install/ --enable-parallel-arch=scalar --enable-sse3 CC=${C_COMP} CXX=${CPP_COMP} CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
make
make install
cd ${INSTALL_DIR}

# Install Chroma
echo "Installing Chroma..."
git clone --recursive https://github.com/JeffersonLab/chroma.git
cd chroma
mkdir scalar-build
mkdir scalar-install
./autogen.sh
cd scalar-build
../configure --prefix=${INSTALL_DIR}/chroma/scalar-install/ --with-qdp=${INSTALL_DIR}/qdpxx/scalar-install/ --with-qmt=${INSTALL_DIR}/qmt/install/ --enable-sse3 CC=${C_COMP} CXX=${CPP_COMP} CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
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
../configure --prefix=${INSTALL_DIR}/chroma_utils/install/ CC=${C_COMP} CXX=${CPP_COMP} CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
make
make install
cd ${INSTALL_DIR}

# Create symbolic links
sudo ln -s ${INSTALL_DIR}/chroma/scalar-install/bin/chroma /usr/local/bin/chroma
sudo ln -s ${INSTALL_DIR}/chroma/scalar-install/bin/purgaug /usr/local/bin/purgaug

echo "Done"

# Alternative flags for CXXFLAGS:  -O3 -finline-limit=50000  -funroll-all-loops
# Alternative chroma options: --enable-sse-wilson-dslash (not sure why it doesn't work)
