INSTALL_DIR=$PWD
CPP_FL="--std=c++11 -march=native -Os"
C_FL="--std=gnu99 -march=native -Os"
CUDA_DIR="/usr/local/cuda/"

# # Install QMT
# echo "Installing QMT..."
# git clone --recursive https://github.com/usqcd-software/qmt.git
# cd qmt
# mkdir build
# mkdir install
# cd build
# ../configure --prefix=${INSTALL_DIR}/qmt/install/ CC=gcc-4.8 CXX=g++-4.8 CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
# make
# make install
# cd ${INSTALL_DIR}

# # Install QMP
# echo "Installing QMP..."
# git clone --recursive https://github.com/usqcd-software/qmp.git
# cd qmp
# mkdir build
# mkdir install
# cd build
# ../configure --prefix=${INSTALL_DIR}/qmp/install/ CC=gcc-4.8 CXX=g++-4.8 CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
# make
# make install
# cd ${INSTALL_DIR}

# # Install LLVM
# echo "Installing LLVM..."
# wget http://llvm.org/releases/3.4/llvm-3.4.src.tar.gz
# tar -zxvf llvm-3.4.src.tar.gz
# rm llvm-3.4.src.tar.gz
# mv llvm-3.4 llvm
# cd llvm
# mkdir build
# mkdir install
# cd build
# ../configure --prefix=${INSTALL_DIR}/llvm/install/ CC=gcc-4.8 CXX=g++-4.8 CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}"
# make -j8
# make install -j4
# cd ${INSTALL_DIR}

# # Install QDP-JIT
# echo "Installing QDP-JIT..."
# git clone -b llvm-nvvm --recursive https://github.com/fwinter/qdp-jit.git
# cd qdp-jit
# autoreconf -f
# mkdir parscalar-build
# mkdir parscalar-install
# cd parscalar-build
# ../configure --prefix=${INSTALL_DIR}/qdp-jit/parscalar-install --with-qmt=${INSTALL_DIR}/qmt/install/ --enable-parallel-arch=parscalar --enable-sse3 CC=gcc-4.8 CXX=g++-4.8\
# 			 --with-qmp=${INSTALL_DIR}/qmp/install/ CXXFLAGS="${CPP_FL} -I${INSTALL_DIR}/qmp/install/" CFLAGS="${C_FL} -I${INSTALL_DIR}/qmp/install/"\
# 			 --enable-largefile --enable-parallel-io --enable-dml-output-buffering --disable-generics --with-llvm=${INSTALL_DIR}/llvm/install/\
# 			 --with-nvvm=${CUDA_DIR}/nvvm/ --with-cuda=${CUDA_DIR} 
# make
# make install
# cd ${INSTALL_DIR}

# # Install QUDA
# echo "Installing QUDA"
# git clone --recursive https://github.com/lattice/quda.git
cd quda
mkdir build
mkdir install
sed -i '401 a dev = 0;' ./lib/interface_quda.cpp
cd build
cmake -j4 .. -DQUDA_GPU_ARCH=sm_50 -DCMAKE_C_COMPILER=gcc-4.8 -DCMAKE_CXX_COMPILER=g++-4.8 -DQUDA_INTERFACE_QDPJIT=ON -DQUDA_QDPJIT=ON
make
cd ${INSTALL_DIR}

# Install Chroma
echo "Installing Chroma..."
git clone --recursive https://github.com/JeffersonLab/chroma.git
cd chroma
mkdir parscalar-build
mkdir parscalar-install
./autogen.sh
cd parscalar-build
../configure --prefix=${INSTALL_DIR}/chroma/parscalar-install/ --with-qdp=${INSTALL_DIR}/qdp-jit/parscalar-install/ --with-qmt=${INSTALL_DIR}/qmt/install/\
 			 --with-quda=${INSTALL_DIR}/quda/build/ --enable-sse3 CC=gcc-4.8 CXX=g++-4.8 CXXFLAGS="${CPP_FL}" CFLAGS="${C_FL}" --enable-jit-clover --enable-quda-deviface
# add  "-ldl -lm" to the last append in the mainprogs makefiles
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
sudo ln -s ${INSTALL_DIR}/chroma/parscalar-install/bin/chroma /usr/local/bin/chroma-jit
sudo ln -s ${INSTALL_DIR}/chroma/parscalar-install/bin/purgaug /usr/local/bin/purgaug-jit

echo "Done"

# Alternative flags for CXXFLAGS:  -O3 -finline-limit=50000  -funroll-all-loops
# Alternative chroma options: --enable-sse-wilson-dslash (not sure why it doesn't work)
