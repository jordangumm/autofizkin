hd $SRC_DIR

autoreconf --install --force
./configure --prefix=$PREFIX
make
make install
