# Ensure repositories are up-to-date
apt-get update
# Install debian packages for dependencies
apt-get install -y \
    autoconf \
    automake \
    cryptsetup \
    fuse2fs \
    git \
    fuse \
    libfuse-dev \
    libglib2.0-dev \
    libseccomp-dev \
    libtool \
    pkg-config \
    runc \
    squashfs-tools \
    squashfs-tools-ng \
    uidmap \
    wget \
    zlib1g-dev #\
#    make

export VERSION=1.23.2 OS=linux ARCH=amd64  # change this as you need

wget -O /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz \
  https://dl.google.com/go/go${VERSION}.${OS}-${ARCH}.tar.gz
tar -C /usr/local -xzf /tmp/go${VERSION}.${OS}-${ARCH}.tar.gz

echo 'export PATH=$PATH:/usr/local/go/bin' >> ~/.bashrc
source ~/.bashrc

git clone --recurse-submodules https://github.com/sylabs/singularity.git
cd singularity

git checkout --recurse-submodules v4.2.1

./mconfig
make -C builddir
make -C builddir install

singularity --version