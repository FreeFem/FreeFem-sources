#!/bin/sh
# Copies Debian packages into a web server directory structure

# Any error stops the script
set -e

# Displays all commands
set -v

# $1 must point to the directory where the files will be copied
if ssh $3 test ! -d $1
    then
    echo "Usage: CopyToServer.sh <directory to copy files to> <testing/unstable/debian64> <remote server location>"
    exit 1
fi

# $2=testing/unstable/debian64
if ssh $3 test "$2" != testing -a "$2" != unstable -a "$2" != debian64
    then
    echo "Usage: CopyToServer.sh <directory to copy files to> <testing/unstable/debian64> <remote server location>"
    exit 1
fi

# Directory structure we want to build
if ssh $3 test -d $1/dists/$2
then
    ssh $3 rm -r $1/dists/$2
fi
if test "$HOSTNAME" = iris
then
    basedir="dists/$2/ff++/binary-amd64"
else
    basedir="dists/$2/ff++/binary-i386"
fi
ssh $3 mkdir -p $1/$basedir

# Server configuration file
scp apt-ftparchive.conf $3:$1

# Copies all packages (potential improvement: we could spread
# architecture-independant packages into a separate directory).

scp ../../freefem++*.{dsc,deb,tar.gz} $3:$1/$basedir

# Create package list
ssh $3 "cd $1 && apt-ftparchive packages ./$basedir > ./$basedir/Packages"
ssh $3 "cd $1 && gzip -c ./$basedir/Packages > ./$basedir/Packages.gz"
