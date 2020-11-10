#!/usr/bin/env bash
#
# This script configures, compiles and package FreeFem into AppImage
#
# This script normally run from .travis.yml
# but can be run locally in Docker as well (see rebuild-and-run-appimage-builder.sh)

set -x
set -u
set -e

# https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

# Sources top directory
SRC_DIR=$(dirname $(dirname $SCRIPTPATH))

cd $SRC_DIR

# make sure we are in the right directory
if [ ! -f "readme/INNOVATION" ] ; then
  echo "Build script error: unable to change dir into FreeFem sources top directory"
  exit 1
fi


TOOLDIR="$HOME/opt"

# Download AppImage Linux Deploy tool
LINUX_DEPLOY_TOOL="$TOOLDIR/linuxdeploy-x86_64.AppDir/AppRun"
APPIMAGE_TOOL="$TOOLDIR/appimagetool-x86_64.AppDir/AppRun"

function appimage_get_and_extract()
{
  local url="$1"

  local appimage=$(basename $url)
  local appname="${appimage%.*}"

  mkdir -p "$HOME/opt"
  pushd "$HOME/opt"
  rm -rf "squashfs-root"
  rm -rf "${appname}.AppDir"
  wget -c $url
    chmod +x "$appimage"
    ./$appimage --appimage-extract
    rm -f "$appimage"
    mv squashfs-root "${appname}.AppDir"
  popd
}

if [ ! -e "$LINUX_DEPLOY_TOOL" ] ; then
  appimage_get_and_extract \
    "https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage"
fi

if [ ! -e "$APPIMAGE_TOOL" ] ; then
  appimage_get_and_extract \
    "https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage"
fi


# TODO: when migration from Autoconf to CMake will be finished in FreeFem project,
# we need to update configure and build commands here to something like:
# cmake . -DCMAKE_INSTALL_PREFIX=/usr
# make
# make install DESTDIR=$SRC_DIR/AppDir
# https://docs.appimage.org/packaging-guide/from-source/native-binaries.html#cmake

autoreconf -i
./configure \
    CFLAGS='-fpermissive' \
    CXXFLAGS='-fpermissive' \
    --prefix=$SRC_DIR/AppDir/usr \
    --enable-download \
    --enable-optim \
    --enable-generic
./3rdparty/getall -a
make
make install

# TODO: ffglut tool is missing now
# To try to add it, add --with-glut to ./configure
# and add 'freeglut3-dev' to apt-get install list.

# Remove some extras that not need in AppImage
rm -rf $SRC_DIR/AppDir/usr/lib/ff++/4.1/include
rm -rf $SRC_DIR/AppDir/usr/lib/ff++/4.1/examples
rm -rf $SRC_DIR/AppDir/usr/share/freefem++/4.1/examples
rm -rf $SRC_DIR/AppDir/usr/share/doc/


# Create additional files need to AppImage
mkdir -p AppDir/usr/share/applications/
cat << EOF > AppDir/usr/share/applications/FreeFem.desktop
[Desktop Entry]
Name=FreeFEM
Type=Application
Exec=freefem-apprun.sh
Icon=freefem
Comment=Easy to use PDE solver
Categories=Science;
Terminal=true
EOF

mkdir -p AppDir/usr/share/icons/hicolor/32x32/apps
cp etc/logo/logo.png AppDir/usr/share/icons/hicolor/32x32/apps/freefem.png

mkdir -p AppDir/usr/bin
cp etc/AppImage/freefem-apprun.sh AppDir/usr/bin/freefem-apprun.sh

# Use it if you want resulting file to be FreeFem-3f71f1f-x86_64.AppImage for ex.
# export VERSION=$(git rev-parse --short HEAD)

$LINUX_DEPLOY_TOOL --appdir AppDir

# the output will be FreeFem-x86_64.AppImage
env ARCH=x86_64 $APPIMAGE_TOOL AppDir
