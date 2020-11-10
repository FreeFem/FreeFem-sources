#!/usr/bin/env bash
#
# The purpose of this custom AppRun script is
# to allow symlinking the AppImage and invoking
# the corresponding binary depending on which
# symlink was used to invoke the AppImage
#

if [ -n "$APPIMAGE_DEBUG" ] ; then
  set -x
fi

# support for running extracted AppImage like this:
# $ ./FreeFem-x86_64.AppImage --appimage-extract
# $ cd squashfs-root
# $ ./AppRun FreeFem++ ../mycode.edp
if [ -z "$APPDIR" ] ; then
    APPDIR=`dirname -- "$0"`
fi

if [ -z "$1" -o "$1" = "--help" ] ; then
cat << EOF

   This is AppImage version of FreeFem.

   To run individual tools specify them as the first argument:

   ./FreeFem.AppImage FreeFem++ mycode.edp

EOF
  exit 0
fi

export PATH=$APPDIR/usr/bin:$PATH

BINARY_NAME="$1"
shift

if [ -e "$APPDIR/usr/bin/$BINARY_NAME" ] ; then
  if [ -n "$APPIMAGE_DEBUG" ] ; then
    exec strace -o log -s 120 -f "$APPDIR/usr/bin/$BINARY_NAME" "$@"
  else
    exec "$APPDIR/usr/bin/$BINARY_NAME" "$@"
  fi
else
  echo "Error: Tool $1 is not a part of FreeFem package"
  exit 2
fi
