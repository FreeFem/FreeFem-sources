#!/bin/bash

set -e
set -x

PACKAGE_DIR=$1
DMG_BASENAME=$2
VERSION=$3

mkdir -p "$HOME/$DMG_BASENAME"
mv "$PACKAGE_DIR/Applications/FreeFem++.app" "$HOME/$DMG_BASENAME"
cp ./bin/script/Install-app.sh "$HOME/$DMG_BASENAME"
FFTESTDIR=$(find . -name Laplace.edp)
sed -e "s/@VV@/$VERSION/" -e "s%@APPFF@%$DMG_BASENAME%" -e s%@FFTEST@%$FFTESTDIR% <./bin/script/README.md.in >"$HOME/$DMG_BASENAME/README.md"

hdiutil create -srcfolder "$HOME/$DMG_BASENAME" -fs HFS+ "$DMG_BASENAME.dmg"
rm -rf "$HOME/$DMG_BASENAME"
