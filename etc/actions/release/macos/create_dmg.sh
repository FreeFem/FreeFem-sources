#!/bin/bash

set -e
set -x

PACKAGE_DIR=$1
DMG_BASENAME=$2

mkdir "$HOME/$DMG_BASENAME"
mv "$PACKAGE_DIR/Applications/FreeFem++.app" "$HOME/$DMG_BASENAME"
cp ./bin/script/Install-app.sh "$HOME/$DMG_BASENAME"

hdiutil create -srcfolder "$HOME/$DMG_BASENAME" -fs HFS+ "$DMG_BASENAME.dmg"
