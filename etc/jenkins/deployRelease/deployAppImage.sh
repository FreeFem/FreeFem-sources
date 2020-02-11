#!/bin/sh

set -x
set -u
set -e

## Parameters
TOKEN=$1
ORGANIZATION="FreeFem"
REPOSITORY="FreeFem-sources"
VERSION=`grep AC_INIT configure.ac | cut -d"," -f2`
RELEASE_TAG_NAME="v$VERSION"
APPIMAGE_NAME="FreeFEM-x86_64-${RELEASE_TAG_NAME}.AppImage"

## AppImage build
CONTAINER_NAME="freefem-appimage-builder"
SOURCES_MOUNT_POINT="/home/ubuntu/FreeFem-sources"

chmod +x etc/AppImage/rebuild-and-run-appimage-builder.sh && ./etc/AppImage/rebuild-and-run-appimage-builder.sh
docker exec $CONTAINER_NAME $SOURCES_MOUNT_POINT/etc/AppImage/build-appimage.sh
# AppImage is created: FreeFEM-x86_64.AppImage

## Change name
mv FreeFEM-x86_64.AppImage $APPIMAGE_NAME

## Deploy in GitHub release
RELEASE=`curl 'https://api.github.com/repos/'$ORGANIZATION'/'$REPOSITORY'/releases/tags/'$RELEASE_TAG_NAME`
UPLOAD_URL=`printf "%s" "$RELEASE" | jq -r '.upload_url'`

if [ -x $UPLOAD_URL ]
then
	echo "Release does not exists"
	exit 1
else
	RESPONSE=`curl --data-binary "@$APPIMAGE_NAME" -H "Authorization: token $TOKEN" -H "Content-Type: application/octet-stream" "$UPLOAD_URL=$APPIMAGE_NAME"`
fi