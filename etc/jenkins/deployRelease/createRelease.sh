#!/bin/sh

## Parameters
TOKEN=$1
ORGANIZATION="FreeFem"
REPOSITORY="FreeFem-sources"

## Release informations
VERSION=`grep AC_INIT configure.ac | cut -d"," -f2 | tr - .`
RELEASE_TAG_NAME="v$VERSION"
RELEASE_TARGET_COMMITISH="master"
RELEASE_NAME="FreeFEM v$VERSION"
RELEASE_BODY="**Warning**: this is an automatic release. If you encounter some trouble with packages, post a message in the [forum](https://community.freefem.org)"
RELEASE_DRAFT="false"
RELEASE_PRERELEASE="false"

echo "Release name: $RELEASE_NAME"

## Check if release exists
RESPONSE=`curl 'https://api.github.com/repos/'$ORGANIZATION'/'$REPOSITORY'/releases/tags/'$RELEASE_TAG_NAME`
RESPONSE_URL=`printf "%s" "$RESPONSE" | jq -r '.url'`

## Create release
if [ -z "$RESPONSE_URL" ]
then
	RELEASE_PARAMETERS=$(printf '{"tag_name": "%s", "target_commitish": "%s", "name": "%s", "body": "%s.", "draft": %s, "prerelease": %s}' "$RELEASE_TAG_NAME" "$RELEASE_TARGET_COMMITISH" "$RELEASE_NAME" "$RELEASE_BODY" "$RELEASE_DRAFT" "$RELEASE_PRERELEASE")
	RELEASE=`curl -H "Authorization: token $TOKEN" --data "$RELEASE_PARAMETERS" 'https://api.github.com/repos/'$ORGANIZATION'/'$REPOSITORY'/releases'`
else
	echo "Release already exists"
	exit 1
fi
