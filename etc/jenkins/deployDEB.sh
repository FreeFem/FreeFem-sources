#!/bin/sh

TOKEN=$1
ORGANIZATION="FreeFem"
REPOSITORY="FreeFem-sources"

LATEST_RELEASE=`curl 'https://api.github.com/repos/'$ORGANIZATION'/'$REPOSITORY'/releases/latest'`
CURRENT_VERSION=`echo ${LATEST_RELEASE}| jq -r '.tag_name'`
UPLOAD_URL=`echo ${LATEST_RELEASE}| jq -r '.upload_url'`
FILE_NAME='freefem_'$VERSION'-1_amd64.deb'

git reset --hard
git clean -fdx

autoreconf -i \
    && ./configure --enable-download --enable-optim --enable-generic \
    && ./3rdparty/getall -a \
    && make -j4 \
    && echo "FreeFEM: Finite Element Language" > description-pak \
    && sudo checkinstall -D --install=no \
        --pkgname "freefem" --pkgrelease "1" \
        --pkgversion "${VERSION}" --pkglicense "LGPL-2+" \
        --pkgsource "https://github.com/FreeFem/FreeFem-sources" \
        --pkgaltsource "https://freefem.org/" \
        --maintener "FreeFEM" --backup=no --default \
    && RESPONSE=`curl --data-binary "@$FILE_NAME" -H "Authorization: token $TOKEN" -H "Content-Type: application/octet-stream" "$UPLOAD_URL=$FILE_NAME"`
