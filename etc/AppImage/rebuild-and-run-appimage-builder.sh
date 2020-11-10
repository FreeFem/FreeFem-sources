#!/usr/bin/env bash
#
# You need to prepare this docker container if you want to buitd AppImage locally.
# Normally it should be built with TravisCI (see .travis.yml)
#
# Install Docker following the instructions at
# https://docs.docker.com/install/linux/docker-ce/ubuntu/
# if you are using Ubuntu, or choose what's appropriate for your system:
# https://www.docker.com/get-started

set -x
set -u
set -e

# https://stackoverflow.com/questions/4774054/reliable-way-for-a-bash-script-to-get-the-full-path-to-itself
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

# Sources top directory
SRC_DIR=$(dirname $(dirname $SCRIPTPATH))

CONTAINER_NAME="freefem-appimage-builder"
IMAGE_NAME="freefem-appimage-builder"
SOURCES_MOUNT_POINT="/home/ubuntu/FreeFem-sources"

if docker inspect -f '{{.State.Running}}' "$CONTAINER_NAME" ; then
  docker stop "$CONTAINER_NAME"
  docker rm "$CONTAINER_NAME"
fi

docker build -t "$IMAGE_NAME" "$SCRIPTPATH/appimage-builder/" --network=host

# Run a container with interactive Bash as the main process (PID 1)
# with the FreeFem sources directory mounted:
docker run --detach \
           --interactive \
           --tty \
           --name "$CONTAINER_NAME" \
           --volume "$SRC_DIR:$SOURCES_MOUNT_POINT" \
           "$IMAGE_NAME"

docker ps

cat << EOF

  Docker container for AppImage build is running.

  To start build use the following command:

  docker exec $CONTAINER_NAME $SOURCES_MOUNT_POINT/etc/AppImage/build-appimage.sh

  If you want to attach interactive Bash shell inside the container use

  docker attach "$CONTAINER_NAME"

  Press CTRL-p CTRL-q to detach and leave the container running.

EOF
