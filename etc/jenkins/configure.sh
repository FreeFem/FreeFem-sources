#!/bin/sh

autoreconf -i \
	&& ./configure --enable-download --enable-optim --prefix=/builds/freefem \
	&& ./3rdparty/getall -a
