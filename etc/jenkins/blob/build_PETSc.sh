#!/bin/bash

cd 3rdparty/ff-petsc \
	&& make petsc-slepc \
	&& cd - \
	&& ./reconfigure
