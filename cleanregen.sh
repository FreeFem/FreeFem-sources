#!/bin/sh -x
# Clean files that can be regenerated, to avoid CVS conflicts

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr
# $Id$

find . -name Makefile.in -exec rm {} \;
rm configure
rm src/fflib/strversionnumber.cpp
rm examples++-tutorial/StokesUzawa.edp
rm config.h.in
rm src/ide/hl_lex.c++
rm src/ide/hl_yacc.h
rm src/ide/hl_yacc.c++
