#!/bin/sh -x
# Clean files that can be regenerated, to avoid CVS conflicts

# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr
# $Id$

find . -name Makefile.in -exec rm {} \;

rm HISTORY
rm config.h.in
rm configure
rm examples++-tutorial/all.edp
rm examples++/all.edp
rm examples++/regtests.edp
rm src/fflib/strversionnumber.cpp
rm src/ide/hl_lex.c++
rm src/ide/hl_yacc.c++
rm src/ide/hl_yacc.h
rm src/lglib/lg.tab.?pp
