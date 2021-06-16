# -----------------------------------------------------------------
# $Revision: 1.9 $
# $Date: 2009/03/25 23:10:50 $
# -----------------------------------------------------------------
# Programmer(s): Radu Serban and Aaron Collier @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# Makefile for FKINSOL module
#
# src/kinsol/fcmix/Makefile.  Generated from Makefile.in by configure.
# -----------------------------------------------------------------

SHELL = /bin/bash

srcdir       = .
builddir     = .
abs_builddir = /home/marta/Documents/3D-PDR-master/sundials-2.5.0/src/kinsol/fcmix
top_builddir = ../../../
prefix       = /usr/local
exec_prefix  = ${prefix}
includedir   = ${prefix}/include
libdir       = ${exec_prefix}/lib

INSTALL        = /usr/bin/install -c
INSTALL_LIB    = ${INSTALL}
INSTALL_HEADER = ${INSTALL} -m 644

LIBTOOL      = $(SHELL) $(top_builddir)/libtool
LIBTOOL_DEPS = config/ltmain.sh

CPP      = cc -E
CPPFLAGS = 
CC       = cc
CFLAGS   = -g -O2
LDFLAGS  = 
LIBS     = -lm 

LAPACK_ENABLED = no

top_srcdir = $(srcdir)/../../..

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/src/kinsol -I$(top_builddir)/include

LIB_REVISION = 0:1:0

FKINSOL_LIB = libsundials_fkinsol.la

FKINSOL_SRC_FILES = fkinsol.c fkindense.c fkinband.c fkinpreco.c fkinjtimes.c fkinbbd.c
FKINSOL_BL_SRC_FILES = fkinlapack.c fkinlapdense.c fkinlapband.c 

FKINSOL_OBJ_FILES = $(FKINSOL_SRC_FILES:.c=.o)
FKINSOL_BL_OBJ_FILES = $(FKINSOL_BL_SRC_FILES:.c=.o)

FKINSOL_LIB_FILES = $(FKINSOL_SRC_FILES:.c=.lo)
FKINSOL_BL_LIB_FILES = $(FKINSOL_BL_SRC_FILES:.c=.lo)

mkinstalldirs = $(SHELL) $(top_srcdir)/config/mkinstalldirs

# ----------------------------------------------------------------------------------------------------------------------

all: $(FKINSOL_LIB)

$(FKINSOL_LIB): $(FKINSOL_LIB_FILES)
	@if test "X${LAPACK_ENABLED}" = "Xyes"; then \
	  make lib_with_bl;                          \
	else                                         \
	  make lib_without_bl;                       \
	fi

lib_without_bl: $(FKINSOL_LIB_FILES)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) -o $(FKINSOL_LIB) $(FKINSOL_LIB_FILES) -rpath $(libdir) $(LDFLAGS) $(LIBS) -static -version-info $(LIB_REVISION)

lib_with_bl: $(FKINSOL_LIB_FILES) $(FKINSOL_BL_LIB_FILES)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) -o $(FKINSOL_LIB) $(FKINSOL_LIB_FILES) $(FKINSOL_BL_LIB_FILES) -rpath $(libdir) $(LDFLAGS) $(LIBS) -static -version-info $(LIB_REVISION)

install: $(FKINSOL_LIB)
	$(mkinstalldirs) $(libdir)
	$(LIBTOOL) --mode=install $(INSTALL_LIB) $(FKINSOL_LIB) $(libdir)

uninstall:
	$(LIBTOOL) --mode=uninstall rm -f $(libdir)/$(FKINSOL_LIB)

clean:
	$(LIBTOOL) --mode=clean rm -f $(FKINSOL_LIB)
	rm -f $(FKINSOL_LIB_FILES)
	rm -f $(FKINSOL_BL_LIB_FILES)
	rm -f $(FKINSOL_OBJ_FILES)
	rm -f $(FKINSOL_BL_OBJ_FILES)

distclean: clean
	rm -f Makefile

fkinsol.lo: $(srcdir)/fkinsol.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkinsol.c
fkinpreco.lo: $(srcdir)/fkinpreco.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkinpreco.c
fkinjtimes.lo: $(srcdir)/fkinjtimes.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkinjtimes.c
fkinbbd.lo: $(srcdir)/fkinbbd.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkinbbd.c
fkindense.lo: $(srcdir)/fkindense.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkindense.c
fkinband.lo: $(srcdir)/fkinband.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkinband.c
fkinlapack.lo: $(srcdir)/fkinlapack.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkinlapack.c
fkinlapdense.lo: $(srcdir)/fkinlapdense.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkinlapdense.c
fkinlapband.lo: $(srcdir)/fkinlapband.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fkinlapband.c

libtool: $(top_builddir)/$(LIBTOOL_DEPS)
	@cd ${top_builddir} ;                 \
	 ${SHELL} ./config.status --recheck ; \
	 cd ${abs_builddir}
