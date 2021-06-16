# -----------------------------------------------------------------
# $Revision: 1.8 $
# $Date: 2009/03/25 23:10:50 $
# -----------------------------------------------------------------
# Programmer(s): Aaron Collier @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2005, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# Makefile for FIDA module
#
# src/ida/fcmix/Makefile.  Generated from Makefile.in by configure.
# -----------------------------------------------------------------

SHELL = /bin/bash

srcdir       = .
builddir     = .
abs_builddir = /home/marta/Documents/3D-PDR-master/sundials-2.5.0/src/ida/fcmix
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

INCLUDES = -I$(top_srcdir)/include -I$(top_srcdir)/src/ida -I$(top_builddir)/include

LIB_REVISION = 1:0:0

FIDA_LIB = libsundials_fida.la

FIDA_SRC_FILES = fida.c fidaband.c fidadense.c fidajtimes.c fidapreco.c fidaewt.c fidaroot.c fidabbd.c 
FIDA_BL_SRC_FILES = fidalapack.c fidalapdense.c fidalapband.c

FIDA_OBJ_FILES = $(FIDA_SRC_FILES:.c=.o)
FIDA_BL_OBJ_FILES = $(FIDA_BL_SRC_FILES:.c=.o)

FIDA_LIB_FILES = $(FIDA_SRC_FILES:.c=.lo)
FIDA_BL_LIB_FILES = $(FIDA_BL_SRC_FILES:.c=.lo)

mkinstalldirs = $(SHELL) $(top_srcdir)/config/mkinstalldirs

# ----------------------------------------------------------------------------------------------------------------------

all: $(FIDA_LIB)

$(FIDA_LIB): $(FIDA_LIB_FILES)
	@if test "X${LAPACK_ENABLED}" = "Xyes"; then \
	  make lib_with_bl;                          \
	else                                         \
	  make lib_without_bl;                       \
	fi

lib_without_bl: $(FIDA_LIB_FILES)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) -o $(FIDA_LIB) $(FIDA_LIB_FILES) -rpath $(libdir) $(LDFLAGS) $(LIBS) -static -version-info $(LIB_REVISION)

lib_with_bl: $(FIDA_LIB_FILES) $(FIDA_BL_LIB_FILES)
	$(LIBTOOL) --mode=link $(CC) $(CFLAGS) -o $(FIDA_LIB) $(FIDA_LIB_FILES) $(FIDA_BL_LIB_FILES) -rpath $(libdir) $(LDFLAGS) $(LIBS) -static -version-info $(LIB_REVISION)

install: $(FIDA_LIB)
	$(mkinstalldirs) $(libdir)
	$(LIBTOOL) --mode=install $(INSTALL_LIB) $(FIDA_LIB) $(libdir)

uninstall:
	$(LIBTOOL) --mode=uninstall rm -f $(libdir)/$(FIDA_LIB)

clean:
	$(LIBTOOL) --mode=clean rm -f $(FIDA_LIB)
	rm -f $(FIDA_LIB_FILES)
	rm -f $(FIDA_BL_LIB_FILES)
	rm -f $(FIDA_OBJ_FILES)
	rm -f $(FIDA_BL_OBJ_FILES)

distclean: clean
	rm -f Makefile

fida.lo: $(srcdir)/fida.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fida.c
fidaewt.lo: $(srcdir)/fidaewt.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidaewt.c
fidaband.lo: $(srcdir)/fidaband.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidaband.c
fidadense.lo: $(srcdir)/fidadense.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidadense.c
fidalapack.lo: $(srcdir)/fidalapack.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidalapack.c
fidalapband.lo: $(srcdir)/fidalapband.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidalapband.c
fidalapdense.lo: $(srcdir)/fidalapdense.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidalapdense.c
fidajtimes.lo: $(srcdir)/fidajtimes.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidajtimes.c
fidapreco.lo: $(srcdir)/fidapreco.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidapreco.c
fidabbd.lo: $(srcdir)/fidabbd.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidabbd.c
fidaroot.lo: $(srcdir)/fidaroot.c
	$(LIBTOOL) --mode=compile $(CC) $(CPPFLAGS) $(INCLUDES) $(CFLAGS) -c $(srcdir)/fidaroot.c

libtool: $(top_builddir)/$(LIBTOOL_DEPS)
	@cd ${top_builddir} ;                 \
	 ${SHELL} ./config.status --recheck ; \
	 cd ${abs_builddir}
