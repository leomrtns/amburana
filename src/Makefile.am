#AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../@BIOMCMCLIB@/lib -I$(srcdir)/../submodules/parasail/ @OPENMP_CPPFLAGS@
#LDADD = $(GTKDEPS_LIBS) ../submodules/parasail/libparasail.la ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)
AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@
LDADD = $(GTKDEPS_LIBS)  ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
## TODO: libparasail should not be hardcoded here, should be set by configure

common_headers = amburana.h simple_heap.h general_heap.h minhash.h
common_src     = simple_heap.c minhash.c
noinst_LTLIBRARIES = libamburana.la   ## change to noinst_LT... if should be linked statically (not installed globally)
libamburana_la_SOURCES = config.h $(common_headers) $(common_src) 
bin_PROGRAMS = amburana 
amburana_SOURCES = main.c config.h 
