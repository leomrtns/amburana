AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)

common_headers = simple_heap.h general_heap.h minhash.h
common_src     = main.c  simple_heap.c minhash.c
bin_PROGRAMS = amburana 

amburana_SOURCES = $(common_src) config.h 

noinst_LTLIBRARIES = libamburana.la   ## change to noinst_LT... if should be linked statically (not installed globally)
libamburana_la_SOURCES = config.h $(common_headers) $(common_src) 
