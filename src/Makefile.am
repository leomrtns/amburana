AM_CPPFLAGS = $(GTKDEPS_CFLAGS)  -I$(srcdir)/../@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@
LDADD = $(GTKDEPS_LIBS)  @ZLIB_LIBS@ ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@  @ZLIB_CFLAGS@

common_headers = amburana.h simple_heap.h general_heap.h minhash.h sketch_distance.h clustering_affineprop.h \
                 clustering_hierarchical.h
common_src     = simple_heap.c minhash.c sketch_distance.c clustering_affineprop.c \
                 clustering_hierarchical.c

noinst_LTLIBRARIES = libamburana.la   ## change to noinst_LT... if should be linked statically (not installed globally)
libamburana_la_SOURCES = config.h $(common_headers) $(common_src) 

# next lines seem redundant, but needed since I wanna use the dynamic library (.la) instead of including all "files.h"
bin_PROGRAMS = amburana 
amburana_SOURCES = main.c config.h 
amburana_LDADD = libamburana.la $(LDADD) 
amburana_CPPFLAGS = -I$(srcdir)/. $(AM_CPPFLAGS)
