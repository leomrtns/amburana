AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../@BIOMCMCLIB@/lib  @OPENMP_CPPFLAGS@
AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@
LDADD = $(GTKDEPS_LIBS) ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)

common_headers = amburana.h 
common_src     = main.c 
bin_PROGRAMS = amburana 

super_distance_SOURCES = $(common_src) $(common_headers) config.h 
