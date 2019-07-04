#AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../src -I$(srcdir)/../@BIOMCMCLIB@/lib \
#	-I$(srcdir)/../submodules/parasail/ @OPENMP_CPPFLAGS@ -DTEST_FILE_DIR='"$(srcdir)/files/"'
#LDADD = $(GTKDEPS_LIBS)  @CHECK_LIBS@   ../submodules/parasail/libparasail.la  ../src/libamburana.la ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)
AM_CPPFLAGS = $(GTKDEPS_CFLAGS) -I$(srcdir)/../src -I$(srcdir)/../@BIOMCMCLIB@/lib @OPENMP_CPPFLAGS@ -DTEST_FILE_DIR='"$(srcdir)/files/"'
LDADD = $(GTKDEPS_LIBS)  @CHECK_LIBS@   ../src/libamburana.la ../biomcmc-lib/lib/libbiomcmc.la  $(AM_LDFLAGS)

AM_CFLAGS = @AM_CFLAGS@  @OPENMP_CFLAGS@  @CHECK_CFLAGS@


EXTRA_DIST = files # directory with fasta etc files (accessed with #define TEST_FILE_DIR above)
LIST_OF_TEST_PROGS= check_unit check_simpleheap
# duplicate, since we want all to be compiled only with 'make check'
TESTS = $(LIST_OF_TEST_PROGS)           # list of test programs
check_PROGRAMS = $(LIST_OF_TEST_PROGS)  # list of programs to be compiled only with 'make check'

check_unit_SOURCES = check_unit.c
check_simpleheap_SOURCES = check_simpleheap.c
