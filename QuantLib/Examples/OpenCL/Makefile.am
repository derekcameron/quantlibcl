
AM_CPPFLAGS = -I${top_srcdir} -I${top_builddir}

if AUTO_EXAMPLES
bin_PROGRAMS = OpenCL
TESTS = OpenCL$(EXEEXT)
else
noinst_PROGRAMS = OpenCL
endif
OpenCL_SOURCES = OpenCL.cpp dcmt.cpp mcsimulation.cpp
OpenCL_LDADD = ../../ql/libQuantLib.la
OpenCL_LDFLAGS = -lpthread

EXTRA_DIST = \
	acklam.cl \
	kernels.cl \
	mersennetwister.cl \
	NVIDIA_mersennetwister.cl \
	simple.cl \
	thomaslukrng.cl \
    OpenCL.dev \
    OpenCL_vc9.vcproj \
    ReadMe.txt

.PHONY: examples check-examples

examples: OpenCL$(EXEEXT)

check-examples: examples
	./OpenCL$(EXEEXT)

dist-hook:
	mkdir -p $(distdir)/bin
	mkdir -p $(distdir)/build

