### Makefile.am
include $(top_builddir)/version.mk

AM_CFLAGS = $(EXTRA_CFLAGS)
AM_CPPFLAGS = -D_POSIX_C_SOURCE=201001L -D_XOPEN_SOURCE=700 -D_BSD_SOURCE
AM_CPPFLAGS += -DHAVE_VERSION_H
AM_LDFLAGS = $(XCCLDFLAGS)

bin_PROGRAMS =
noinst_PROGRAMS =
lib_LTLIBRARIES =
noinst_LTLIBRARIES =
pkglib_LTLIBRARIES =
noinst_HEADERS =
BUILT_SOURCES =
EXTRA_DIST = $(BUILT_SOURCES)
CLEANFILES = 
SUFFIXES =

CLEANFILES += version.c
EXTRA_DIST += version.c.in

noinst_LIBRARIES = libdfp.a
libdfp_a_SOURCES =
libdfp_a_SOURCES += dfp754_d32.c dfp754_d32.h
libdfp_a_SOURCES += dfp754_d64.c dfp754_d64.h
libdfp_a_CPPFLAGS = $(AM_CPPFLAGS)
libdfp_a_CPPFLAGS += $(dfp754_CFLAGS)

bin_PROGRAMS += mcfit
mcfit_SOURCES = mcfit.c mcfit.yuck
mcfit_SOURCES += version.c version.h
mcfit_CPPFLAGS = $(AM_CPPFLAGS)
mcfit_CPPFLAGS += $(dfp754_CFLAGS)
mcfit_CPPFLAGS += -DHAVE_VERSION_H
mcfit_LDFLAGS = $(AM_LDFLAGS)
mcfit_LDFLAGS += $(dfp754_LIBS)
mcfit_LDADD = libdfp.a
BUILT_SOURCES += mcfit.yucc

bin_PROGRAMS += mcfit2
mcfit2_SOURCES = mcfit2.c mcfit2.yuck
mcfit2_SOURCES += version.c version.h
mcfit2_CPPFLAGS = $(AM_CPPFLAGS)
mcfit2_CPPFLAGS += $(dfp754_CFLAGS)
mcfit2_CPPFLAGS += -DHAVE_VERSION_H
mcfit2_LDFLAGS = $(AM_LDFLAGS)
mcfit2_LDFLAGS += $(dfp754_LIBS)
mcfit2_LDADD = libdfp.a
BUILT_SOURCES += mcfit2.yucc

bin_PROGRAMS += mcmc
mcmc_SOURCES = mcmc.c mcmc.yuck
mcmc_SOURCES += version.c version.h
mcmc_SOURCES += pcg_basic.c pcg_basic.h
mcmc_CPPFLAGS = $(AM_CPPFLAGS)
mcmc_CPPFLAGS += $(dfp754_CFLAGS)
mcmc_CPPFLAGS += -DHAVE_VERSION_H
mcmc_LDFLAGS = $(AM_LDFLAGS)
mcmc_LDFLAGS += $(dfp754_LIBS)
mcmc_LDADD = libdfp.a
BUILT_SOURCES += mcmc.yucc

bin_PROGRAMS += mcmc2
mcmc2_SOURCES = mcmc2.c mcmc2.yuck
mcmc2_SOURCES += version.c version.h
mcmc2_SOURCES += pcg_basic.c pcg_basic.h
mcmc2_CPPFLAGS = $(AM_CPPFLAGS)
mcmc2_CPPFLAGS += $(dfp754_CFLAGS)
mcmc2_CPPFLAGS += -DHAVE_VERSION_H
mcmc2_LDFLAGS = $(AM_LDFLAGS)
mcmc2_LDFLAGS += $(dfp754_LIBS)
mcmc2_LDADD = libdfp.a
BUILT_SOURCES += mcmc2.yucc


## version rules
version.c: version.c.in $(top_builddir)/.version
	$(AM_V_GEN) PATH="$(top_builddir)/build-aux:$${PATH}" \
		yuck$(EXEEXT) scmver --ignore-noscm --force -o $@ \
			--use-reference --reference $(top_builddir)/.version $<

## yuck rule
SUFFIXES += .yuck
SUFFIXES += .yucc
.yuck.yucc:
	$(AM_V_GEN) PATH=$(top_builddir)/build-aux:"$${PATH}" \
		yuck$(EXEEXT) gen -o $@ $<

## Create preprocessor output (debugging purposes only)
.c.i:
	$(COMPILE) -E -o $@ $<

## Create assembler output (debugging purposes only)
.c.s:
	$(COMPILE) -S -c $(AM_CFLAGS) $<
