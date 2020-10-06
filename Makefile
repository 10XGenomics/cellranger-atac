#!/bin/bash

#
# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
#
# Build cellranger-atac-public.
#

.PHONY: all clean

LIBPY=lib/python
VERSION=$(shell git describe --tags --always --dirty)

#
# Targets for development builds.
#
all: lib/bin tenkit-all louvain plsa cython
clean: louvain-clean plsa-clean clean-cython tenkit-clean

# Targets for cython builds
# Find python libs, includes, and default flags.
PYTHON_ROOT:=$(shell python-config --exec-prefix)
PYTHON_LIB_DIR:=$(dir $(lastword $(wildcard $(PYTHON_ROOT)/lib*/libpython*.so*)))
PYTHON_SITE_PKG:=$(lastword $(wildcard $(PYTHON_LIB_DIR)/python*/site-packages))
PYTHON_CFLAGS:=$(shell python-config --cflags)
PYTHON_LDFLAGS:=$(shell python-config --ldflags)

# Find stuff that needs to be cythonized.
CYTHON_SRCS=$(shell find $(PWD)/mro/atac/stages $(PWD)/lib/cellranger/lib/python -type f -name '*.pyx')
CYTHON_LIBS=$(patsubst %.pyx, %.so, $(CYTHON_SRCS))
CYTHON_BUILDPATH=$(shell pwd)/lib/cellranger/lib/cython
CYTHON_FLAGS?=--line-directives $(EXTRA_CYTHON_FLAGS)

# Prevent make from automatically deleting intermediate files.
.PRECIOUS: $(CYTHON_BUILDPATH)/%.c $(CYTHON_BUILDPATH)/%.o

.PHONY: cython

$(CYTHON_BUILDPATH)/%.c: $(PWD)/%.pyx
	mkdir -p $(@D) && cython $(CYTHON_FLAGS) -w $(<D) -o $(abspath $@) $(<F)

$(CYTHON_BUILDPATH)/%.o: $(CYTHON_BUILDPATH)/%.c
	$(CC) $(PYTHON_CFLAGS) $(CFLAGS) -g -O3 -c -fPIC -fopenmp \
	    -I$(PYTHON_SITE_PKG)/numpy/core/include \
	    -o $@ \
	    $<

$(PWD)/%.so: $(CYTHON_BUILDPATH)/%.o
	$(CC) -L$(PYTHON_LIB_DIR) $(PYTHON_LDFLAGS) $(LDFLAGS) -shared -fopenmp -fPIC $< -o $@

cython: $(CYTHON_LIBS)

clean-cython:
	rm -rf $(CYTHON_BUILDPATH)
	rm -f $(CYTHON_LIBS)

tenkit-all:
	make -C tenkit all

tenkit-clean:
	make -C tenkit clean

lib/bin:
	mkdir lib/bin

louvain:
	make -C lib/louvain
	cp lib/louvain/convert lib/bin
	cp lib/louvain/louvain lib/bin

louvain-clean:
	make -C lib/louvain clean
	rm -f lib/bin/convert
	rm -f lib/bin/louvain

plsa:
	make -C lib/plsa
	cp lib/plsa/plsa lib/bin

plsa-clean:
	make -C lib/plsa clean
	rm -f lib/bin/plsa

