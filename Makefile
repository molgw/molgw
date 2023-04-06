# This file is part of MOLGW
# Author: Fabien Bruneval

PREFIX?=.

.PHONY: all test clean archive tarball archive install

all:
	cd src && $(MAKE)

test:
	cd tests && python3 ./run_testsuite.py

clean:
	cd src && $(MAKE) clean

tarball:
	cd src && $(MAKE) tarball

archive:
	cd src && $(MAKE) archive

install:
	cp -u molgw $(PREFIX)/molgw

