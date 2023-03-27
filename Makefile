# This file is part of MOLGW
# Author: Fabien Bruneval

-include ./my_machine.arch

.PHONY: all clean archive tarball archive

all:
	cd src && $(MAKE)

test:
	cd tests && python ./run_testsuite.py

clean:
	cd src && $(MAKE) clean

tarball:
	cd src && $(MAKE) tarball

archive:
	cd src && $(MAKE) archive


