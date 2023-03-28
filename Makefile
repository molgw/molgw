# This file is part of MOLGW
# Author: Fabien Bruneval


.PHONY: all clean archive tarball archive

all:
	cd src && $(MAKE)

test:
	cd tests && python ./run_testsuite.py

github-test:
	cd tests && python ./run_testsuite.py --exclude benzene_he_rt-tddft.in

clean:
	cd src && $(MAKE) clean

tarball:
	cd src && $(MAKE) tarball

archive:
	cd src && $(MAKE) archive


