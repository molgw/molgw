# This file is part of MOLGW
# Author: Fabien Bruneval

-include ./src/my_machine.arch
-include ./my_machine.arch

PREFIX?=.

.PHONY: test clean archive tarball archive install

molgw:
	cd src && $(MAKE)

test:
	cd tests && python3 ./run_testsuite.py

clean:
	cd src && $(MAKE) clean

tarball:
	cd src && $(MAKE) tarball

archive:
	cd src && $(MAKE) archive

install: molgw
	mkdir -p $(PREFIX)/bin
	cp -u molgw $(PREFIX)/bin/molgw
	cp -rp basis $(PREFIX)/basis

#uninstall:
#	$(RM) -r $(PREFIX)

