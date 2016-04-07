MAKEFLAGS += --no-builtin-rules                                                                 
.DEFAULT_GOAL := SRC/GridGen
.DELETE_ON_ERROR:
.SUFFIXES:
.PHONY: install deb

INSTALL := cp
INSTALL_PROGRAM := $(INSTALL)
DEBUILD := debuild
prefix := /usr/local
exec_prefix := $(prefix)
bindir := $(exec_prefix)/bin

deb:
	$(DEBUILD) -i -us -uc -b

install: SRC/GridGen | $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $< $(DESTDIR)$(bindir)/

$(DESTDIR)$(bindir):
	mkdir -p $@

SRC/GridGen:
	$(MAKE) -C SRC

clean:
	find . -name *.mod -exec $(RM) {} \;
	$(RM) SRC/GridGen
	$(RM) -r SRC/LIBGRID
