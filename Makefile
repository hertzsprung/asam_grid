MAKEFLAGS += --no-builtin-rules                                                                 
.DEFAULT_GOAL := SRC/GridGen
.DELETE_ON_ERROR:
.SUFFIXES:
.PHONY: install

INSTALL := cp
INSTALL_PROGRAM := $(INSTALL)
prefix := /usr/local
exec_prefix := $(prefix)
bindir := $(exec_prefix)/bin

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
