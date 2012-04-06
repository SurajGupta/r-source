## These are used for 'tests' dependencies: 
R	= ../bin/R
Rbinary = ../bin/R.binary 
Rbase   = ../library/base

install-help: help-files
install-latex: latex-files
install-html: html-files
test-Examples: Examples-test/All-Ex.Rout

##-- Here I only do 'base';  
##--  with 'help', should do it for ALL libraries
##--  with 'html' / 'latex' we don't yet have anything -- ANYWAY ! ---
help-INDEX: ../help/base/AnIndex 
../help/base/AnIndex : ../src/library/base/man/*
##	install-base-Index...

help-files: help-INDEX
	@echo installing ALL manual entries
	./install-allhelp
	@touch $@

latex-files: help-INDEX
	@echo 'creating latex manual version (in dvi)'
	./install-latex
	@touch $@

html-files: help-INDEX
	@echo 'creating html manual version (in html)'
	./install-html
	@touch $@

Examples-test/All-Ex.Rout: Examples-test/All-Ex.R $(Rbinary) $(Rbase)
	@if [ -f $@.bak ]; then mv $@.bak $@.bakk ; fi
	@if [ -f $@ ]; then mv $@ $@.bak ; fi
	@echo Running all help-EXAMPLES
	(cd Examples-test;rm -f .RData; ../$(R) -v 4 < All-Ex.R)> $@

Examples-test/All-Ex.R: help-INDEX
	@if [ ! -d `dirname $@` ]; then mkdir `dirname $@` ; fi
	@if [ -f $@.bak ]; then mv $@.bak $@.bakk ; fi
	@if [ -f $@ ]; then mv $@ $@.bak ; fi
	@echo Extracting all help-EXAMPLES into $@
	./extract-Examples > $@

clean:
	@rm -f COMPILE INSTALL SHLIB
	@rm -f Examples-test/All-Ex.R

realclean: clean
	@rm -fr Examples-test
