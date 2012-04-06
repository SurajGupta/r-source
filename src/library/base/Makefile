base:
	@echo building system base library
	@cat funs/[a-z]* |sed '/^[ 	]*$$/d'  > ../../../library/base
#				      ^^ = TAB + Space
##-- do NOT delete all '#' lines -- this BREAKS strings containing '#' !!
##-- Could still delete all  '/^[\t ]*#/ lines, but why not leave comments ?
##	@cat $(FUNS) STARTUP |sed  '/^$$/d;/^[ 	]*#/d'  > ../../../library/base
##					      ^^ = TAB + Space
