base:
	@echo building system base library
	@cat `ls funs/[a-z]* | sed '/CVS/d; /~$$/d; /\.orig$$/d; /\.bak$$/d'` \
	 | sed '/^[ 	]*$$/d'  > ../../../library/base
#		   ^^ = TAB + Space
##-- do NOT delete all '#' lines -- this BREAKS strings containing '#' !!
##-- Could still delete all  '/^[\t ]*#/ lines, but why not leave comments ?
##	@cat $(FUNS) STARTUP |sed  '/^$$/d;/^[ 	]*#/d'  > ../../../library/base
##					      ^^ = TAB + Space
