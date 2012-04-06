install: R.sh
	@echo creating front-end shell script
	@cat R.sh | (cd ../.. ; sed "s@RHOME=.*@RHOME=`pwd`@" > bin/R)
	@chmod +x ../../bin/R
