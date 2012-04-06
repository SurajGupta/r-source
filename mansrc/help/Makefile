all:
	@echo Initializing help directory
	-mkdir -p ../../help; cd ../../help; rm -f `ls | sed '/CVS/d'`
	-help.format ../man
	-help.index ../man < /dev/null

clean:
	@echo Cleaning in manual/help
