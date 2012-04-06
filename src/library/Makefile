all:
	@cd base; make
	@cd profile; make

libs:
	@echo building application libraries
	@sh LIBSETUP

clean:
	@sh LIBCLEAN
	@rm LIBSETUP LIBCLEAN
	@echo cleaning in library
