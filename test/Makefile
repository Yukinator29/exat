SHELL=/bin/bash

test: logclean test.exat test.spectrum

verbose: test.verbose

test.verbose: logclean test.G16.verbose test.external.verbose

exat: test.exat

test.exat: logclean test.G16 test.external

spectrum: logclean test.spectrum

plot: logclean test.plot

clean:
	cd G16      && make clean
	cd ext      && make clean
	cd spectrum && make clean

test.G16:
	(cd G16 && make test) | tee -a test.log

test.G16.verbose:
	(cd G16 && make test.verbose) | tee -a test.log

test.external:
	(cd ext      && make test) | tee -a test.log

test.external.verbose:
	(cd ext      && make test.verbose) | tee -a test.log

logclean:
	rm -f test.log	

test.spectrum:
	(cd spectrum && make test) | tee -a test.log

test.plot:
	(cd spectrum && make test.plot) | tee -a test.log

wipe: clean logclean

