SHELL=/bin/bash
export SPECTRUM := ../../src/spectrum.py

test: test.spectrum

test.spectrum: clean spectrum_executable
	./runtest_spectrum.sh

test.plot: clean spectrum_executable
	./runtest_spectrum.sh plot

####

clean:
	rm -f test*

spectrum_executable: 
	@(if  [ ! -x "${SPECTRUM}" ] ; then \
	    echo "ERROR: could not find exat.py"; exit 2; fi;\
	) 	
