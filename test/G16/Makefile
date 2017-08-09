SHELL=/bin/bash
export EXAT := ../../src/exat.py
export SPECTRUM := ../../src/spectrum.py

test: clean exat_executable test.exat

test.exat: clean exat_executable 
	./runtest_exat.sh

test.verbose: clean exat_executable
	./runtest_exat.sh verbose

####

clean:
	rm -f test*

exat_executable:  
	@(if  [ ! -x "${EXAT}" ] ; then \
	    echo "ERROR: could not find exat.py"; exit 2; fi;\
	) 	

spectrum_executable: 
	@(if  [ ! -x "${SPECTRUM}" ] ; then \
	    echo "ERROR: could not find exat.py"; exit 2; fi;\
	) 	
