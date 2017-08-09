#!/bin/bash

# change this if needed
#SPECTRUM='../../src/spectrum.py'

tstplot=$1

# Diffing function
#
diffing()
{
 tst=$1
 ref=$2
 # Always diff test ref (ignoring # comments) 
 diff -I '^#' $1 $2 > diff_${1}_${2}
 result=$?
 if test $result -eq 0; then
   rm diff_${1}_${2}
 fi
 echo $result
}

# Test function 
difftest()
{
  ntest=$1; shift
  nerr=0
  diffs=()
  echo ' ==================================================== '
  echo " Diffing results for test ${ntest}: "
  for arg in $@; do
    tstfile="test${ntest}.${arg}"
    reffile="ref${ntest}.${arg}"
    diffile="diff_${tstfile}_${reffile}"
    echo " ... diffing $tstfile $reffile"
    res=`diffing $tstfile $reffile`
    if test $res -ne 0; then
      nerr=`expr $nerr + 1`
      diffs="$diffs \n $diffile" 
    fi
  done
  # collect errors
  if test $nerr -gt 0; then 
    echo ' POSSIBLE FAILURE ... check these files: ' 
    echo -e $diffs
    echo
    result=1
  else
    echo ' PASSED'
    result=0
  fi
  echo ' ==================================================== '
# sleep 1
  return $result
}



##MAIN

## Run all tests for spectrum.py
ntests=0
nerrors=0


echo ' ==================================================== '
echo '                  SPECTRUM TESTS                      '
echo ' ==================================================== '

echo " Running test 1: "
# Normal
$SPECTRUM Test1.results.out --HWHH 300 -o test1 > test1.SPEC.out 
echo
difftest 1 OD.dat LD.dat CD.dat 
a=$? ; nerrors=`expr $nerrors + $a`
ntests=`expr $ntests + 1`

echo " Running test 2:"
# Autorange
$SPECTRUM Test1.results.out --auto  --sigma 300 -o test2 > test2.SPEC.out
echo
difftest 2 OD.dat LD.dat CD.dat 
a=$? ; nerrors=`expr $nerrors + $a`
ntests=`expr $ntests + 1`


echo
echo ' Test completed'
echo ""
printf " %3d tests completed  \n" $ntests
printf " %3d tests failed \n"  $nerrors

echo 
echo -e 'END\n'

# Plot
if [  "$tstplot" == "plot" ]; then
 echo
 echo " Test 1 try with matplotlib plotting: "
 $SPECTRUM ref1.results.out --HWHH 300 -o test_spec_plot --plot --min 300 --max 180
 rm -f test_spec_plot.*
fi

exit 0

