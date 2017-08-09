#!/bin/bash

# change this if needed
#EXAT='/home/lore/molecolab_tools/exat/src/exat.py --nchrom 3 --ntran 3'

if [ "$1" == "verbose" ]; then
  outstream='/dev/stdout'
else
  outstream='/dev/null'
fi

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

## Run all tests for gdvh23
ntests=0
nerrors=0


echo ' ==================================================== '
echo '                   EXAT-EXT TESTS                     '
echo ' ==================================================== '
echo -e  " Running test 0: "
# No seltran, no other things
$EXAT -vv Test1.in --indipo Test1.dipo.in --insite Test1.site.in\
 --incoup Test1.coup.in --incent Test1.cent.in\
 -e -o test0 | tee test0.EXAT.out > $outstream
echo
difftest 0 results.out diag.dat matrix.dat
a=$? ; nerrors=`expr $nerrors + $a`
ntests=`expr $ntests + 1`

### END 



echo
echo ' Test completed'
echo ""
printf " %3d tests completed  \n" $ntests
printf " %3d tests failed \n"  $nerrors


echo 
echo -e 'END\n'

exit 0 

