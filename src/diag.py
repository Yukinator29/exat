#!/usr/bin/env python
#
# Copyright 2013, Sandro Jurinovich
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  
#
__description__ = \
"""
This script compute eigenvalues and eigenvectors of the excitonic matrix.
The Eigenvalues corresponds to the Excitonic Energies and the eigenvectors
are the contributes of the old state i to the new excitonic state k.
All the values in the excitonic matrix must be expressed in cm-1.

 Excitonic Matrix (M)          eigenval                 eigenvec
-----------------------      ------------       -------------------------
     E1 V12 V13                   E1'            c1(E1') c1(E2') c1(E2')
    V12  E2 V23                   E2'            c2(E1') c2(E2') c2(E2')
    V13 V23  E3                   E3'            c3(E1') c3(E2') c3(E2')

Eigenvalues are in crescent order of energies, as well as the corresponding
eigenvectors. An external formatted file is also saved containing the
eigenvalue vector, eginvector matrix.
"""
__author__ = "Sandro Jurinovich"
__date__   = "101213"

import os, sys
import numpy as np
import common as c
import argparse as arg

# *****************************************************************************

def eighsort(matrix):
  e,v = np.linalg.eigh(matrix)
  index = e.argsort()
  e = e[index]
  v = v[:,index]
  return e,v


def diagonalize(M):

#
# Use the routine linalg.eigh to compute the eigenvalues and eigenvectors of a 
# Hermitian or symmetric matrix. The routine returns a 1-D array containing the 
# eigenvalues and a 2-D square matrix of the corresponding eigenvectors (in columns).
# > eigenval = eigenvalues array,not necessarily ordererd.
# > eigenvec = eigenvectors, the coloumn eigenvec[:,i] is the normalized eigenvector
#              corresponding to the eigenvalue eigenval[i].
#

  eigenval, eigenvec = eighsort(M)

  eigenvec = eigenvec.T

  dim = len(eigenval)

  # Compute the c2
  prob = np.array(eigenvec)**2

  # Compute the partecipation ratio
  # TO BE DONE

  # Save the output into a formatted external file

  TblCoeff = np.column_stack((eigenval/c.PhyCon['eV2wn'],eigenval,eigenvec))
  TblProb  = np.column_stack((eigenval/c.PhyCon['eV2wn'],eigenval,prob))
  outfile = c.OutFiles['diag']
  out = open(outfile,'w')
  np.savetxt(out,TblCoeff,fmt="%10.4f ",delimiter='',newline='\n')
  out.write("\n")
  np.savetxt(out,TblProb,fmt='%10.4f ',delimiter='',newline='\n')
  out.close()

  # possibly save accurate coeff file
  if c.OPT['savecoeff']:
    if c.v(0):
      print " ... excitonic coefficients will be saved to file %s"\
       %  (c.OutFiles['coeff'])
    np.save(c.OutFiles['coeff'],eigenvec)
        

  # Convert the results into array-type objects
  eigenvec = np.asarray(eigenvec).reshape((dim,dim))
  eigenval = np.asarray(eigenval).reshape(-1)


  return(eigenval,eigenvec)


# *****************************************************************************

#
# Stand-alone version
#

if __name__ == "__main__" :

  OPT = { 'eunit'   : 'eV',
          'diagout' : 'diag.dat' }

  # Read input line
  parser = arg.ArgumentParser(description="This module diagonalize the excitonic matrix")
  parser.add_argument('-v',help='Increase the verbosity of the output',action="count")
  parser.add_argument('-i',help='Matrix file containing site energies and couplings',default='matrix.dat')
  parser.add_argument('-o',help='Output file containing the eigenvalues and eigenvectors',default='diag.dat')
  parser.add_argument('--eunit',help='Specify the unit for site energies (default is cm-1)',default='cm-1',choices=['eV','cm-1'])
  args = parser.parse_args()

  # Set the options
  MatrixFile = args.i
  if args.v > 0  : OPT['verbosity'] = args.v
  OPT['eunit']   = args.eunit
  OPT['diagout'] = args.o

  # Print welcome message
  c.welcome(sys.argv[0])

  # Read the excitonic matrix
  if os.path.isfile(MatrixFile) == True:
    print(" > reading the excitonic matrix from %s" % MatrixFile )
    M = np.loadtxt(MatrixFile,dtype="float")
    print("   ... matrix dimension are: %s" % str(np.shape(M)))
    if OPT['eunit'] ==  "eV" :
      print("   ... convert site energies from eV to cm-1")
      for i in range(len(M)):
        M[i][i] = M[i][i]*c.PhyCon['eV2wn']
  else:
    print("File %s is missing!" % MatrixFile)
    print("\nABORT!\n")
    sys.exit()

  # Call the diagonalization function
  print(" > diagonalizing the %s matrix" % MatrixFile )
  diagonalize(OPT,M)
  print(" > Output file %s has been saved!" % OPT['diagout'] )

  # END
  print("\n Done! \n")

