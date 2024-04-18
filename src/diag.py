#!/usr/bin/env python
#

# EEEEEE   XX    XX         AAA        TTTTTTTTTT
# EE        XX  XX         AA AA       TT  TT  TT
# EE         XXXX         AA   AA          TT
# EEEEEE      XX         AAA   AAA         TT
# EE         XXXX       AAAAAAAAAAA        TT
# EE        XX  XX     AA         AA       TT
# EEEEEE   XX    XX   AA           AA      TT   
#
# EXcitonic Analysis Tool         @MoLECoLab 
# https://molecolab.dcci.unipi.it/tools/
#

#
# *************************************
# EXAT - EXcitonic Analysis Tool
# diag.py MODULE
# *************************************
#

# Copyright (C) 2014-2017 
#   S. Jurinovich, L. Cupellini, C.A. Guido, and B. Mennucci
#
# This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found in LICENSE or at
#   <http://www.gnu.org/licenses/>.
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
      print(" ... excitonic coefficients will be saved to file %s"\
       %  (c.OutFiles['coeff']))
    np.save(c.OutFiles['coeff'],eigenvec)
        

  # Convert the results into array-type objects
  eigenvec = np.asarray(eigenvec).reshape((dim,dim))
  eigenval = np.asarray(eigenval).reshape(-1)


  return(eigenval,eigenvec)


# *****************************************************************************

