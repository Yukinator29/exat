#!/usr/bin/env python

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
# buildmatrix.py
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

import sys, os
import argparse as arg
import numpy    as np
import common   as c


# *****************************************************************************

def matrixbuilder(site,coup):
  ntran = c.NTran

# Check the dimensions
  dimen = sum(ntran)
  ncoup = 0
  for i in range(c.NChrom):
    for j in range(i+1,c.NChrom):
      ncoup += ntran[i]*ntran[j]

  lcoup = len(coup)

  if lcoup != ncoup :
    print("Couplings Found     : %4d" % lcoup)  
    print("Couplings Requested : %4d" % ncoup)  
    c.error("Confused in the Dimension!","matrixbuilder")

  if c.v(1):
    print(" ... Matrix dimension       : %4d" % dimen) 
    print(" ... Number of chromophores : %4d" % c.NChrom) 
    print(" ... Number of COUPLINGS    : %4d" % ncoup) 

# Convert site energies in cm-1
  site = site*c.PhyCon['eV2wn']

# Set all matrix elements to zero:
  mat = np.zeros((dimen,dimen))

# Write the diagonal elements:
  for i in range(dimen):
    mat[i][i] = site[i]/2
 
# Write the off-diagonal blocks:
  if c.OPT['CleanCoup'] > 0.0:  
    if c.v(): 
      print("   ... applying treshold of %4.1f" % c.OPT['CleanCoup'])
  if c.OPT['ScaleCoup'] != 1.0: 
    if c.v(): 
      print("   ... applying scaling factor of %4.1f" % c.OPT['ScaleCoup'])
  L=0
  for igi in range(c.NChrom):
    for igj in range(igi+1,c.NChrom):
      for i in range(ntran[igi]):
        for j in range(ntran[igj]):
          I = sum(ntran[0:igi])+i
          J = sum(ntran[0:igj])+j
          if (coup[L] >= c.OPT['CleanCoup']) or (coup[L] <= -c.OPT['CleanCoup']): 
            mat[I][J] = coup[L]*c.OPT['ScaleCoup']
          else:
            mat[I][J] = 0.0
          L = L+1

# Build the whole matrix
  M1 = np.matrix(mat)
  M2 = np.matrix.transpose(M1)
  M=M1+M2


# Save the matrix into a file
  np.savetxt(c.OutFiles['matrix'],M,fmt='%10.1f',delimiter='',newline='\n')
# print("   ... File "+c.OPT['OutMatrix']+" has been saved!")
#  
  return M

# *****************************************************************************

