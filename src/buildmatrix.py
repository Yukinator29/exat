#!/usr/bin/env python
#
# Copyright 2014, Sandro Jurinovich
# This program is distributed under General Public License v. 3.  
# See the file LICENCE for a copy of the license.  
#
# > EXAT vs. 2.0
# > Module: buildmatrix.py
#
# This module builds the excitonic matrix using site energies and couplings.
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
    print "Couplings Found     : %4d" % lcoup  
    print "Couplings Requested : %4d" % ncoup  
    c.error("Confused in the Dimension!","matrixbuilder")

  if c.v(1):
    print  " ... Matrix dimension       : %4d" % dimen 
    print  " ... Number of chromophores : %4d" % c.NChrom 
    print  " ... Number of COUPLINGS    : %4d" % ncoup 

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
      print "   ... applying treshold of %4.1f" % c.OPT['CleanCoup']
  if c.OPT['ScaleCoup'] != 1.0: 
    if c.v(): 
      print "   ... applying scaling factor of %4.1f" % c.OPT['ScaleCoup']
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

#
# Running matrix builder from command line
#
if __name__ == "__main__":

  # Set default options 
  OPT = { 'verbosity' : 2,
          'external'  : True,
          'OutMatrix' : 'matrix.dat',
          'CleanCoup' : 0.0}

  # Parse inpute line
  parser = arg.ArgumentParser(description="Build the excitonic matrix")
  parser.add_argument('-v',help='Increase the verbosity of the output',action="count")
  parser.add_argument('-s',help='File containing the site energies, default is site.in',default="site.in")
  parser.add_argument('-c',help='File containing the couplings, default is coup.in',default="coup.in")
  parser.add_argument('-o',help='Output file where excitonic matrix will be saved, default is matrix.dat',default="matrix.dat")
  parser.add_argument('-nt',help='Number of transition per chromophores (supposed to be equal for each chromophore)',type=int)
  parser.add_argument('-nc',help='Number of chromophores)',type=int)
  parser.add_argument('--sunit',help='Specify the unit for the site energies',default='eV',choices=['eV','cm-1'])
  parser.add_argument('--cunit',help='Specify the unit for the couplings',default='cm-1',choices=['eV','cm-1'])
  parser.add_argument('--tresh',help='Indicate the threshold for couplings (cm-1)',type=float,default=0.0)
  parser.add_argument('--scalesite',help='Shift the site energies in cm-1',type=float,default=0.0)
#  parser.add_argument('-p',help='Activate split triangular part modality',default=False)
#  parser.add_argument('--triup',help='Values printed in low triangular part',type=float,default=0.0)
#  parser.add_argument('--tridw',help='Values printed in low triangular part',type=float,default=0.0)
  args = parser.parse_args()

  if ( args.nc == None ) : print ("\n Specifify the number of chromophores!\n ") ; sys.exit()
  if ( args.nt == None ) : print ("\n Specifify the number of transitions!\n ")  ; sys.exit()


  # Set options
  OPT['OutMatrix'] = args.o
  if args.tresh > 0.0 : OPT['CleanCoup'] = args.tresh
  SITEFILE   = args.s
  COUPFILE   = args.c
  NTranChrom = args.nt  # We suppose to have the same number of transitions for each chromophore!
  NChrom     = args.nc

  # Print welcome message
  c.welcome(sys.argv[0])

  # Load data from external files
  c.checkfile(SITEFILE)
  c.checkfile(COUPFILE)
  site      = np.loadtxt(SITEFILE,dtype="float")+args.scalesite/c.PhyCon['eV2wn']
  coup      = np.loadtxt(COUPFILE,dtype="float")
  NSite     = len(site)
  NTran     = map(int,[NTranChrom]*NChrom)
  NTotTran  = sum(NTran)
  
  if ( NSite != NTotTran) :
    print("\nInconsitency between NTran and dimension of site energies!\n")
    sys.exit()

  print(" NSite  : %3d " % NSite)
  print(" NTran  : %s "  % str(NTran))
  print

# Call the matrixbuilder function
  matrixbuilder(site,coup,NTran,NChrom)

# END
  print(" \nDone!\n")
