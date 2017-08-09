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
# exat.py MAIN PROGRAM
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

# Import python modules
import sys, os, glob 
import numpy as np

# Import COMMON modules
import common as c   # Common
import util   as u   # Utilities

# Import WORKING modules
import readdata      # Read data (from Gaussian output)
import buildmatrix   # Build the excitonic matrix
import diag          # Diagonalize the excitonic matrix
import trans         # Compute the excitonic dip and R


# *****************************************************************************
#
# (0) Read the options from input line and set the OPT dictionary
#    

c.welcome()
u.GetOpts()

# *****************************************************************************
#
# (1) Read the input parameters from external files provided
#     by the user or extract data from different Gaussian versions
#

print "\n (1) Read data" 

# System-independent call (including seltran)
Cent,DipoLen,DipoVel,Mag,Site,Coup,Kappa = readdata.Read()


# Modify data

# If requested reorient TrDipo and change the sign to corresponding couplings
# This option works only with 1 transition per chromophore
if ( c.OPT['reorient'] is not None ) : 
  if c.OPT['read'] == 'external': 
    c.error(' Reorient Dipoles is not possible with read external')
  else:
    print " ... The direction of transition dipole moments will be changed according to selected axis"
    if c.OPT['forcedipo'] == True : 
      print " ... The orientation of transition dipole moments will forced parallel to selected axis"
    # Reorient DipoLen and DipoVel, on the basis of DipoLen (!!)
    DipoLen,DipoVel,Mag,Coup = u.reorientdipo(DipoLen,DipoVel,Mag,Coup)


# If requested scales dipoles and couplings
if c.OPT['ScaleTran'] == True: u.scaletran(DipoLen,Coup,Mag) 

if c.OPT['read'] != 'external':
  # If requested, modify electric transition dipoles
  if c.OPT['ModDipoLen'] == True: DipoLen = u.moddipo(DipoLen)
  # If requested performs the transition dipole angle analysis
  if c.OPT['anadipo'] == True : u.dipoanalysis(DipoLen)
  # If requested, modify site energies
  if c.OPT['ModSite'] == True: u.modsite(Site)

# If requested, scale dipole (Length)
if  c.OPT['ScaleDipo'] != 1.0 : 
  print " ... All transition dipoles (length) will be scaled by factor %8.4f"\
     % c.OPT['ScaleDipo']
  DipoLen  *= c.OPT['ScaleDipo']              # Scale transition dipole moments 

# If verbosity is > 1 , print site energies and save site.out file
if c.OPT['verbosity'] > 1 : u.prtsite(Site)

# *****************************************************************************
#
# 2. Build the Excitonic Matrix
#
print "\n (2) Build the excitonic hamiltonian" 

if c.OPT['read'] != 'external' and  c.OPT['ModCoup'] == True : 
  Coup = u.modcoup(len(Coup),Coup)

M = buildmatrix.matrixbuilder(Site,Coup)

# If verbosity is > 1 , print couplings and distances and save coup.out file
if c.v(1): u.prtcoup(Coup,Cent,Kappa)

# *****************************************************************************
#
# 3. Diagonalize the Excitonic Hamiltonian
#

print "\n (3) Diagonalize the excitonic hamiltonian" 
energy,coeff = diag.diagonalize(M)

# Convert some quantities to A.U.
EEN  =  energy/c.PhyCon['Town']  # Excitonic Energies (Hartree)

# *****************************************************************************
#
# 4. Perform excitonic calculations
#

RxDel  = np.cross(Cent/c.PhyCon['ToAng'],DipoVel)

# Compute internal magnetic moment (DipoMag is gauge-dependent)
MagInt = Mag - RxDel

print "\n (4) Compute the excitonic properties" 

EXCDipoLen = trans.EXCalc(coeff,DipoLen)

# Compute Linear Absorption Spectrum
print "\n ... Compute the Linear Absorption Spectrum" 
EXCDipo2   = np.sum(EXCDipoLen**2,axis=1)

# Compute Linear Dichroism Spectrum
print "\n ... Compute the Linear Dichroism Spectrum" 
LD = trans.LinDichro(EXCDipoLen)

# Compute Rotational Strength ...
print "\n ... Compute the Circular Dichroism Spectrum" 
EXCRot = trans.RotStrength(EEN,Cent,coeff,DipoLen,EXCDipoLen,
                             DipoVel,MagInt,RxDel,Site)

#
# 5. Print out the reuslts
#

# Visualization options

# Save file for DIPOLE visualization in VMD
# This is here to print also excitonic dipo
# Print intrinsic magnetic moments
if c.OPT['read'] != 'external': 
  u.savegeom()
  u.savevisudipo(Cent,DipoLen,EXCDipoLen,-MagInt)
 #np.savetxt(c.OutFiles['cent'],Cent,fmt='%12.5f',
 #  header='Coordinates of centers (ang)')

# Print out results.out
u.resout(energy,EXCDipo2,LD,EXCRot)

print("\n Done! \n")

