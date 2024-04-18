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
# common.py COMMON MODULE
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
import sys, os 
import numpy as np

# VERSION
VERSION  = "1.0.1"
PROGVERS = "Exat - EXcitonic Analysis Tool - Version %s" % VERSION

# ******************************************************************************
#
# Fill OPT dictionary with default values. If you want to change the default
# options you can modify here.
#
# ******************************************************************************

OPT =    {
            'workdir'   : os.getcwd(),      # Specify the working directory
            'verbosity' : 0 ,               # Verbosity of the output on screen
            'RCalc'     : 'mu',             # Select type of calculation for R
            'Cent'      : 'geom',           # How to compute the center of the chromophore
            'external'  : False,            # Read all data from external files
            'read'      : 'g16',            # Select Gaussian version for EET coupling 
            'logfile'   : None,             # Gaussian LogFile
            'CleanCoup' : 0.0,              # Treshold for couplings
            'ScaleCoup' : 1.0,              # Scaling factor for coupling
            'ScaleDipo' : 1.0,              # Scaling factor for dipoles
            'reorient'  : None,             # Axis to reorient transition dipole moments (length)
            'forcedipo' : False,            # Force dipo to be aligned with --reorient axis
            'anadipo'   : None,             # Request dipole orientation analysis
            'ScaleTran' : None,             # Request dipole and coupling scaling
            'ModSite'   : False,            # Whether to modify site energies
            'ModDipoLen': False,            # Whether to modify transition electric dipoles
            'ModCoup'   : False,            # Whether to modify couplings
            'ModCent'   : False,            # Whether to modify centers
            'OutPrefix' : False,
            'LDAxis'    : "z",              # LD Axis
            'ctype'     : 'geom' ,          # How define the center of the chromophore
            'coup'      : 'trden',          # Type of coupling calculation
            'refrind'   : 1.0,              # Refraction index for computing screening in point-dipole coupling calculation
            'seltran'   : False,            # Select only particular transitions.
            'CouPrtThr' : False,            # Threshold for coupling printout
            'savetprod' : False,            # Save triple product of (mu) CD
            'savecoeff' : False             # Save hi-acc coefficients in numpy file
          }

# ******************************************************************************
#
# Physical constants and conversion factors
#
# ******************************************************************************

PhyCon =   {
            'Planck'    : 6.62606957E-34  , # Plank constant ( J * sec)
            'Slight'    : 2.99792458E10   , # Speed of ligth (m/s)
            'Avog'      : 6.02214129E23   , # Avogadro number
            'ToAng'     : 0.52917721092   , # Angstroms per Bohr
            'ToKG'      : 1.660538921E-27 , # Kilograms per AMU
            'ToE'       : 1.602176565E-19 , # ESU per electron charge (Coloumb*electron)
            'Town'      : 220000          , # cm-1 per hartree
            'ToDeb'     : 2.54158         , # Debye per electroncharge per Bohr
            'Hartre'    : 4.35974434E-18  , # Joules per Hartre
            'eV2wn'     : 8065.5446811132 , # eV per wavenumber 
            'ToeV'      : 27.211396132    , # eV per Hartree
            'EMKG'      : 4.35974434E-18*10000  # ????
           }


# ******************************************************************************
#
# Atomic Masses
#
# ******************************************************************************

AtMass =   {   1  : 1.00797  ,
               5  : 11.0093053,
               6  : 12.01115 ,
               7  : 14.0067  ,
               8  : 15.9994  ,
               9  : 18.9984033,
              12  : 24.312  }


# ******************************************************************************
#
# Standard file name
#
# ******************************************************************************

ExtFiles = {
   # Input Files
   'incoup'       : 'coup.in',        # Coupling values (cm-1)
   'insite'       : 'site.in',        # Site energies (eV)
   'incent'       : 'cent.in',        # Postion of the centers (Ang)
   'indipo'       : 'dipo.in',        # Electric transtion dipoles (lenght formulation) in a.u.
   'indipovel'    : 'dipovel.in',     # Electric transtion dipoles (velocity formulation) in a.u.
   'indipomag'    : 'dipomag.in',     # Magnetic transtion dipoles in a.u.
   'crlist'       : 'chromlist.in',   # List of chromophores and selected tranitions
   'refaxis'      : 'reference.in',   # List of chromophores, refaxes and angles (for select tr)
   'ScaleTran'    : 'scaletran.in'    # List of chromophores and scaling factors (for select tr)
   }   

OutFiles = {
   # Output Files
   'results'      : 'results.out',    # Excitonic energies, mu^2,LD,CD
   'matrix'       : 'matrix.dat',     # Excitonic matrix 
   'diag'         : 'diag.dat',       # Excitonic coefficients and squared coefficients
   'coeff'        : 'coeff.npy',      # Excitonic coefficients (numpy file) 
   'tprod'        : 'tprod.dat',      # Triple product matrix 
   'xyz'          : 'geometry.xyz',   # Complete geometry
   'visudipo'     : 'visudipo.vmd',   # VMD script to visualize tr dipoles
   'dipo'         : 'dipo.out',       # Site dipoles
   'magdipo'      : 'magdipo.out',    # Magnetic dipoles
   'rstrength'    : 'component.out',  # Rotational strength components
   'coup'         : 'coup.out',       # List of couplings
   'site'         : 'site.out',       # List of site energies
   'exatdata'     : 'exat.npz',       # Exat data in npz binary format
   # What?
   'dipolen'      : 'dipolen.out',    # List of transition dipole moments (length)
   'dipovel'      : 'dipovel.out',    # List of transition dipole moments (velocity)
   'cent'         : 'cent.out',
   }


# ******************************************************************************
#
# Common variables 
#
# ******************************************************************************

# These may or may not have a duplicate in the global namespace
NChrom = 2
NTran  = []
ChromList = []
CoupList = [] # List of indices k->ij (i,j chromophores)
NAtom = []    # List of atom numbers, for each chromophore
anum  = []    # List of atom numbers
xyz   = []    # Geometry

# ******************************************************************************
#
# Print the welcome message
#
# ******************************************************************************

def welcome():
   print('                                                                      ')
   print('                                        mm                            ')
   print('                                     mMMm                             ')
   print('                                   mMMMMm         m                   ')
   print('                                  mMMMMm          mMm                 ')
   print('                                  mMMMMm          mMm                 ')
   print('                                  mMMMMMm        mMMm                 ')
   print('                                  MMMMMMMMMMMMMMMMMMm                 ')
   print('                                 mMMMMMMMMMMMMMMMMMm                  ')
   print('        __  ___      __    ____________      __MMMm     __            ')
   print('       /  |/  /___  / /   / ____/ ____/___  / /  ____ _/ /_           ')
   print('      / /|_/ / __ \/ /   / __/ / /   / __ \/ /  / __ `/ __ \          ')
   print('     / /  / / /_/ / /___/ /___/ /___/ /_/ / /__/ /_/ / /_/ /          ')
   print('    /_/  /_/\__________/_____/\____/_____/_____|__,_/_.___/           ')
   print('            /_  __/ __ \/ __ \/ /  / ___/                             ')
   print('             / / / / / / / / / /   \__ \                              ')
   print('            / / / /_/ / /_/ / /___ __/ /                              ')
   print('           /_/  \____/\____/_____/____/                               ')
   print('             mMMMMMMMMMMMMMMMm                                        ')
   print('           mMMMMMMMMMMMMMMMm                                          ')
   print('         mMMMMMMMMMMMMMMMMM   + ------------------------------------ +')
   print('        mMMMMMMMMMMMMMMMMm    |             E  X  A  T               |')
   print('       mMMMMMMMMMMMMMMMMMm    + ------------------------------------ +')
   print('       mMMMMm       mMMMMMm   | S. Jurinovich, L. Cupellini          |')
   print('       mMMMm       mMMMMMMm   | C.A. Guido                           |')
   print('        mMm       mMMMMMMm    |                            ver %-6.5s|' % VERSION)
   print('         m       mMMMMMMm     |              molecolab.dcci.unipi.it |')
   print('                mMMMMMm       + ------------------------------------ +')
   print('                                                                      ')


# ******************************************************************************
#
# Standard error:
#
# ******************************************************************************

def error(string,where=None):
  # Determine calling function
  if where is None: where = sys._getframe().f_back.f_code.co_name
  print(("\n------ ERROR in %s ------\n%s\n" % (where,string)))
  sys.exit()


# ******************************************************************************
#
# Debug message:
#
# ******************************************************************************

def debug():
  # Determine calling function
  where = sys._getframe().f_back.f_code.co_name
  print(("\n------ EXIT FOR DEBUGGING -------\nI am in: %s\n\n" % where))
  sys.exit()


# ******************************************************************************
#
# Progression Bar
#
# ******************************************************************************

def ShowProgressBar(istep,ntot):
  percentage = np.round(float(istep)/float(ntot)*100)
  sys.stdout.write('\r')
  sys.stdout.write("[%-50s] %d%%" % ('='*int(percentage/2), percentage))
  sys.stdout.flush()


# ******************************************************************************
#
# Print lists by using a specifc format
#
#   List   : input list to print on screen (list of elements of the same type)
#   NER    : maximum number of element to be print in a row (integer)
#   fmt    : format specification (string)
#   margin : (optional) number of white spaces from left margin (integer)
#   header : (optional) text header (string)
#
# ******************************************************************************

def printout(List,NER,fmt,margin=2,header=''):
  try:
    N     = len(List)          # Total number of elements in List
  except:
    N     = 1
  if N > 1:
    NER   = int(NER)           # Number of element per row
    NFR   = N/NER              # Number of full row
    NLR   = N-NFR*NER          # Number of element in the last row
    if ( N >= NER):
      Fmt = (" "*margin+(fmt)*NER+"\n")*(N/NER)+(" "*margin+(fmt)*NLR+"\n")
    else:
      Fmt = (" "*margin+(fmt)*N+"\n")
    print("\n"+" "*margin+header+"\n")
    print(Fmt % tuple(List))
  else:
    print("\n"+" "*margin+header+"\n")
    print(fmt % List)
  pass


# ******************************************************************************
#
# StringConverter
#
# ******************************************************************************
def stringconverter(In):

  Out = []
  blocks = In.split(',')
  for block in blocks:
    blk = block.split('-')
    if len(blk) > 1:
      sel =  list(map(int,blk))
      Out += list(range(sel[0],sel[1]+1,1))
    else:
      Out += [int(block)]
  return Out


# ******************************************************************************
#
# Check if the file exists:
#
# ******************************************************************************

def checkfile(FILENAME):
  if (os.path.isfile(FILENAME) == False):
    print(("\n File %s not found!\n" % FILENAME))
    sys.exit()

# ******************************************************************************
#
# Set prefix to output files
#
# ******************************************************************************
def setoutfiles():
  if OPT['OutPrefix'] is not False: 
    for k in OutFiles:
      OutFiles[k] = OPT['OutPrefix']+'.'+OutFiles[k]
  
# ******************************************************************************
#
# Print out a dictionary or a list (OPT, etc) 
#
# ******************************************************************************

def PrintHeader(title=None):
  length = 70
  sep    = " "+length*'-'
  if title is not None:
    # Find out where to put the title 
    fblank = length/2 - len(title)/2
    fmt     = sep+'\n'+' '*fblank+'%s\n'+sep
    print(fmt % (title))
  else:  print(sep)
  pass

def PrintDict(Stuff,title):
  PrintHeader(title)
  if type(Stuff) is dict:
    for key,val in list(Stuff.items()):  print("   %-14s : %-10s" % (key,val))
  elif type(Stuff) is list:
    for key,val in enumerate(Stuff):  print("   %14d : %-10s" % (key,val))
  pass
  PrintHeader(None)

def PrintOpt():
  PrintDict(OPT,'Selected options:')
  pass

# Return if to print the message!
def v(level=-1):
  v = OPT['verbosity']
  if v is None: v = 0
  return v > level

# ******************************************************************************
#
# VMD Header
#
# ******************************************************************************
def VMDHeader():
  s  = '## Generated with %s\n\n' % PROGVERS 
  s += '## VMD script to view transition dipoles... follow the instructions\n'
  s += '## below the VMD functions to select the dipole(s) to draw         \n\n\n'
  s += '## General Settings                                                \n'
  s += 'menu main on                                                       \n'
  s += 'display projection orthographic                                    \n'
  s += '                                                                   \n'
  s += '## VMD functions to draw a vector                                  \n'
  s += 'proc vmd_draw_arrow {mol start end} {                              \n'
  s += '  set length [veclength [vecsub $end $start]]                      \n'   
  s += '  set conelen [expr max(0.4,0.2*$length) ]                         \n'   
  s += '  set scale [expr max(0.5,(1.0-$conelen/$length))]                 \n'   
  s += '                                                                   \n'   
  s += '  set middle [vecadd $start [vecscale $scale [vecsub $end $start]]]\n'   
  s += '  # You can change the radius of the cylinder/cone                 \n'   
  s += '  graphics $mol cylinder $start $middle radius 0.05 resolution 24  \n'   
  s += '  puts [list cone $middle $end radius 0.15]                        \n'   
  s += '  graphics $mol cone $middle $end radius 0.15 resolution 24        \n'   
  s += '}                                                                  \n'
  s += '                                                                   \n'
  s += 'proc vmd_draw_vector { mol pos val } {                             \n'
  s += '    # Change +1 with any value to scale all vectors!               \n'   
  s += '    set end   [ vecadd $pos [ vecscale +1 $val ] ]                 \n'
  s += '    vmd_draw_arrow $mol $pos $end                                  \n'
  s += '}                                                                  \n\n\n'
  # Instructions
  s += '## Under the comment CHROM X TRAN J there is the command to draw the transition dipole.  \n'
  s += '## The first command prints the electric dipole, the second prints the magnetic dipole.  \n'
  s += '## The electric first transition of each chromophore is displayed by default\n '
  s += '## Below the site transition dipoles, there are also the excitonic (electric) dipoles\n '
  s += '## Uncomment the ones you want to display before running vmd.  \n\n\n'
  return s



# --- END of MODULE ---
