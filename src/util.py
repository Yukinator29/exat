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
# util.py UTILITIES MODULE
# *************************************
#
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
#
import os,sys,glob
import numpy  as np
from scipy.spatial import distance
import common as c
import argparse as arg


# *****************************************************************************
#
# Parse input line
#

class CustomHelpFormatter(arg.HelpFormatter):
  # 
  def _split_lines(self, text, width):
  #return all but include blanklines 
    return ['']+arg.HelpFormatter._split_lines(self, text, width)+['']

# Get Options with argparse
def GetOpts():

  # Read input line
  parser = arg.ArgumentParser(description="""EXAT - EXcitonic Analysis Tool: 
  Performs an excitonic calculation from the 
  output of Gaussian g16 log file(s). """,
  formatter_class=lambda prog: CustomHelpFormatter(prog,
   max_help_position=36,indent_increment=3))

  # Specify which chromlist file should be used
  parser.add_argument('chromlist',metavar='ChromListFile',
    help='''File containing the chromophore list (Default: chromlist.in).
    The file should contain one chromophore per line.''',
    default="chromlist.in",nargs='?')

  # Version
  parser.add_argument('--version','-V',action='version',version=c.PROGVERS,
     help="Show program's full version number, and exit")  

  # Input options
  grpinp = parser.add_argument_group(' Input options')
  read_group = grpinp.add_mutually_exclusive_group()
  read_group.add_argument('--readexternal','-e',action="store_true",
     help='Read all quantities from external files')
  grpinp.add_argument('--log','-i',metavar='LogFile',
     help='Gaussian logfile name (for g16 only)')

  # Output options
  grpout = parser.add_argument_group(' Output options')
  grpout.add_argument('--out','-o',metavar='Prefix',
    help='Prefix for all the output files.')
  grpout.add_argument('--savecoeff',action='store_true',
    help='Save excitonic coefficient in a numpy file for later use.')
  grpout.add_argument('--savetprod','--tprod',action='store_true',
    help='Save matrix of triple products R_ij*(mu_i x mu_j).')

  # Site Energies
  grpsit = parser.add_argument_group(' Options for site energies')
  grpsit.add_argument('--insite','--site',metavar='SiteEnergyFile',
    help='Modify site energies according to SiteEnergyFile',default=None)

  # Centers options
  grpcnt = parser.add_argument_group(' Center options')
  read_grpcnt = grpcnt.add_mutually_exclusive_group()
  read_grpcnt.add_argument('--cent',metavar='Center',nargs="*",
    help='How to compute the chromophore center. \
    Default is the geometric center of the chromophore.',default="geom")
  read_grpcnt.add_argument('--incent',metavar='CenterFile',
    help='Modify chromophore centers according to CenterFile',default=None)


  # Couplings
  grpcou = parser.add_argument_group(' Options for couplings') 
  grpcou.add_argument('--coup', 
   choices=['trden','PDA','coulomb'],    
   help='''Choose the coupling calculation type. 
   trden includes all terms. 
   Coulomb includes only the coulomb term. 
   PDA computes coupling with point-dipole approxmation.
   Default is trden.''')
  grpcou.add_argument('--cleancoup',metavar='Threshold',type=float,
    help='Set to zero all the couplings below Threshold')
  grpcou.add_argument('--seltran',action="store_true",
    help="""Activate selection of transitions on the basis of ChromListFile. 
    ChromListFile should contain, for each line, the transitions to consider. """)
  grpcou.add_argument('--refrind',  metavar='n',type=float,
    help='''Specify the refraction index n for PDA coup calculations. 
    Default is 1.0 (vacuum).''')
  grpcou.add_argument('--scalecoup',metavar='Factor',type=float,
    help='Scale all couplings by Factor')
  grpcou.add_argument('--incoup','--modcoup',metavar='CoupFile',
    help='Modify electronic couplings according to the CoupFile',default=None)

  # Transition Dipoles
  grpdip = parser.add_argument_group(' Options for transition dipoles ')
  grpdip.add_argument('--scaledipo',type=float,
    help='Scale all transition dipole moments')
  grpdip.add_argument('--reorient',default=None,nargs=2,
    type=int,metavar=('ID1','ID2'),
    help='''Id of two atoms that define the axis 
    with respect to which the TrDip will be reoriented''')
  grpdip.add_argument('--forcedipo',action="store_true",
    help='Force TrDip to be aligned along the axis specified by --reorient keyword')
  grpdip.add_argument('--anadipo',metavar='FILE',
    help='''Analysis of the transition dipole orientation with respect to the
    reference frame indicated in FILE.''')
  grpdip.add_argument('--scaletran',metavar='ScaleFile',
    help='Request scaling of all dipoles and couplings as indicated in ScaleFile.')
  grpdip.add_argument('--indipo','--dipolen',metavar='DipoFile',default=None,
    help='Modify electric dipole moments according to DipoFile')

  # CD related
  rsttitle = ' CD related options: \n  Options to compute the rotational strength'
  grprst = parser.add_argument_group(rsttitle)
  rcalc_group = grprst.add_mutually_exclusive_group()
  rcalc_group.add_argument('--mu',action='store_true',
    help='Use approximate treatment (electric dipoles only). This is the default.')
  rcalc_group.add_argument('--mag',action='store_true',
    help='''Use full treatment (electric and magnetic dipoles), 
    splitting the electric and magnetic contributions''')
  
  # Advanced options
  grpadv = parser.add_argument_group(' Advanced and expert options ')
  grpadv.add_argument('--ldaxis',
    help='Choose the reference axis for LD calculation')

  # Verobisity of the output and printing options
  grpprt = parser.add_argument_group(' Printout and verbosity options')
  grpprt.add_argument('-v',action='count',
    help='Increase the verbosity of the output')
  grpprt.add_argument('-q',action="store_true",
    help='Quiet output: minimum printout')
  grpprt.add_argument('--prtcoup-threshold',metavar='Threshold',
    type=float,default=100.0,
    help='Threshold for coupling printout with -vv option enabled (default: 100cm^-1)')

  # Parse args
  args = parser.parse_args()

  # Set OPT dictionary
  if args.v > 0   : c.OPT['verbosity'] = args.v
  elif args.q     : c.OPT['verbosity'] = -1

  c.ExtFiles['crlist'] = args.chromlist

  # Check I/O options
  if args.readexternal:
    c.OPT['read'] = 'external'
    # With --readexternal, --seltran option is activated by default
    c.OPT['seltran'] = True  
    if args.insite is not None: c.ExtFiles['insite'] = args.insite
    if args.incoup is not None: c.ExtFiles['incoup'] = args.incoup
    if args.indipo is not None: c.ExtFiles['indipo'] = args.indipo
    if args.incent is not None: c.ExtFiles['incent'] = args.incent

  else:
    if not args.log:
      c.error("You have to indicate the Gaussian .log file!","Input parse")
    else:
      c.OPT['logfile']  = args.log
      c.OPT['read'] = 'g16'

  if args.mu == True       : c.OPT['RCalc'] = "mu"
  elif args.mag == True    : c.OPT['RCalc'] = "mag"
  else                     : c.OPT['RCalc'] = "mu"  # Default 

  if args.cent : c.OPT['Cent'] =  args.cent

  if args.out      : 
    c.OPT['OutPrefix'] = args.out
    c.setoutfiles()
  if args.savecoeff:
    c.OPT['savecoeff'] = True
  if args.savetprod:
    c.OPT['savetprod'] = True

  if args.seltran  : c.OPT['seltran']  = args.seltran
  if args.scalecoup is not None: c.OPT['ScaleCoup']  = args.scalecoup
  if args.scaledipo : c.OPT['ScaleDipo']  = args.scaledipo
  if args.reorient  : c.OPT['reorient']  = args.reorient

  if args.forcedipo  : 
    if c.OPT['reorient'] != None:
      c.OPT['forcedipo']  = args.forcedipo
    else:
      c.error("To force a dipole you have also to specify --reorient keyword",
      "argparse")

  if args.anadipo  : 
    c.OPT['anadipo'] = True
    c.ExtFiles['refaxis'] = args.anadipo
  if args.scaletran : 
    c.OPT['ScaleTran']  = True 
    c.ExtFiles['ScaleTran']  = args.scaletran

  if args.insite is not None:
    c.OPT['ModSite']   = True
    c.ExtFiles['insite'] = args.insite

  if args.indipo is not None:
    c.OPT['ModDipoLen']  = True
    c.ExtFiles['indipo'] = args.indipo

  if args.incoup is not None:
    c.OPT['ModCoup']  = True
    c.ExtFiles['incoup'] = args.incoup

  if args.incent is not None:
    c.OPT['ModCent']  = True
    c.ExtFiles['incent'] = args.incent

  if args.coup      : c.OPT['coup']   = args.coup
  if args.refrind   : c.OPT['refrind']   = args.refrind
  if args.cleancoup : c.OPT['CleanCoup']   = args.cleancoup
  if args.ldaxis    : c.OPT['LDAxis'] = args.ldaxis
  setldaxis()

  c.OPT['CouPrtThr'] = args.prtcoup_threshold

  if c.v(1): c.PrintOpt()



# *****************************************************************************
#
# Compute the chromophore centers
#
def calchromcent():

  Cent = []
  Init = 0

  if "geom" in c.OPT['Cent'] : ICalc = "geom"
  elif "mass" in c.OPT['Cent'] : ICalc = "mass"
  else: 
    ICalc = "select"
    ISel  = c.stringconverter(c.OPT['Cent'][0])
    ISel  = (np.array(ISel)-1).tolist()
    if c.v():
      print(" The following atoms will be selected to compute chromophore center: ")
      print(ISel)

  for i in range(c.NChrom):
    End  = Init + c.NAtom[i]
    IAnu = c.anum[Init:End]
    IXYZ = c.xyz[Init:End]
    Init = End

    if    ICalc == "geom" : 
      ICent = np.average(IXYZ,axis=0)
    elif  ICalc == "mass" : 
      IMas = anu2awg(IAnu)
      ICent = np.average(IXYZ,axis=0,weights=IMas)
    elif  ICalc == "select" :
      SAnu = list(IAnu[i] for i in ISel)
      SXYZ = list(IXYZ[i] for i in ISel)
      SMas = anu2awg(SAnu)
      ICent = np.average(SXYZ,axis=0,weights=SMas)

    # Replicate the center for every transiton...
    for j in range(c.NTran[i]):
      #Cent.extend(calcent(IAnu,IXYZ))
      Cent.extend(ICent)

  Cent = np.reshape(Cent,(sum(c.NTran),3)).tolist()
  return Cent

# *****************************************************************************
#
# Change sign of couplings involving trans Tran an chrom Chrom
#
def changesign(Chrom,Tran,NChrom,NTran,Coup):
  z = 0
  for i in range(NChrom):
    for j in range(i+1,NChrom):
      for k in range(NTran[i]):
        for l in range(NTran[j]):
          if (Chrom == i and Tran == k ) or ( Chrom == j and Tran == l):
            if c.v(0) : 
              print("   ... ChgSgn for Chrom: %3d   tr: %3d    z=%3d  coup=%10.4f  (V: %2d %2d  -- %2d %2d)"\
                % (Chrom+1,Tran+1,z,Coup[z],i+1,k+1,j+1,l+1)) 
            Coup[z] = -Coup[z]
          z += 1
  return Coup

# *****************************************************************************
#
# Reorient transition dipoles and change the coupling sign
#
def reorientdipo(DipoLen,DipoVel,Mag,Coup):
  Init = 0
  At1,At2 = c.OPT['reorient'] # direction: from 1 to 2
  K=0
  for I in range(c.NChrom):
    # Extract xyz of Chrom I
    End  = Init + c.NAtom[I]
    IXYZ = np.array(c.xyz[Init:End])
    Init = End
    # Define axis 
    axis  = IXYZ[At2-1]-IXYZ[At1-1]
    # Normalize
    axis /= np.linalg.norm(axis)
    for J in range(c.NTran[I]):
      sign = np.sign(np.dot(axis,DipoLen[K]))

      if sign < 0 : 
        if c.v():
          print("Chrom: %3d  Trans: %3d  Transition electric dipole moment will be reoriented"\
             % (I+1,J+1))
        DipoLen[K] *= -1 
        # Change signs to DipoVel and Mag as well
        DipoVel[K] *= -1 
        Mag[K]     *= -1 
        Coup = changesign(I,J,c.NChrom,c.NTran,Coup)

      if c.OPT['forcedipo'] == True:
        if c.NTran[I] > 1 : c.error("NTran > 1: you cannot force all dipoles to be parallel to a specific axis!","reorientdipo")
        DipoLen[K] = axis*np.linalg.norm(DipoLen[I])

      K += 1

  return DipoLen,DipoVel,Mag,Coup

# *****************************************************************************
#
# Transition dipole analysis
#
def buildrefframe(A,B,C):
  # Creates a reference frame from three points.
  # x is the A->B axis, and y is in the ABC plane
  # z = x cross y is defined by the couterclockwise ABC rotation

  ux = (B-A)/np.linalg.norm(B-A)  # Versor defined A ----> B
  AC = (C-A)
  ACprojx = np.dot(AC,ux) # projection on ux
  # check quasi-linear dependence
  if abs(ACprojx/np.linalg.norm(AC)) > np.cos(np.radians(10.0)) and c.v():
    print(" WARNING: choosen atoms are almost aligned!")
    print(" Angle between vectors: %8.2f degrees; threshold: %8.2f degrees" \
    % (np.degrees(np.arccos(ACprojx/np.linalg.norm(AC))),10))  
  # subtract the projection on x and normalize
  AC -= ACprojx*ux
  uy  = AC/np.linalg.norm(AC)
  # create z axis (x cross y)
  uz = np.cross(ux,uy)
  return ux,uy,uz

def changerefframe(oldframe,newframe,V):
  # Gives components of V in the new reference frame
  # works only for orthonormal ref frames
  transform = np.dot(oldframe,newframe.T)
  U = np.dot(transform,V)
  return U

def diporient(IChrom,ITranStrt,Chrom,NTran,DipoLen,frame,Theta,Phi,ChangeTheta,ChangePhi):
  # Orient dipoles for one chromophore, return new dipole 
  # frame: reference frame for this chromophore
  i = IChrom
  k = ITranStrt
  for j in range(NTran[i]):
    if ChangeTheta[k] is None and ChangePhi[k] is None:
      k+=1; continue
    # Dipole magnitude
    MDipo    = np.linalg.norm(DipoLen[k])
    # Change theta or not
    # newtheta in standard definition
    if ChangeTheta[k] is None: newtheta = np.radians(90.0 - Theta[k])
    else: newtheta = np.radians(90.0 - ChangeTheta[k])
    # Change phi or not
    if ChangePhi[k] is None: newphi = np.radians(Phi[k])
    else: newphi = np.radians(ChangePhi[k])

    # Coordinates wrt "new" reference frame (ux,uy,uz)
    DipoU    = np.empty(3)
    DipoU[0] = MDipo*np.sin(newtheta)*np.cos(newphi) # ux
    DipoU[1] = MDipo*np.sin(newtheta)*np.sin(newphi) # uy
    DipoU[2] = MDipo*np.cos(newtheta)                # uz
    # Go back to global reference frame
    DipoLen[k] = changerefframe(frame,np.identity(3),DipoU)

    if c.v():
      print(" Chrom : %10s   Tran : %3d   theta = %10.4f  phi = %10.4f   Dipole = %8.4f %8.4f %8.4f" \
      % (Chrom,j+1,90-np.degrees(newtheta),np.degrees(newphi),DipoLen[k,0],DipoLen[k,1],DipoLen[k,2]))
    k += 1

  return DipoLen
 
def dipoanalysis(DipoLen):
    
  NCols = 4
  THETA = []
  PHI   = []

  if c.v():
    print("\n ==== Transition Dipole Analysis Function ==== \n")
    print(" Computing angle between transition dipole moments and a given reference frame")
    print(" The reference frame is defined from the positions of three (groups of) atoms ")
    print(" Note that the sign of the angle is related to the third (group of) atoms") 
    print() 

  c.checkfile(c.ExtFiles['refaxis'])  
  if c.v():
    print(" > Reference file: %s \n" % c.ExtFiles['refaxis'])

  # Read reference file
  InFile = open(c.ExtFiles['refaxis'],'r')
  InData = InFile.readlines()
  InFile.close()

  # Process reference file
  NLines = len(InData)
  if NLines != c.NChrom: 
    c.error("Reference file must contains the same number of chromphores of the input","dipoanalysis")
  RefChrom = []
  AA = [] ; BB = [] ; CC = []
  ChangeTheta = []
  ChangePhi   = []


  Init = 0
  for i in range(c.NChrom):
    # Line for chromophore i
    Data = InData[i].split()
    if len(Data) < NCols : c.error("Check reference input! Number of coulmns for at least one data is wrong!","dipoanalysis")
    RefChrom.append(Data[0])
    RefGrpA = c.stringconverter(Data[1])
    RefGrpA = (np.array(RefGrpA)-1).tolist()
    RefGrpB = c.stringconverter(Data[2])
    RefGrpB = (np.array(RefGrpB)-1).tolist()
    RefGrpC = c.stringconverter(Data[3])
    RefGrpC = (np.array(RefGrpC)-1).tolist()
    End  = Init + c.NAtom[i]
    IXYZ = c.xyz[Init:End]
    AA.append(np.average(list(IXYZ[j] for j in RefGrpA),axis=0))
    BB.append(np.average(list(IXYZ[j] for j in RefGrpB),axis=0))
    CC.append(np.average(list(IXYZ[j] for j in RefGrpC),axis=0))

    # Read whether to change the dipole orientation
    for j in range(c.NTran[i]):
      # theta[Chrom i, Tr j] in Data[col]
      # phi    "        "    in Data[col+1] 
      col = NCols - 2 + (j+1)*2
      if len(Data) < col: 
        #fill ChangeTheta and ChangePhi lists
        nfill = c.NTran[i]-j
        ChangeTheta.extend([None]*nfill)
        ChangePhi.extend(  [None]*nfill)
        break #no need to go further
      # go on (there may be something)
      try: ChangeTheta.append(float(Data[col]))
      except: ChangeTheta.append(None)
      # see phi
      try:  ChangePhi.append(float(Data[col+1]))
      except: ChangePhi.append(None)
      
    Init = End

  ChangeTheta = np.array(ChangeTheta)
  ChangePhi   = np.array(ChangePhi)

  Theta = []; Phi = []
  
  k = 0
  for i in range(c.NChrom):
    # Build new reference frame
    ux,uy,uz = buildrefframe(AA[i],BB[i],CC[i])
    for j in range(c.NTran[i]):
      # Evaluate the angle
      Dipo  = DipoLen[k]
      UDipo = Dipo/np.linalg.norm(Dipo) 
  
      theta = np.arccos(np.dot(UDipo,uz)) 
      theta = 90.0 - np.degrees(theta) # angle wrt the plane
  
      phi   = np.arctan2(np.dot(UDipo,uy),np.dot(UDipo,ux))
      phi   = np.degrees(phi) # angle wrt x axis (AB)
  
      THETA.append(theta)
      PHI.append(phi)
      if c.v():
        print(" Chrom : %10s   Tran : %3d   theta = %10.4f  phi = %10.4f"\
        % (c.ChromList[i],j+1,theta,phi))
      if c.v(0): 
        print(" Dipole: " + 12*" " + " X= %10.6f   Y= %10.6f   Z= %10.6f " % tuple(Dipo))
      # save values
      Theta.append(theta)
      Phi.append(phi)
      k += 1

  Theta = np.array(Theta); Phi = np.array(Phi)
  
  if c.v(): print() 

  # Changing dipoles
  ChangeAnyTheta = not (ChangeTheta == [None]*len(ChangeTheta)).all()
  ChangeAnyPhi   = not (ChangePhi == [None]*len(ChangePhi)).all()

  if ChangeAnyTheta or ChangeAnyPhi: 
    if c.v(): 
      print(" ... Dipoles of the following transitions will be oriented: ")
    K = 0
    for i in range(c.NChrom): 
      refframe = np.array(buildrefframe(AA[i],BB[i],CC[i])).T
      DipoLen = diporient(i,K,c.ChromList[i],c.NTran,DipoLen,refframe,Theta,Phi,ChangeTheta,ChangePhi)
      K += c.NTran[i] # Go ahead to next chrom
  elif c.v():
    print(" ... Dipoles will not be changed")

  # End
  if c.v():
    print("\n ====  End of Transition Dipole Analysis  ==== \n")
  return (np.array(THETA),np.array(PHI))


# *****************************************************************************
#
# Transition dipoles and coupling scaling (ScaleTran)
#

# Helper function for couplings
#Scales by Fact all couplings that involve transition Tran of chromophore Chrom

def scalecouptr(Fact,Chrom,Tran,NChrom,Coup):
  z = 0
  for i in range(NChrom):
    for j in range(i+1,NChrom):
      for k in range(c.NTran[i]):
        for l in range(c.NTran[j]):
          if (Chrom == i and Tran == k ) or ( Chrom == j and Tran == l):
            if c.v(0) : print(("   ... ScaleCoup for Chrom: %3d   tr: %3d    z=%3d  coup=%10.4f  (V: %2d %2d  -- %2d %2d)" % (Chrom+1,Tran+1,z,Coup[z],i+1,k+1,j+1,l+1)))
            Coup[z] *= Fact
          z += 1
  return Coup

# scaletran function
# Scales electric and magnetic dipoles of transition as defined in c.ExtFiles['ScaleTran']
# Scales all couplings accordingly

def scaletran(DipoLen,Coup,Mag):
  reffile = c.ExtFiles['ScaleTran']
  c.checkfile(reffile)  
  if c.v():
    print() 
    print(" ... scaling transition densities using scaling factors from file: %s" % reffile) 

  # Read reference file
  InFile = open(reffile,'r')
  InData = InFile.readlines()
  InFile.close()

  # Process reference file
  NLines = len(InData)
  if NLines != c.NChrom: 
    c.error("Reference file must contains the same number of chromphores of the input","scaletran")

  k = 0
  for i in range(c.NChrom):
    # Line for chromophore i
    Data = InData[i].split()

    for j in range(c.NTran[i]):
      # read whether to scale something in tran j (column j+1)
      col = j+1
      if len(Data) < col: break # no need to go further 
      try:    scalefact = float(Data[col])
      except: scalefact = None
      # Scale dipo and couplings of tran j 
      if scalefact is not None and scalefact != 1.0:
        if c.v():
          print("     Chrom: %10s   Tran: %3d  scale transition density by %10.6f" % (c.ChromList[i],j+1,scalefact))
        # scale Dipole 
        DipoLen[k] *= scalefact
        Mag[k]     *= scalefact
        # scale Coupling
        Coup = scalecouptr(scalefact,i,j,c.NChrom,Coup)
      k += 1
  pass

# *****************************************************************************
#
# Modify site energies 
#
def modsite(Site):
  reffile = c.ExtFiles['insite']
  c.checkfile(reffile)  

  # Read reference file
  with open(reffile,'r') as InFile:
    InData = InFile.readlines()
  # Process reference file
  NLines = len(InData)
  if NLines != c.NChrom : c.error("Reference file must contain the same number of chromphores of the input")
  k = 0
  for i in range(c.NChrom):
    shift = False
    if "+" in InData[i] or "-" in InData[i]: shift = True  # check if we want to shift or not
    Data = InData[i].split()
    for j in range(c.NTran[i]):
      # read whether to change site energy of tran j  (column j+1)
      col = j+1
      if len(Data) < col: break # no need to go further 
      try:     Ene = float(Data[col])
      except:  Ene = None
      if Ene is not None:
        if shift == False:
          if c.v():
            print("     Chrom: %10s   Tran: %3d  site energy set to %8.5f eV"\
                  % (c.ChromList[i],j+1,Ene))
          Site[k] = Ene
        elif shift == True:
          if c.v():
            print("     Chrom: %10s   Tran: %3d  site energy sfhit of %8.5f eV"\
                  % (c.ChromList[i],j+1,Ene))
          Site[k] += Ene
      k += 1
  return Site

# *****************************************************************************
#
# Modify electric transition dipole moments
#
def moddipo(Dipo):
  reffile = c.ExtFiles['indipo']
  c.checkfile(reffile)
  if c.v():
    print(" Dipoles will be read from %s" % (reffile))

  # Read reference file (expected the same format as the exat output)
  InData = np.loadtxt(reffile)
  InData /= c.PhyCon['ToDeb']
  NLines = len(InData)
  if NLines != sum(c.NTran) : c.error("Reference file must contain the same number of transitions of the input")

  # Process reference file
  if c.v():
    n = 0
    for i in range(c.NChrom):
      for k in range(c.NTran[i]):
        fmt = "Chrom: %3d  Tran: %3d  Orig dipo = %8.4f %8.4f %8.4f    New dipo = %8.4f %8.4f %8.4f"
        print(fmt % (i+1,k+1,Dipo[n][0],Dipo[n][1],Dipo[n][2],InData[n][0],InData[n][1],InData[n][2]))
        n += 1

  return InData

# *****************************************************************************
#
# Modify electronic couplings
#
def modcoup(NCoup,Coup):
  reffile = c.ExtFiles['incoup']
  c.checkfile(reffile)
  if c.v():
    print(" ... couplings will be subsituted with those in %s file" % reffile)
  Coup = np.loadtxt(reffile)
  return Coup


# *****************************************************************************
#
# Assign the atomic weigths
#

def anu2awg(In):
  # Return atomic masses 
  Out = []
  for I in In:
    Out.append(c.AtMass[I])
  return Out

# Save the complete geometry in xyz file
def savegeom():
  if c.v(1): print(" ... save geometry to file %s " %  c.OutFiles['xyz']) 
  NAtoms = len(c.anum)
  OutF = open(c.OutFiles['xyz'],'w')
  OutF.write("%d\n" %  NAtoms)
  OutF.write("Generated by %s\n" %  c.PROGVERS)
  for i in range(NAtoms):
    OutF.write("%3d %10.5f %10.5f %10.5f \n" % tuple([c.anum[i]]+c.xyz[i])) 
  return


# *****************************************************************************
#
# Save a .vmd file to visualize dipoles in VMD
# Usage: vmd -e file.vmd
#
def savevisudipo(C,Dipo,ExcDipo,MagDipo=None):

  if c.v(1): 
    print(" ... save dipole visualization script to %s " %  c.OutFiles['visudipo']) 
    print("     Edit the file following the instructions therein, then type ")
    print("     vmd -e %s " % c.OutFiles['visudipo'])

  OutFN = c.OutFiles['visudipo']
  with  open(OutFN,'w') as OutF:
    OutF.write(c.VMDHeader())
  
    Scale    = 2.0*c.PhyCon['ToDeb']
    ScaleMag = 2.0
    
    OutF.write("mol new %s type xyz\n" % c.OutFiles['xyz'])
  
    # Write down site dipoles
    OutF.write("\n## SITE DIPOLES ...  \n\n")
    
    k = 0
    for i in range(c.NChrom):
      OutF.write(" \n")
      for j in range(c.NTran[i]):
        OutF.write(" \n")
        OutF.write("## CHROM %6s TRAN %3d \n" % (c.ChromList[i],j+1))
        # Default: visualize only electric dipole of first transition
        OutF.write("## Electric Transition Dipole \n")
        if j != 0: OutF.write('#') 
        OutF.write("graphics 0 color red; ")
        OutF.write("vmd_draw_vector 0 { %10.4f %10.4f %10.4f } { %10.4f %10.4f %10.4f }\n"\
        % (tuple(C[k])+tuple(Dipo[k]*Scale)) )
        if MagDipo is not None:
          OutF.write("## Magnetic Transition Dipole \n")
          OutF.write("#graphics 0 color blue; ") # Default: comment this
          OutF.write("vmd_draw_vector 0 { %10.4f %10.4f %10.4f } { %10.4f %10.4f %10.4f }\n"\
          % (tuple(C[k])+tuple(MagDipo[k]*ScaleMag)) )
        
        k += 1
    #end for

    # Write down excitonic (electric) dipoles 
    totcent = np.average(C,axis=0)
    OutF.write("\n## EXCITONIC DIPOLES ...  \n\n")
  
    # All excitonic transitions

    for k in range(sum(c.NTran)): 
      OutF.write("## EXC STATE %3d \n" % (k+1))
      if k != 0: OutF.write('#') 
      OutF.write("graphics 0 color green; ")
      OutF.write("vmd_draw_vector 0 { %10.4f %10.4f %10.4f } { %10.4f %10.4f %10.4f }\n"\
      % (tuple(totcent)+tuple(ExcDipo[k]*Scale)) )


    # End
    OutF.write("display resetview\n")
  # file closed

  # Save dipoles
  OutFN = c.OutFiles['dipo']
  formt = "%10.4f %10.4f %10.4f"
  np.savetxt(OutFN,Dipo*c.PhyCon['ToDeb'],fmt=formt,header="Dipoles in debye")

  # Save excitonic dipoles
  # TODO

  if MagDipo is not None:
    # Save magnetic dipoles
    OutFN = c.OutFiles['magdipo']
    np.savetxt(OutFN,MagDipo,fmt=formt,header="Magnetic Dipoles")


# *****************************************************************************
#
# Save coupling values and distances between chromophore pairs
# 
#
def savecoupdist(Coup,Cent):
  OutCoup = open(c.ExtFiles['out_coup'],'w')
  n = 0
  for i in range(c.NChrom):
    for j in range(i+1,c.NChrom):
      for k in range(c.NTran[i]):
        for l in range(c.NTran[j]):
          ii = sum(c.NTran[0:i])+k
          jj = sum(c.NTran[0:j])+l
          dist = np.linalg.norm(np.array(Cent[jj])-np.array(Cent[ii]))
          OutCoup.write("%5s (%2d) %5s (%2d)  %10.4f %10.4f\n" % (c.ChromList[i],k+1,c.ChromList[j],l+1,dist,Coup[n]))
          n += 1

  OutCoup.close()



# *****************************************************************************
#
# Fill the NTran array with NTranIn transition for chromophore
#

def fillNTran(NTranIn,NChrom):

  # Check if number of chromophores is greater than 1
  if NChrom < 2:
    c.error(" We need at least two chromophores ",'fillNTran')

  # Fill the NTran list
  N = len(NTranIn)
  if N == 1:
    NTran = NTranIn*NChrom
  elif N == NChrom:
    NTran = NTranIn
  else:
    c.error(" Inconsitency between NChrom and NTran ",'fillNTran')

  return NTran


# *****************************************************************************
#
# Compute the oscillator strenghts
#

# TO BE WRITTEN!


#
# =============================================================================
# CONVERSION UTILITIES
# =============================================================================
#

# List of physic constants
me       = 9.10938291E-31    # mass of the electron (Kg)
ee       = 1.60217656E-19    # charge of the electron (C)
cc        = 299792458         # speed of light (m/s)
h        = 6.62606957E-34    # plank constant (J s)
heV      = 4.13566752E-15    # plank constant (eV s)
hbar     = h/(2*np.pi)          # h/2pi (J s)
heVbar   = heV/(2*np.pi)        # h/2pi (eV s)
r0       = 0.52917721092     # Bohr Radius (Ang)

# List of conversion factors
Deb2Cm   = 3.33564E-30       # Debye (statC m)   ---> C m
au2Deb   = 2.54154           # a.u.              ---> Debye


# ***********************************************

#
# Save a results.out file
#

def resout(energy,dip,LD,R):
  EeV = energy/c.PhyCon['eV2wn']
  n = len(energy)
  outfile = c.OutFiles['results']
  with open(outfile,'w') as out:
    out.write('# Created with %s\n' % (c.PROGVERS))
    out.write('# Initial command: \n')
    out.write('# '+' '.join(sys.argv)+'\n')
    out.write('# State Energy     mu^2      LD       RotStr \n')
    for i in range(n):
      out.write('%3d %10.4f %10.4f %10.4f %10.4f \n'\
       % (i+1,EeV[i],dip[i],LD[i],R[i]))
    if c.v():
      print("   ... file %s containing the excitonic results has been saved!\n"\
       % outfile)
  pass

# ***********************************************


def coupforster(Cent,Dipo,IKappa=False):
  NChrom = c.NChrom
  Coup = [] ; Kappa = []
  n = 0
  for i in range(NChrom):
    for j in range(i+1,NChrom):
      for k in range(c.NTran[i]):
        for l in range(c.NTran[j]):
          ii = sum(c.NTran[0:i])+k
          jj = sum(c.NTran[0:j])+l
          Ci = np.array(Cent[ii])
          Cj = np.array(Cent[jj])
          Di = np.array(Dipo[ii])
          Dj = np.array(Dipo[jj])
          tcoup,tkappa = forster(Di,Dj,Ci,Cj)
          Coup.append(tcoup)
          Kappa.append(tkappa)
          n += 1
  if IKappa is True: return Coup,Kappa
  else: return Coup


#
# Compute FORSTER Couplings
# ----------------------------------------------
# D1 ... transition dipole 1 (a.u.)
# D2 ... transition dipole 2 (a.u.)
# C1 ... center of D1 (Ang)
# C2 ... center of D2 (Ang)
# ECXOpts['refrind']    -> refraction index
# MainOpts['verbosity'] -> verbosity output
# -----------------------------------------------
#
def forster(D1,D2,C1,C2):
  Ang2au      = c.PhyCon['ToAng']
  hartree2cm  = c.PhyCon['Town']
  s = 1.0/(c.OPT['refrind']**2)
  C1 = C1/Ang2au
  C2 = C2/Ang2au
  dist     = C2-C1
  normdist = np.linalg.norm(dist)
  versdist = dist/normdist
  normD1   = np.linalg.norm(D1)
  normD2   = np.linalg.norm(D2)
  versD1   = D1/normD1
  versD2   = D2/normD2
  fact1    = np.dot(versD1,versD2)
  kappa    = np.dot(versD1,versD2)-3*(np.dot(versD1,versdist)*np.dot(versD2,versdist))
  coupvac  = (kappa * normD1 * normD2) / (normdist**3)
  coup     = coupvac*s

  #if c.OPT['verbosity'] > 2:
  #  print "Kappa         : %8.3f"      % kappa
  #  if ( ECXOpts['refrind'] == 1.0 ):
  #    print "Coup     : %8.3f cm-1" % (coupvac*hartree2cm)
  #  else:
  #    print "Coup vacuo    : %8.3f cm-1" % (coupvac*hartree2cm)
  #    print "Coup (s=%4.2f) : %8.3f cm-1" % (s,coup*hartree2cm)

  return(coup*hartree2cm,kappa)



#
# Print out couplings and center distance
#
def prtcoup(Coup,Cent,Kappa=False):
  # Print out couplings larger than a threshold
  NTran  = c.NTran
  NChrom = c.NChrom
  Dist   = distance.cdist(Cent,Cent, 'euclidean')
  thresh = c.OPT['CouPrtThr']
  if thresh is False: return
  print("\n Couplings larger than %7.2f cm^-1: " % thresh)
  print(" Chrom Chrom  Tran  Tran   V (cm-1)  Dist (Ang)        Kappa")
  print(" -------------------------------------------------------------")
  z = 0
  OutCoupFile = open(c.OutFiles['coup'],'w')
  if c.OPT['coup'] != 'forster': Kappa = np.zeros((len(Coup)))
  for i in range(NChrom):
    for j in range(i+1,NChrom):
      for k in range(NTran[i]):
        for l in range(NTran[j]):
          if abs(Coup[z]) > thresh:
            # Print on screen
            print("%6s %6s  %3d %3d   %10.4f   %10.4f   %10.4f" \
            % (c.ChromList[i],c.ChromList[j],k+1,l+1,Coup[z],Dist[sum(NTran[0:i])+k,sum(NTran[0:j])+l],Kappa[z]))
          # Save all on output file
          OutCoupFile.write("%6s %6s  %3d %3d   %10.4f   %10.4f   %10.4f \n" \
          % (c.ChromList[i],c.ChromList[j],k+1,l+1,Coup[z],Dist[sum(NTran[0:i])+k,sum(NTran[0:j])+l],Kappa[z]))
          z += 1
  OutCoupFile.close()
  return


def prtsite(Site):
  OutFile = open(c.OutFiles['site'],'w')
  z = 0
  for i in range(c.NChrom):
    fmt   = "%3d   " + ("%10.4f"*c.NTran[i])+"\n"
    isite = tuple(Site[z:z+c.NTran[i]])
    prt   = tuple([i+1])+isite
    OutFile.write(fmt % prt)
    z += c.NTran[i]
  OutFile.close()
  return


def setldaxis():
  if c.OPT['LDAxis'] not in ('x','y','z'):
    axis = c.OPT['LDAxis'].split(',')
    if len(axis) != 3:
      c.error(' Wrong definition of LD Axis: %s' % c.OPT['LDAxis'])
    c.OPT['LDAxis'] = axis 
    


