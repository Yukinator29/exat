#!/usr/bin/env python

#
# Copyright 2013, Sandro Jurinovich
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.  
#

# Standard Python Modules
import sys,os,glob,fileinput,re,copy
import subprocess as sp
import numpy  as np
import argparse as arg

# Import COMMON Modules
import common as c
import util   as u

# *****************************************************************************
# 
# MAIN function to read data
#

def Read():
  
  if c.OPT['read'] == 'g16':
    if c.v():
      print  " ... excpected gaussian version: g16"
      print  " ... reading in %s .log file" % c.OPT['logfile']
    Site,Dipo,DipoVel,Mag,Coup = readgaulog36(c.OPT['logfile'])
  
  elif c.OPT['read'] == 'external':
  
    if c.v(): 
      print " ... all parameters are read from external files" 
      print " ... reading %s file to set the chromophore list" % c.ExtFiles['crlist'] 
    ReadChromList()
    Cent,Site,Dipo,DipoVel,Mag,Coup = ReadExternal()
  

  if c.OPT['read'] != 'external':
    if c.v():
      print " ... compute the center of each transtion" 
    Cent = u.calchromcent()
  
  # Compute dipole couplings
  if c.OPT['coup'] == 'PDA' :
    if c.v():
      print " ... using dipole-dipole couplings based on electric transition dipole moments" 
      print "     REFRACTION INDEX = %5.2f" % c.OPT['refrind'] 
    Coup,Kappa = u.coupforster(Cent,Dipo,IKappa=True)
  else:
    Kappa=False
  
  #print(" ... save couplings and distances between chromophres")
  #u.savecoupdist(Coup,Cent)
  
  # Applying selection of transitions 
  if c.OPT['read'] != 'external' and c.OPT['seltran'] == True: 
    Site,Dipo,DipoVel,Mag,Cent,Coup,Kappa = seltran(Site,Dipo,DipoVel,Mag,Cent,Coup,Kappa)
  
  Cent     =  np.array(Cent)
  DipoLen  =  np.array(Dipo)                    # Electric tranistion dipole moments (a.u.)
  DipoVel  =  np.array(DipoVel)                 # Nabla (a.u.)
  Mag      =  np.array(Mag)
  Site     = np.array(Site)
  Coup     = np.array(Coup)
  Kappa    = np.array(Kappa)

  if c.OPT['read'] != 'external' and c.OPT['ModCent'] == True:
    if c.v: print " ... centers are replaced with external data" 
    c.checkfile(c.ExtFiles['incent'])
    NewCent = np.loadtxt(c.ExtFiles['incent'],dtype="float")
    if c.v(1):
      for i in range(len(NewCent)):
        print "     (%3d) ReadCent - LogCent = %12.4f Ang" % (i+1,np.linalg.norm((NewCent[i]-Cent[i])))
    Cent = NewCent  
  return Cent,DipoLen,DipoVel,Mag,Site,Coup,Kappa


# *****************************************************************************
# 
# Read data from external files
#
def ReadExternal():

  # Read External files:
  if c.v():
    print "     > CHROMLIST         : %s" % c.ExtFiles['crlist']
    print "     > site energies     : %s" % c.ExtFiles['insite']
    print "     > coupling          : %s" % c.ExtFiles['incoup']
    print "     > elec. tr. dipoles : %s" % c.ExtFiles['indipo']
    print "     > chrom centers     : %s" % c.ExtFiles['incent']

  # Read chromlist to identify number of chromophores 
  # and number of transitions
  ChromList,NChrom,NTran  =  ReadChromList()

  c.NTran = []
  for i in NTran:
    c.NTran.append(len(i))

  if c.v():  
    print " ... number of chromophores         : %3d" % NChrom 
    print " ... N. transitions per chromophore : %s"  % c.NTran 

  # Read Site Energies
  Site = np.zeros(sum(c.NTran))
  Site = u.modsite(Site)

  # Check if the file exists and read the files
  c.checkfile(c.ExtFiles['incoup'])
  c.checkfile(c.ExtFiles['incent'])
  c.checkfile(c.ExtFiles['indipo'])
  Coup = np.loadtxt(c.ExtFiles['incoup'],dtype="float",ndmin=1)
  Cent = np.loadtxt(c.ExtFiles['incent'],dtype="float")
  Dipo = np.loadtxt(c.ExtFiles['indipo'],dtype="float")

  DipoVel = np.copy(Dipo)
  Mag     = np.copy(Dipo)
  Cent    = np.array(Cent)
  Site    = np.array(Site)
  Dipo    = np.array(Dipo)/c.PhyCon['ToDeb']
  DipoVel = np.array(DipoVel)
  Mag     = np.array(Mag)
  Coup    = np.array(Coup)


  return Cent,Site,Dipo,DipoVel,Mag,Coup



# *****************************************************************************
# 
# Read ChromList file
#
def ReadChromList():

  InFile = c.ExtFiles['crlist']
  c.checkfile(InFile)
  File = open(InFile,'r').readlines()
  NChrom = len(File)

  ChromList = [] ; NTran = []
  for i in range(NChrom):
    Tmp    = File[i].split()
    IChrom = Tmp[0]
    ChromList.append(IChrom)
    if c.OPT['seltran'] == True:
      ITran  = map(int,Tmp[1:])
      if not ITran :
        ErrMsg = "Chromophore %s has no transition selected" % Tmp
        c.error(ErrMsg,"ReadChromList")
      else:
        NTran.append(ITran)
  #Copy variables to common
  c.ChromList = ChromList
  c.NChrom    = NChrom
  return ChromList,NChrom,NTran



# *****************************************************************************
#
# Read into gdvh36 .log files to extract site energies, dipoles ...
# Also extract couplings
#

def readgaulog36(logfile):

  # Check if the file exists
  c.checkfile(logfile)

  # Inizialize lists
  xyz  = [] ; anum = [] ; site = []; dipo = [] ; dipovel = [] ; mag = [] ; NAtoms = []

  # Initialize NTran to Gaussian default
  # to avoid errors later
  NTran = 3

  # What coupling to look for
  if   c.OPT['coup'] == 'trden' or c.OPT['coup'] == 'total':  
    if c.OPT['read'] == 'gdvh36': 
      FindCoup = "> TOTAL" 
    else:
      FindCoup = "Total coupling"
    DoCoup   = True 
  elif c.OPT['coup'] == 'coulomb':  
    if c.OPT['read'] == 'gdvh36': 
      FindCoup = "> Coulomb"
    else: 
      FindCoup = "Coulomb"
    DoCoup   = True 
  elif c.OPT['coup'] == 'PDA':
    print("   ... skipping couplings reading: coupling will be computed using point dipole approximation\n")
    DoCoup   = False 
  else:
    c.error("Confused in coupling choice!","readgaulog36")


  # If the file exists then open it and read line by line
  with open(logfile,'r') as f:
    # We split the file reading in 2 loops to improve performance
    # Read properties
    while True:  
      line = f.readline()
      if not line: break
      if 'Electronic Coupling for Excitation Energy Tranfer' in line: break

      # Looks for NTran:
      if "9/41=" in line:
        try:
          NTran = int(line.split("9/41=")[1].split(",")[0])
        except:
          NTran = int(line.split("9/41=")[1].split("/")[0])
          

      # Looks for NChrom:
      if "62=" in line:
        NChrom = int(line.split("62=")[1].split(",")[0])
        NAtoms = [0]*NChrom 
        FragAt = [None]*NChrom


      # Extract number of atoms in gdvh36
      if line[1:19] == 'Symbolic Z-matrix:':
        kk = 0
        while True:
          line = f.readline()
          if line[0:9] == ' Charge =':
            pass 
          elif '(fragment=' in line.lower():
            d = line.split()[0]
            IFrag = int(d.split('=')[1].split(')')[0].split(',')[0])
            NAtoms[IFrag-1] += 1
            # Assign atom to fragment
            try:    FragAt[IFrag-1].append(kk)
            except: FragAt[IFrag-1] = [kk]
            kk += 1
          elif "NAtoms=" in line:
            NAtmTot = int(line[11:16])
            break

      # Looks for the atomic coordinates 
      elif "Input orientation:"  in line:
        atom = True
        j = 0
        while atom == True :
          line = f.readline()
          # Skip 4 lines
          if j < 4: j += 1;  continue
          tempcen = line.split()
          if len(tempcen) != 6:
            atom = False
          else:
            xyz.append(map(float,tempcen[3:6]))
            anum.append(int(tempcen[1]))
            j = j + 1

      # Extracts the dipole moments (length)
      elif "electric dipole" in line:
        atom = True
        j = 0
        while atom == True :
          pos  = f.tell()
          line = f.readline()
          # Skip 1 line
          if j < 1: j += 1;  continue
          tempdip = line.split()
          if len(tempdip) != 6 : 
            f.seek(pos) #go back one line
            break
          else:
            dipo.append(map(float,tempdip[1:4]))
            j = j + 1

      # Extracts the dipole moments (velocity)
      elif "transition velocity dipole" in line:
        atom = True
        j = 0
        while atom == True :
          pos  = f.tell()
          line = f.readline()
          if j < 1: j += 1;  continue
          tempdip = line.split()
          if len(tempdip) != 6 : 
            f.seek(pos) #go back one line
            break
          else:
            dipovel.append(map(float,tempdip[1:4]))
            j = j + 1

      # Extracts the magnetic moment
      elif "transition magnetic dipole" in line:
        atom = True
        j = 0
        while atom == True :
          pos  = f.tell()
          line = f.readline()
          if j < 1: j += 1;  continue
          tempdip = line.split()
          if len(tempdip) != 4 : 
            f.seek(pos) #go back one line
            break
          else:
            mag.append(map(float,tempdip[1:4]))
            j = j + 1

      # Extracts the site energy values
      elif "nm " in line:
        site.append(float(line.split()[4]))
    #end while

    # read couplings
    if DoCoup: 
      if c.v(): print  " ... reading couplings in %s file using %s values" % (c.OPT['logfile'],c.OPT['coup'])
      Cdtyp = [('Ch1',int),('Tr1',int),('Ch2',int),('Tr2',int),('Coup',float)]
      Coup = []

      if c.OPT['read'] == 'gdvh36':
        # read couplings for gdvh36-molecolab 
        if c.OPT['coup'] == 'total': ExplCoup = []
        while True:  
          line = f.readline()
          if not line: break
          elif FindCoup in line: 
            Coup += [(int(line[25:29]),int(line[31:35]),int(line[43:47]),int(line[49:53]),float(line[55:69]))]
          elif c.OPT['coup'] == 'total' and '> Explicit MMPol'  in line:
            ExplCoup += [(int(line[25:29]),int(line[31:35]),int(line[43:47]),int(line[49:53]),float(line[55:69]))]
      else:
        # read couplings for gdvh36 (plain)
        while True:
          line = f.readline()
          if not line: break
          elif 'Frag=' in line and 'State=' in line:
            Ch1,St1,Ch2,St2 = (int(line[6:9]),int(line[16:19]),int(line[45:48]),int(line[55:58]))
            # Read other lines
            while True:
              line = f.readline()
              if not line: break
              if FindCoup in line: 
                CoupValue = float(line[31:43])*c.PhyCon['eV2wn']
                break
            # swap chrom 1 and 2 for compatibility with order of gdvh36-molecolab
            Coup += [(Ch2,St2,Ch1,St1,CoupValue)]

      # Transform in np.array and sort 
      Coup =  np.array(Coup,dtype=Cdtyp)
      Coup.sort(order=['Ch1','Ch2','Tr1','Tr2','Coup'],axis=0)
      Coup = Coup['Coup']# Get rid of indices
      if c.OPT['coup'] == 'total' and len(ExplCoup) == len(Coup): 
      # Compute and subtract explicit term
        ExplCoup =  np.array(ExplCoup,dtype=Cdtyp)
        ExplCoup.sort(order=['Ch1','Ch2','Tr1','Tr2','Coup'],axis=0)
        Coup -= ExplCoup['Coup']# Get rid of indices
         
    # End read couplings
    else:
      print("   ... skipping couplings reading: coupling will be computed using point dipole approximation\n")
      Coup = None
  # File is closed
  
  #Transform NTran in list for compatibility
  NTran = [NTran]*NChrom

  # Reorder anum and xyz frag-wise order
  anum1 = np.array(anum)
  xyz1  = np.array(xyz)
  anum = []; xyz = []
  for i in range(NChrom):
    anum = anum + anum1[FragAt[i]].tolist()
    xyz  = xyz  +  xyz1[FragAt[i]].tolist()

  # Save variables in common
  c.ChromList = range(1,NChrom+1)
  c.NAtom     = NAtoms
  c.NChrom    = NChrom 
  c.NTran     = NTran
  c.anum      = anum
  c.xyz       = xyz
  

  return site,dipo,dipovel,mag,Coup


# 
#  Select properties
#

def delist(List,NChrom,IndChrom,Sel):
  k = 0
  NewList = []
  for IChrom in range(NChrom):
    for ITran in range(c.NTran[IChrom]):
      if IChrom in IndChrom:
        Ind = IndChrom.index(IChrom)
        if ITran+1 in Sel[Ind]:
          NewList.append(List[k])
      k += 1
  return NewList

def seltran(Site,Dipo,DipoVel,Mag,Cent,Coup,Kappa=False):

  NChrom = c.NChrom
  # Read the external file containg the selected transtions
  IFile = c.ExtFiles['crlist']
  if c.v():
    print "   ... transition of interests will be selected on the basis of %s file" % IFile  
  SelChromList,SelNChrom,Sel = ReadChromList()

  if c.OPT['read'] == 'gdvh36': 
    ChromList = map(str,range(1,NChrom+1))
  else:
    ChromList = c.ChromList

  IndChrom = [ ChromList.index(x) for x in SelChromList ]


  Site      = delist(Site,NChrom,IndChrom,Sel)
  Dipo      = delist(Dipo,NChrom,IndChrom,Sel)
  DipoVel   = delist(DipoVel,NChrom,IndChrom,Sel)
  Mag       = delist(Mag,NChrom,IndChrom,Sel)
  Cent      = delist(Cent,NChrom,IndChrom,Sel)

  #
  # Select couplings
  #

  k = 0
  NewCoup = [] ; NewKappa = []
  for I in range(NChrom):
    for J in range(I+1,NChrom):
      for it in range(c.NTran[I]):
        for jt in range(c.NTran[J]):
          if I in IndChrom and J in IndChrom:
            IndI = IndChrom.index(I)
            IndJ = IndChrom.index(J)
            if it+1 in Sel[IndI] and jt+1 in Sel[IndJ]:
              NewCoup.append(Coup[k])
              if Kappa is not False: NewKappa.append(Kappa[k])
          k += 1

  #
  # Select geometry (for later usage)
  #
  if c.OPT['read'] == 'gdvh36' and SelNChrom != NChrom:
    if c.v():
      print " "*7+"cutting geometry, leaving only selected chromophores"
    xyz = []; anum = []; NAtom = []
    End = 0
    for i in range(NChrom):
      Start = End
      End   = Start + c.NAtom[i] 
      if i in IndChrom: 
        # Append to new list
        xyz  = xyz  +  c.xyz[Start:End]
        anum = anum + c.anum[Start:End]
        NAtom.append(c.NAtom[i])
    c.NAtom = NAtom
    c.anum  = anum
    c.xyz   = xyz 
  #
  # Reset NTran and NChrom properly
  #
  NTran = []
  for i in Sel:
    NTran.append(len(i))
  NChrom = len(NTran)

  c.NChrom = NChrom
  c.NTran = NTran # Save NTran to Common 

  return Site,Dipo,DipoVel,Mag,Cent,NewCoup,NewKappa

# *****************************************************************************
 
def printlocal(Site,Dipo,Dip2,Cent):

  for n in range(c.OPT['NChrom']):
    print("\nChromophore     : %10s" % n )
#    print("Center of trans : %10.4f %10.4f %10.4f" % (Cent[n][0],Cent[n][1],Cent[n][2]))
    print'----------------------------------------------------------'
    print' #   E (eV)      mx      my     mz        Dip2   (a.u.) '
    print'----------------------------------------------------------'
#
    for i in range(c.OPT['NTran'][n]):
      k = sum(c.OPT['NTran'][0:n])+i
      print("%2d %8.4f   %6.3f  %6.3f  %6.3f    %6.3f "%(i+1,Site[k],Dipo[k][0],Dipo[k][1],Dipo[k][2],Dip2[k]))

# *****************************************************************************

