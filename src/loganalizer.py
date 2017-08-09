#!/usr/bin/env python

#
# loganalizer.py is a program to analize your gaussian output
#

# Standard Python Modules
import sys,os,glob
import numpy    as np
import argparse as arg
import common   as c
import util     as u
import readdata as r




#  # *****************************************************************************
#  #
#  # Compute the angle between two vectors.
#  #
#  
#  def calcangle(RefVec,Vec):
#  
#    DotProd  = np.dot(RefVec,Vec)
#    NormProd = np.linalg.norm(RefVec)*np.linalg.norm(Vec)
#    Angle    = np.degrees(np.arccos(DotProd/NormProd))
#    Sign     = np.sign(DotProd)
#  
#    return Angle,Sign
#  
#  
#  
#  
#  # *****************************************************************************
#  #
#  # Read into Gaussian.log files to extract site energies, dipoles ...
#  # This function is compatible with both the version of Gaussian.
#  # In the case of gdvh23 should be called for each file.
#  #
#  
#  def readgaulog(logfile):
#  
#    # If the file exists then open it, read all data and close it
#    infile = open(logfile,'r')
#    data = infile.read().split("\n")
#    infile.close()
#  
#    # Inizialize lists
#    xyz  = [] ; anum = [] ; site = []; dipo = [] ; dipovel = [] ; mag = [] ; NAtoms = []
#    rotstr = []
#  
#    # Read the loaded file line by line
#    i = 80
#    while i < len(data):
#  
#      # Looks for the atomic coordinates 
#      if data[i].find("Input orientation:") != -1:
#      #if data[i].find("Standard orientation:") != -1:
#        atom = True
#        j = 0
#        while atom == True :
#          tempcen = data[i+5+j].split()
#          if len(tempcen) != 6:
#            atom = False
#          else:
#            xyz.append(map(float,tempcen[3:6]))
#            anum.append(int(tempcen[1]))
#            j = j + 1
#  
#      # Extracts the dipole moments (length)
#      elif data[i].find("electric dipole") != -1:
#        atom = True
#        j = 0
#        while atom == True :
#          tempdip = data[i+2+j].split()
#          if len(tempdip) != 6 :
#            atom = False
#          else:
#            dipo.append(map(float,tempdip[1:4]))
#            j = j + 1
#  
#      # Extracts the dipole moments (velocity)
#      elif data[i].find("transition velocity dipole") != -1:
#        atom = True
#        j = 0
#        while atom == True :
#          tempdip = data[i+2+j].split()
#          if len(tempdip) != 6 :
#            atom = False
#          else:
#            dipovel.append(map(float,tempdip[1:4]))
#            j = j + 1
#  
#      # Extracts the magnetic moment
#      elif data[i].find("transition magnetic dipole") != -1:
#        atom = True
#        j = 0
#        while atom == True :
#          tempdip = data[i+2+j].split()
#          if len(tempdip) != 4 :
#            atom = False
#          else:
#            mag.append(map(float,tempdip[1:4]))
#            j = j + 1
#  
#      # Extracts the rotational strenght
#      elif data[i].find("R(velocity)") != -1:
#        atom = True
#        j = 0
#        while atom == True :
#          temprot = data[i+1+j].split()
#          if len(temprot) != 6 :
#            atom = False
#          else:
#            rotstr.append(float(temprot[4]))
#            j = j + 1
#  
#  
#      # Extracts the site energy values
#      elif data[i].find("nm ") != -1:
#        site.append(float(data[i].split()[4]))
#  
#      # Compute the center of the chromophore
#      elif data[i].find("NActive=") != -1:
#        NAtoms.append(data[i].split("NActive=")[1].split()[0])
#      
#      i += 1
#  
#    NAtoms = map(int,NAtoms)[1:]
#   
#  
#  
#    return (anum,xyz,site,dipo,dipovel,mag,rotstr)
#  
#  
#  # *****************************************************************************
#  #
#  # Write an output file to be processed with spectrum.py 
#  #
#  def exspect(logfile,site,dipolen,rotstr):
#  
#    PrefixOutFile = logfile.split(".log")[0]
#    OutFile = "%s.out" % PrefixOutFile
#    Out = open(OutFile,'w')
#    N = len(site)
#    for i in range(N):
#      dipo2 = (np.linalg.norm(dipolen[i]))**2
#      Out.write("%3d %12.4f %12.4f %12.4f %12.4f\n" % (i+1,site[i],dipo2,rotstr[i],rotstr[i]))
#    Out.close()
#  
#  




#  
#  # *****************************************************************************
#  
#  def calccom(Z,XYZ):
#  
#    dat     = {  1  : 1.00797,
#                 5  : 11.0093053,
#                 6  : 12.01115,
#                 7  : 14.0067,
#                 8  : 15.9994,
#                 9  : 18.9984033,
#                12  : 24.312,
#                15  : 30.973762,
#                16  : 32.06}
#  
#    C = np.array(XYZ)
#    N = len(C)
#    M = np.zeros(N)
#  
#    k = 0
#    for i in Z:  M[k] = dat[i] ; k+=1
#  
#    COM = np.average(C,axis=0,weights=M)
#  
#    return COM
#  
#  
#  def keepreference(IDRef,XYZ):
#  
#    IDRef = np.array(IDRef)-1
#    if len(IDRef) != 2 :
#      print "You have to specify at least 2 atoms"
#      sys.exit()
#    ID1 = IDRef[1]  
#    ID2 = IDRef[0]
#    Ref = np.array(XYZ[ID1])-np.array(XYZ[ID2])
#    RefVer = Ref/np.linalg.norm(Ref)
#  
#    return RefVer
   
   
#****************************************************
# Compute the rotational strenght in the velocity form
# Gaussian util : ecdten.F
#
# EEN,Del,RxDel can be floats or np.arrays
def RVelIso(EEN,Del,RxDel):

  FactRV = 471.4436078822227
  Fact   = -FactRV/EEN
  RVel   = np.sum(Del*RxDel,axis=-1)
  Angle  = np.arccos(RVel/(np.linalg.norm(Del,axis=-1)*np.linalg.norm(RxDel,axis=-1)))
  RVel   = RVel*Fact/2
  Angle  = np.degrees(Angle)

  return RVel,Angle



# *****************************************************************************

#
# If called by command line, read the Gaussian log file provided
#

if __name__ == "__main__" : 

  # Read input line
  parser = arg.ArgumentParser(description="Read Gaussian .log file and extract some useful data...")
  parser.add_argument('--version','-V',action='version',version=c.PROGVERS,
                      help="Show program's full version number, and exit")  

  parser.add_argument('-v',help='Increase the verbosity of the output',action="count",default=0)
  parser.add_argument('logfile',help='Gaussian output log file to process')
  parser.add_argument('-o','--out',help='Prefix for output files')
  parser.add_argument('--seltran',help='Activate seltran keyword',action="store_true")
  parser.add_argument('--cent',help='How to compute the chormophore center',nargs="*",default="geom")
  parser.add_argument('--anadipo',help='Analysis of the transition dipole orientation with respect to a certain axis')
  parser.add_argument('--mu',help='Prefix for output files',action="store_true")
  parser.add_argument('--ang',help='Compute angle',action="store_true")
  parser.add_argument('--spec',help='Extract data for produce uv and cd spectrum',action="store_true")
  parser.add_argument('--reorient',help='Indicate 2 atoms as reference direction',nargs="*",type=int,default=[1,2])
  args = parser.parse_args()

  c.welcome()

  # Set options
  InLogFile   = args.logfile
  if args.out:
    c.OPT['OutPrefix'] = args.out
    c.setoutfiles()

  c.OPT['verbosity'] = args.v
  c.OPT['seltran']   = args.seltran
  c.OPT['reorient']  = args.reorient
  c.OPT['Cent']      = args.cent
  if args.anadipo  :
    c.OPT['anadipo'] = True
    c.ExtFiles['refaxis'] = args.anadipo


  # Check if the .log files
  c.checkfile(InLogFile)

  # Read Gaussian Log files
  out = r.readgaulog(InLogFile)
  anum,xyz,site,dipo,dipovel,mag,NAtom,NTran,NChrom = r.readgaulog(InLogFile)

  c.anum   = anum
  c.xyz    = xyz
  c.NChrom = NChrom
  c.NAtom  = [NAtom]
  c.NTran  = [NTran]

  dipo     = np.array(dipo)
  dipovel  = np.array(dipovel)
  mag      = np.array(mag)
  site     = np.array(site)

  # Tweak for couplings
  coup = np.zeros(c.NTran)

  # Compute the chromophore center
  Cent = u.calchromcent()


  # Select transtion (if requested by user)
  #  if c.OPT['seltran'] == True:
  #    Coup = 0.0
  #    Site,Dipo,DipoVel,Mag,Cent,NewCoup = readdata.seltran(site,dipo,dipovel,mag,Cent,Coup)

  # Analyze electric transition dipole moment
  c.ChromList = [InLogFile.split('.')[0]]
  u.reorientdipo(dipo,dipovel,mag,coup)


# Compute internal magnetic moment (DipoMag is gauge-dependent)
  Cent = np.array(Cent)
  RxDel = np.cross(Cent/c.PhyCon['ToAng'],dipovel)
  MagInt = mag - RxDel
  MagTot = mag
  MagExt = MagTot-MagInt


  print '\n'
  print " Electric transition dipoles: "
  for i in range(c.NTran[0]):
    NormDip =  np.linalg.norm(dipo[i])
    print " Trans %2d  %8.4f  mu  = %8.4f %8.4f %8.4f  ->  %8.4f " % (i+1,site[i],dipo[i][0],dipo[i][1],dipo[i][2],NormDip)

  u.dipoanalysis(dipo)

  print '\n'
  print " Magnetic transition dipoles (intrinsic): "
  for i in range(c.NTran[0]):
    NormMagInt =  np.linalg.norm(MagInt[i])
    print " Trans %2d  %8.4f  mag = %8.4f %8.4f %8.4f  ->  %8.4f " % (i+1,site[i],MagInt[i][0],MagInt[i][1],MagInt[i][2],NormMagInt)

  # Dipole analysis with internal magnetic moment
  u.dipoanalysis(MagInt)

# Save visudipo
  u.savegeom()
  u.savevisudipo(Cent,dipo,dipovel,mag)

# Recompute rotatory strenght
  R,angle = RVelIso(site/c.PhyCon['ToeV'],dipovel,-MagTot)
  R = np.array(R)

# Save results.out
  SqDip = np.sum(dipo**2,axis=1)
  u.resout(site*c.PhyCon['eV2wn'],SqDip,site*0.0,R)

  print("\n Done! \n")
  

