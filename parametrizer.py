import numpy as np
import ROOT
import sys,os
from array import array

#Feed in root file as command line argument
#Should be formatted as NTuple from L1TrackNtupleMaker

#Command line arguments are input filename, and output file name (without the .root)

infilename=sys.argv[1]

infile = ROOT.TFile(infilename, "r")
tree = infile.Get("L1TrackNtuple/eventTree")

outfilename = "OutputTrees/" + sys.argv[2] + ".root"
outfile = ROOT.TFile(outfilename, "RECREATE")
out_tree = tree.CloneTree(0)

p_slope_values = ROOT.std.vector('float')()
p_pt_values = ROOT.std.vector('float')()
out_tree.Branch("param_slope", p_slope_values)
out_tree.Branch("param_pt", p_pt_values)

def padded_bin(i):
    #Convert decimal to binary, always 7 digit output
    s = bin(i)
    return  s[2:].zfill(7)

def getR(z,eta):
    return z*np.tan( 2*np.arctan( np.exp(-1*eta)))

def ComputeRadius(eta, hp):
    #Compute radius of track from eta and hitpattern
    
    eta = abs(eta)
    rhp = hp[::-1]
    ohit = rhp.rfind("1") #Outermost hit (reverse hit pattern so it goes inner->outer, then find last hit)

    irmean = {"L1":851, "L2":1269, "L3":1784, "L4":2347, "L5":2936, "L6":3697}
    izmean = {"D1":2239, "D2":2645, "D3":3163, "D4":3782, "D5":4523}
    rmaxdisk = 120.0
    zlength = 120.0
    # https://github.com/cms-L1TK/cmssw/blob/L1TK-dev-12_6_0_pre5/L1Trigger/TrackFindingTracklet/interface/Settings.h#L171 

    # https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/L1Trigger/TrackTrigger/src/TrackQuality.cc#L58
    if(eta <= 1.26):
        trklist = ["L1", "L2", "L3", "L4", "L5", "L6"]
        if(ohit==6):ohit=5
    elif(eta > 1.26 and eta <= 1.68):
        trklist = ["L1", "L2", "L3", "D1", "D2", "D3", "D4"]
    elif(eta > 1.68 and eta <= 2.08):
        trklist=["L1", "L2", "D2", "D3", "D4", "D5"]
        if(ohit==6):ohit=5
    elif(eta > 2.08):
        trklist=["L1", "D1", "D2", "D3", "D4", "D5"]
        if(ohit==6):ohit=5

    layer = trklist[ohit]
    if(layer.startswith("L")):
        radius = irmean[layer] * rmaxdisk / 4096
        OL = 1
    else:
        zz = izmean[trklist[ohit]] * zlength / 2048
        radius = getR(zz,eta)
        OL = 0

    return radius

def GetFromParametrization(rr,b1,b6):
  hpn = b1+b6
  pfile = open("ParamFiles/radius.csv","r")
  for lin in pfile.readlines():
    sl = lin.split(",")
    rl=float(sl[0])
    rh=float(sl[1])
    bb=sl[2]
    this_s=float(sl[3])

    if(rr >= rl and rr < rh and hpn==bb):
      return this_s

  print("Problem computing parametrization. Returning False")
  return False

for tcount in range(tree.GetEntries()):
  tree.GetEntry(tcount)

  if(tcount % 10000 == 0): print(tcount)  #Watch the status

  etas = [float(ee) for ee in tree.trk_eta]
  hps_raw = [int(hp) for hp in tree.trk_hitpattern]
  hps_binary = [padded_bin(hp) for hp in hps_raw]
  rads = [ComputeRadius(ee,hh) for (ee,hh) in zip(etas,hps_binary)]
  pts = [float(pp) for pp in tree.trk_pt]
  ipts = [1/pp for pp in pts]
  d0s = [float(dd) for dd in tree.trk_d0]

  p_slopes = []
  p_ipts = []
  for (rr, hh, ipt, dd) in zip(rads, hps_binary, ipts, d0s):
    b1 = hh[1] #Ignore first bit
    b6 = hh[-1]

    slope = GetFromParametrization(rr,b1,b6)
    p_slopes.append(slope)
    p_ipts.append(ipt - slope*dd)

  p_slope_values.clear()
  p_pt_values.clear()
  for (ss,ipp) in zip(p_slopes,p_ipts):
    p_slope_values.push_back(ss)
    p_pt_values.push_back(1/ipp)
  out_tree.Fill()

outfile.Write()
outfile.Close()
infile.Close()
