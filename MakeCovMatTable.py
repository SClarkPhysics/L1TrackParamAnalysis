import uproot 
import numpy as np
import ROOT
import pandas as pd
import awkward as ak
from ROOT import *

"""
This script takes in a ROOT file from L2Trigger/TrackFindingTracklet/test/L1TrackNtupleMaker_cfg.py

First it "unpacks" the root file by taking each event with multiple tracks and writing each track as a row in a Pandas Dataframe, then saving as a CSV

The CSV is then read in as a Pandas Dataframe. Simple analysis is performed to remove any infinities, define helpful variables (Hit Pattern, t, deltas). 
Then for the correlated variables, cuts are made based on eta and hitpattern, thena  2D histogram is filled, a profile plot is created and fit. The slope of the fit is then extracted and saved in a
csv file as the covariance matrix element. 

The output file CSV_Files/matrix.txt then contains the matrix error elements for the correlated and uncorrelated pairs separated by the eta/hp cuts.

2D Histograms and Fit Profile plots are all saved in CovPlots/
"""

gROOT.SetBatch()
RDF = ROOT.ROOT.RDataFrame

#Filename from NTuplemaker to run on
inFile = "SingleMu_pt2to100_100k_D49.root" #Full path to input file

fname = inFile.split("/")[-1] #filename only
fname = fname[0:fname.find(".root")]

#Toggle between trk_variable and tp_variable.  
tktp = "trk_"
if(tktp == "trk_"): mtch="trk_matchtp_"  #trk_vars match up with trk_matchtp_vars
elif(tktp == "tp_"): mtch="matchtrk_"

def MakeFolder(N):
    import os
    if not os.path.exists(N):
     os.makedirs(N)

def unpackFile(filename):
  """
  Function to take root file from ntuplemaker and unpack it so that each track takes one row in a dataframe.
  """
  #Get event tree from file
  chain = TChain("L1TrackNtuple/eventTree")
  chain.Add(inFile)

  fname = filename.split("/")[-1]
  fname = fname[:fname.find(".root")]

  #Create RDataFrame
  Rdf = RDF(chain)


  #Track parameters we are interested in
  keepnames = ['pt','eta','phi','d0','z0']
  trk_names = [tktp + t for t in keepnames] #add the prefix to each name
  if(tktp == "tp_"): trk_names.append("tp_charge")
  match_names = [mtch+t for t in keepnames] #Get list of match names

  want_names = trk_names + match_names
  #Extra variables we want
  want_names.append("trk_hitpattern")
  want_names.append("trk_nstub")
  want_names.append("trk_seed")
  want_names.append("trk_matchtp_pdgid")

  #Convert RDataFrame to Numpy arrays
  events = Rdf.AsNumpy(want_names)
  #Convert Numpy arrays to Pandas dataframe
  odf = pd.DataFrame(events)

  newdict = {}
  for cname in odf.columns:
    newarray = []
    for vv in odf[cname]:
      for dd in vv:
        if(cname in ["trk_hitpattern","trk_nstub", "trk_matchtp_pdgid"]):
          newval = int(dd)
        else:
          newval = float(dd)
        newarray.append(newval)

    newdict[cname] = newarray

  df = pd.DataFrame(newdict)

  #Removing rows with infinite values
  df.replace([np.inf, -np.inf], np.nan, inplace=True)
  df = df[df['{}pt'.format(tktp)].notna()]

  MakeFolder("UnpackedFiles")
  df.to_csv("UnpackedFiles/{}.csv".format(fname))
  print("Saving unpacked file as: {}".format("UnpackedFiles/{}.csv".format(fname)))
  
  return

def getTheta(eta):
    theta = 2*np.arctan( np.exp( -1 * eta))
    if (theta > np.pi/2):
        theta = theta - np.pi
        
    return theta

def getLambda(eta):
    return np.pi/2 - getTheta(eta)

def getT(eta):
    return np.tan(getLambda(eta))

def getEta(t):
    return -1*np.log(np.tan(np.pi/4 - 1/2 * np.arctan(t)))


#Open file and unpack
print("Opening Root file: {}".format(inFile))
print("Unpacking file")
unpackFile(inFile)
print("File Unpacked")

#Read unpacked file
df = pd.read_csv("UnpackedFiles/{}.csv".format(fname), index_col=0)

df["{}t".format(tktp)] = df["{}eta".format(tktp)].apply(getT)
df["{}t".format(mtch)] = df["{}eta".format(mtch)].apply(getT)

#Names of variables to keep
keepnames = ['pt','eta','phi','d0','z0']

df["{}invpt".format(tktp)] = 1/df["{}pt".format(tktp)]
df["{}invpt".format(mtch)] = 1/df["{}pt".format(mtch)]
keepnames.append("invpt")

keepnames = ['t' if item=='eta' else item for item in keepnames]

for n in keepnames:
    tname = tktp + n
    mname = mtch+n
    nmname = "match_"+n
    dname = "delta_"+n
    df[dname] = df[mname]-df[tname]
    df.rename(columns = {mname:nmname}, inplace = True)

#Remove tracks with nonsensical delta
df = df[df.match_pt > -999]

# ## Get Hit Pattern in Binary

def padded_bin(i):
    s = bin(i)
    return  s[2:].zfill(7)

def getBit(i,bb):
    return int(padded_bin(i)[bb])

def sumHit(mystring):
    mysum=0
    for i in mystring:
        mysum += int(i)
    return mysum

df["BiHitPattern"] = df["trk_hitpattern"].apply(padded_bin)
pd.options.display.float_format = '{:7}'.format

#####
# b6 is INNERMOST hit
# b1 is OUTER
# b0 is "extra" outermost bit, can only be ON for 1.68 < eta < 2.08

for ii in range(len(df["BiHitPattern"][0])):
    df["b{}".format(ii)] = df["trk_hitpattern"].apply(getBit, args=[ii])

df["n_hits"] = sum([df["b{}".format(i)] for i in range(0,7)])
df["nps_hits"] = sum([df["b{}".format(i)] for i in range(4,7)])

# ## Get $\mu^{+/-}$ Only
pos_df = df[df.trk_matchtp_pdgid > 0]

#####################################################


plist = ["invpt", "t", "phi", "d0", "z0"]

#Maxes contains the |max| values to fill the histograms for each delta_variable
maxes = {"delta_invpt": 0.005, "delta_t": 0.005, "delta_phi": 0.001, "delta_d0": 0.03, "delta_z0": 0.5}
#Plot titles
ltitles = {"invpt": "p_{T}^{-1}",
           "phi": "#phi",
           "t": "t",
           "d0": "d_{0}",
           "z0": "z_{0}",
          }

def smartCutter(dframe, arg):
    """
    Function to apply the eta and hit pattern cuts
    Inputs pandas dataframe and list of cuts separated by comma
    output is dataframe with cuts implemented
    """
    
    myframe = dframe
    if arg=="": 
        print("No cuts")
        print(len(dframe))
        return dframe
    
    for ag in arg.split(","):
#         print(ag)
        
        if(ag[0]==" "): ag = ag[1:]

        col = ag.split(" ")[0]
        val = ag.split(" ")[1]
    
        if col[0]=="b":
            myframe = myframe.loc[myframe[col] == int(val)]
        elif col=="nps_hits":
            myframe = myframe.loc[myframe[col] == int(val)]
        elif col=="etal":
            myframe = myframe.loc[abs(myframe["trk_eta"]) < float(val)]
        elif col=="etag":
            myframe = myframe.loc[abs(myframe["trk_eta"]) > float(val)]
        else: print("Cut Not Recognized")
        
        print(col, val)
        print(len(myframe))
    print("_______________")
        
    return myframe

def makeTH2(dframe, xcol, ycol, sname):
  """
  Function to make and return a TH2D and Profile plot
  """
  
  mx = maxes["delta_{}".format(xcol)]
  my = maxes["delta_{}".format(ycol)]
  print("Max x: {}".format(mx))
  print("Max y: {}".format(my))
  ltx, lty = ltitles[xcol], ltitles[ycol]
  th2 = ROOT.TH2D("{}_{}_{}".format(xcol, ycol, sname),"#Delta {} vs. #Delta {}; #Delta {}; #Delta {}".format(lty,ltx,ltx,lty), nbins, -mx, mx, nbins, -my, my)

  prf=ROOT.TProfile("prof","Profile",nbins,-my,my,-my,my)

  for xx, yy in zip(dframe["delta_{}".format(xcol)], dframe["delta_{}".format(ycol)]):
    th2.Fill(xx,yy)
    prf.Fill(yy,xx)

  return th2.Clone(), prf.Clone()

def AddLineText(pad, slope, intercept):
  """
  Add slope and intercept text to the plots
  """
  sText = "Slope = {:.3f}".format(slope)
  iText = "Inter = {:.3E}".format(intercept)
  tsize = 0.5
  TextFont   = 61  
  TextSize      = 0.5
  TextOffset    = 0.15
  H = pad.GetWh()
  W = pad.GetWw()
  l = pad.GetLeftMargin()
  t = pad.GetTopMargin()
  r = pad.GetRightMargin()
  b = pad.GetBottomMargin()
  e = 0.025
  pad.cd()
  latex = TLatex()
  latex.SetNDC()
  latex.SetTextAngle(0)
  latex.SetTextColor(kRed)	
  extraTextSize = 0.76*TextSize
  pad.cd()
  latex.SetTextFont(TextFont)
  latex.SetTextSize(TextSize*t)
  latex.SetTextAlign(11)
  latex.DrawLatex(l + 0.025, b + 0.08, sText)
  latex.DrawLatex(l + 0.025, b + 0.03, iText)
  pad.Update()

def AddNTrkText(pad, ntrks):
  """
  Add Number of Tracks text to each plot
  """
  trkText = "{} Tracks".format(ntrks)
  tsize = 0.5
  TextFont   = 61  
  TextSize      = 0.5
  TextOffset    = 0.15
  H = pad.GetWh()
  W = pad.GetWw()
  l = pad.GetLeftMargin()
  t = pad.GetTopMargin()
  r = pad.GetRightMargin()
  b = pad.GetBottomMargin()
  pad.cd()
  latex = TLatex()
  latex.SetNDC()
  latex.SetTextAngle(0)
  latex.SetTextColor(kRed)	
  extraTextSize = 0.76*TextSize
  pad.cd()
  latex.SetTextFont(TextFont)
  latex.SetTextSize(TextSize*t)
  latex.SetTextAlign(11)
  latex.DrawLatex(0.68,0.850, trkText)
  pad.Update()

#####################################################

# ## Make ROOT TH2's of the stuff I want

#Number of bins to use for the TH2D
nbins = 30

def GetNEmptyBins(hist):
  """
  Calculate number of empty bins in a histogram
  """
  ne = 0
  for bb in range(hist.GetNbinsX()):
    bc = hist.GetBinContent(bb)
    if(bc == 0):
      ne += 1
  return ne

def GetElementCorr(xcol, ycol, scut, tcut):
    """
    Function to get the Covariance Matrix element for a set of correlated variables
    Inputs variables as xcol, ycol
    scut is string with cuts to perform
    tcut is cuts formatted as text for plot titles
    """

    #Make output folder if it doesn't exist
    PlotFolder = "CovPlots/{}_{}/".format(xcol,ycol)
    MakeFolder(PlotFolder)

    #Make cuts on dataframe, save as new dataframe
    plotdf = smartCutter(pos_df, scut)
    
    #Text formatting for cuts
    if(scut[0]==''): cut_text = ", No #eta, HP Cuts"
    else:
      sel,seh,sb1,sb6=tcut[0],tcut[1],tcut[2],tcut[3]
      cut_text = ", {} < #eta < {} , HP = *{}****{}".format(sel,seh,sb1,sb6)
      if(seh==9999): cut_text = ", #eta > {} , HP = *{}****{}".format(sel,sb1,sb6)
      
    savename = scut.replace(" ","_")
    savename = savename.replace(",","_")
    if(scut==""): savename = "nocuts"

    #Get 2D hist and profile plot
    myhist,prof = makeTH2(plotdf, xcol, ycol, savename)
    myhist.SetTitle(myhist.GetTitle() + cut_text)
    myhist.SetStats(0)
    
    c1 = ROOT.TCanvas()
    c1.cd()
    myhist.Draw("colz")
    ntrk = int(myhist.GetEntries())
    AddNTrkText(gPad, ntrk)

    ##################

    #Use root to make the profile plot, better than the manually created one
    py = myhist.ProfileY() #This one is better
    py.SetName(myhist.GetName() + "Profile Y")

    #Establish fit range from maxes
    fmy = maxes["delta_{}".format(ycol)]

    #Fit the profile plot to a line
    cy = ROOT.TCanvas()
    cy.cd()
    myfunc = "pol1(0)"
    fit = TF1("linefit", myfunc, -fmy, fmy)
    py.Fit(fit, "EMR0")
    nempty = GetNEmptyBins(py)
    #If many empty bins, plot with different options
    if(nempty <= py.GetNbinsX()/3 + 1): py.Fit(fit, "EMR0")
    else: py.Fit(fit, "WR0")

    #Plot formatting and saving
    fit.SetLineColor(ROOT.kRed)
    fit.SetLineWidth(2)
    py.SetLineWidth(2)

    intercept = fit.GetParameter(0)
    slope = fit.GetParameter(1)

    py.SetStats(0)
    py.Draw("E0")
    fit.Draw("same")
    AddLineText(gPad, slope, intercept)
    AddNTrkText(gPad, ntrk)

    c1.Print(PlotFolder+savename+".png")
    cy.Print(PlotFolder+savename+"_profY.png")

    return slope

kf_vars = ["t","z0","phi","invpt","d0"]

#Pairs of uncorrelated variables
uncorr_pairs = [
                ("phi","t"),
                ("phi","z0"),
                ("invpt","t"),
                ("invpt","z0"),
                ("d0","t"),
                ("d0","z0"),
                ]

#Pairs of correlated variables
corr_pairs = [
               ("invpt","d0"),
               ("invpt","phi"),
               ("d0","phi"),
               ("z0","t")
              ]

#Establish cuts
eta_cuts = [0.0, 0.9, 1.26, 1.68, 2.08, 9999] #Eta regions
hp_cuts = ["b1 1, b6 1", "b1 1, b6 0", "b1 0, b6 1", "b1 0, b6 0"] #HP Cuts. b1 outer, b6 inner

#Open output folder
CSVFolder = "CSV_Files/"
MakeFolder(CSVFolder)
outFileName="{}matrix.txt".format(CSVFolder)
CovFile=open(outFileName,"w")

for idx_eta in range(len(eta_cuts)-1):
  #if(idx_eta > 0): break #For debugging
  eta_cut = eta_cuts[idx_eta]
  hcount = 0
  for hp_cut in hp_cuts:
    #if (hcount > 0): break #For Debugging
    hcount += 1
    b1=hp_cut.split(",")[0][-1]
    b6=hp_cut.split(",")[1][-1]
    CovFile.write("eta{}_{},hp{}_{}\n".format(eta_cut, eta_cuts[idx_eta+1],b1,b6))

    scut="etag {}, etal {}, {}".format(eta_cut, eta_cuts[idx_eta+1], hp_cut)
    tcut=[eta_cut, eta_cuts[idx_eta+1], b1,b6] 

    for v1 in kf_vars:
      for v2 in kf_vars:
        cov_element=0

        #Diagonal Elements
        if(v1==v2):
          cov_element=1

        #Uncorrelated pairs, covariance=0
        elif( (v1,v2) in uncorr_pairs or (v2,v1) in uncorr_pairs ):
          cov_element=0

        #Correlated pairs, compute covariance from slope of profile
        elif( (v1,v2) in corr_pairs):
          cov_element=GetElementCorr(v1,v2,scut,tcut)
        elif( (v2,v1) in corr_pairs):
          cov_element=GetElementCorr(v2,v1,scut,tcut)

        CovFile.write("{},{}:{}\n".format(v1,v2,cov_element))

print("Output Saved as: {}".format(outFileName))
