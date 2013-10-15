#!/usr/bin/python

import sys
#sys.argv.append('-b')
import os, commands
import math
import ROOT

from ROOT import *

gROOT.ProcessLine(
    "struct RpcRollEffFit {\
    Int_t rollidtree;\
    Float_t p0;\
    Float_t p0err;\
    Float_t p1;\
    Float_t p1err;\
    Float_t chi2;\
    Float_t ndof;\
    Float_t refp;\
    Float_t refperr;\
    Float_t deltahv;\
    Float_t deltaslope;\
} ;");


ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gStyle.SetOptFit(1)


print "------------------------------------------------"
print "INITIALIZING"

# ARGUMENTS
# python rollFitHistory.py eff RB3 -> fit eff plots no normalization 
# python rollFitHistory.py eff RB4 norm -> fit eff plots with normalization 

variab = sys.argv[1]
variabCap = variab.capitalize()
if variab[0].isupper() :
    variab = variab.lower()

isbarrel = False
isendcap = False
rollstring = sys.argv[2]
if rollstring.find("RB") != -1  :
    isbarrel = True
else :
    isendcap = True


norm = ""
postfix = ""
pf1 = ""
pf2 = ""
donorm = False

if len(sys.argv) > 3 :
    if sys.argv[3] == "norm" :
        donorm = True
        norm = sys.argv[3]
        postfix = "AllRB3"
        pf1 = "_"+norm
        pf2 = "_"+(norm.capitalize())

# fill the list of bad chambers from input blacklist
badChambers = []
blackListFileName = "blacklist.txt"
blackListFile = open(blackListFileName, "r")
for chamber in blackListFile:
    chamber = (chamber.rstrip().split(' '))[1]
#    print chamber
    badChambers.append(chamber)

# dictionaries ID -> string
barrelRollsDict = {}
endcapRollsDict = {}
rb3RollsDict = {}
re1Ring3RollsDict = {}

# dictionary with list of RB3 ID for each Wheel
wheels_names = ["W-2", "W-1", "W+0", "W+1", "W+2" ]
rb3Rolls_InW = {}
for w in wheels_names :
    rb3Rolls_InW[w] = []

# open chamberID definition file and create a dictionary
#chambDict = {}
chambCounter = 0
chambDictFile = open("mapRoll.ascii", "r")
for chambEntry in chambDictFile:
#    print chambEntry
    chambCounter = chambCounter + 1
    chambID = chambEntry.rstrip().split(' ')
#    print chambID[0], chambID[1]
    if not( (chambID[1]) in badChambers ) :
        if (chambID[1]).startswith("W"):
            barrelRollsDict[ chambID[0] ] = chambID[1]
            if ((chambID[1]).find(rollstring)) > -1 :
                rb3RollsDict[ chambID[0] ] = chambID[1]
                if (chambID[1]).startswith("W-2"):
                    (rb3Rolls_InW[ "W-2" ]).append( chambID[0] )
                elif (chambID[1]).startswith("W-1"):
                    (rb3Rolls_InW[ "W-1" ]).append( chambID[0] )
                elif (chambID[1]).startswith("W+0"):
                    (rb3Rolls_InW[ "W+0" ]).append( chambID[0] )
                elif (chambID[1]).startswith("W+1"):
                    (rb3Rolls_InW[ "W+1" ]).append( chambID[0] )
                else:
                    (rb3Rolls_InW[ "W+2" ]).append( chambID[0] )
        else :
            endcapRollsDict[ chambID[0] ] = chambID[1]
            if (((chambID[1]).find("RE+1_R3")) > -1) or (((chambID[1]).find("RE-1_R3")) > -1)  :
                re1Ring3RollsDict[ chambID[0] ] = chambID[1]

#    else :
#        print " BAD ROLL FOUND !!!!!!!!!!!!!!"

rb3ids = rb3RollsDict.keys()
re1Ring3ids = re1Ring3RollsDict.keys()
barrelids = barrelRollsDict.keys()
endcapids = endcapRollsDict.keys()

# NOW WE HAVE THE LISTS OF BARREL / ENDCAP / RB3_ROLLS (reference) / RE1Ring3_ROLLS (reference)

# deltaHV50 and deltaslope info for rolls
dhv50_dict = {}
dslope_dict = {}
plateauInfoFile = open("deltaHV_2012.txt", "r")
for p_entry in plateauInfoFile:
    readline = p_entry.rstrip().split(' ')
#    print readline[1], readline[2], readline[3]
    dhv50_dict[ readline[1] ] = readline[2]
    dslope_dict[ readline[1] ] = readline[3]


# CREATE the run-date dictionary
runs_time = {}
runTimeFile = open("rpcrun_epoch.txt", "r")
for r_entry in runTimeFile:
    ent = r_entry.rstrip().split('\t')
    runs_time[ ent[0] ] = ent[1]
#print runs_time

# GET THE RUN LIST
# open the official input run list from a file and make a vector from it
rpcRuns = []
rpcRunsTime = []

runListFile = open("runslist.txt", "r")
for eachRun in runListFile:
    rpcRuns.append(int(eachRun.rstrip()))
    rpcRunsTime.append( float(runs_time[ eachRun.rstrip() ]) )

ntotruns = len(rpcRuns) 

############################################################

firstRun = "190679"
lastRun = "209151"

print firstRun, runs_time[firstRun]
print lastRun, runs_time[lastRun]

#runNumIndexFirst = rpcRuns.index( firstRun )
#runNumIndexLast = rpcRuns.index( lastRun )

print "OPENING HISTORY PLOTS"

# open root file
rfile_in = ROOT.TFile.Open(variab+"_history_Graph_nofit_newnorm.root")

outfilename = rollstring+"Rolls"+variabCap+(norm.capitalize())+"Fit.root"

rfile_out = ROOT.TFile.Open(outfilename,"RECREATE")

rpcfit = RpcRollEffFit()

mytree = ROOT.TTree('T', 'T')
mytree.Branch('rpcfit',rpcfit,'rollidtree/I:p0/F:p0err/F:p1/F:p1err/F:chi2/F:ndof/F:refp/F:refperr/F:deltahv/F:deltaslope/F')

for myids in rb3ids :
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')
    h1_name = variab+"_"+tmpstring

    h_name = variab+"history"+pf1+"/"+h1_name+pf1
    myh = ROOT.TGraphErrors()
    rfile_in.GetObject(h_name,myh)

    mchi2 = 0
    mndf = 0
    mp0 = 0
    mp0err = 0
    mp1 = 0
    mp1err = 0
    mrefp = 0
    mrefperr = 0
    
    tf1 = ROOT.TF1("tf1","[0]+x*[1]")
    fitstat = myh.Fit(tf1,"S","",2012.26891223,2012.95816911)

    mchi2 = fitstat.Chi2()
    mndf = fitstat.Ndf()
    mp0 = fitstat.Parameter(0)
    mp0err = fitstat.ParError(0)
    mp1 = fitstat.Parameter(1)
    mp1err = fitstat.ParError(1)

    mcov = (fitstat.GetCovarianceMatrix())(0,1);     

# time ref.set to 1st Run 2012

    mrefp = mp0 + 2012.26891223 * mp1

    errsqref1 = mp0err*mp0err
    errsqref2 = 2012.26891223 * 2012.26891223 * mp1err * mp1err
    covtermref = 2.0 * 2012.26891223 * mcov
    mrefperr = math.sqrt(errsqref1 + errsqref2 + covtermref)

        
    print "roll = ", tmpstring
    print "fit chi2 = ", mchi2, " / ", mndf
    print mp0, " +/- ", mp0err
    print mp1, " +/- ", mp1err
    print "cov =", mcov
    print variabCap, " @ 1st Run 2012 =", mrefp, " +/- ", mrefperr

    rpcfit.rollidtree = int(myids)
    rpcfit.p0 = mp0
    rpcfit.p0err = mp0err
    rpcfit.p1 = mp1
    rpcfit.p1err = mp1err
    rpcfit.chi2 = mchi2
    rpcfit.ndof = mndf
    rpcfit.refp = mrefp
    rpcfit.refperr = mrefperr
    rpcfit.deltahv = float(dhv50_dict.get(myids,-10000))
    rpcfit.deltaslope = float(dslope_dict.get(myids,-10000))

    mytree.Fill()
    myh.GetXaxis().SetLimits(2012.2,2013)
    myh.Write()


if donorm :
    myh = ROOT.TGraphErrors()
    myh1 = ROOT.TGraphErrors()
    if isbarrel :
        rfile_in.GetObject(variab+"history"+pf1+"/"+"rb3ref_eff",myh)
        rfile_in.GetObject(variab+"history"+pf1+"/"+"rb3ref_ngood",myh1)
    else : 
        rfile_in.GetObject(variab+"history"+pf1+"/"+"re1ring3ref_eff",myh)
        rfile_in.GetObject(variab+"history"+pf1+"/"+"re1ring3ref_ngood",myh1)

    tf1 = ROOT.TF1("tf1","[0]+x*[1]")
    fitstat = myh.Fit(tf1,"S","",2012.26891223,2012.95816911)
    myh.GetXaxis().SetLimits(2012.2,2013)
    myh.Write()
    myh1.GetXaxis().SetLimits(2012.2,2013)
    myh1.Write()


mytree.Print()
mytree.Write()

rfile_out.Close()

print "DONE"
print "variable considered for fit: ", variab
print "SELECTED ROLLS: ", rollstring
print "normalization ?", donorm








