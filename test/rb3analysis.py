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
} ;");


ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gStyle.SetOptFit(1)


print "------------------------------------------------"
print "INITIALIZING"

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
            if ((chambID[1]).find("RB3")) > -1 :
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
rfile_in = ROOT.TFile.Open("rb3RollsEffFit.root")
fittree = rfile_in.Get("T")

rpcfit = RpcRollEffFit()
fittree.SetBranchAddress("rpcfit", AddressOf(rpcfit, "rollidtree"))

nroll_sel = 0

roll_sel = []

for i in xrange(fittree.GetEntries()):
    fittree.GetEntry(i)
    if rpcfit.chi2 > 10000 :
        print rpcfit.rollidtree
        nroll_sel = nroll_sel + 1
        roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
print nroll_sel

outplotdir = "rb3plots"

if not (os.path.exists(outplotdir)) :
    os.mkdir(outplotdir)

for rollname in roll_sel :
    myh = ROOT.TGraphErrors()
    rfile_in.GetObject("eff_"+rollname,myh)
    myh.Draw("AP")
    c1.SaveAs(outplotdir+"/"+rollname+".png")


print "DONE"







