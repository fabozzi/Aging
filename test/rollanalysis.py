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
# python rollanalysis.py eff RB3 [inbox/outofbox/largechi/chi2gt10k/outofbox_negp1/inbox_largeNegDeltahv/...] 
# python rollFitHistory.py eff RB4 outofbox norm -> fit eff plots with normalization 

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

chi2cut = 10000
if rollstring == "RB2" :
    chi2cut = 4000

cut_dict = {}

cut_dict["inbox"] = 0
cut_dict["outofbox"] = 1
cut_dict["largechi"] = 2
cut_dict["chi2gt10k"] = 3
cut_dict["outofbox_negp1"] = 4
cut_dict["inbox_largeNegDeltahv"] = 5
cut_dict["chi2lt10k_p1gt01"] = 6
cut_dict["inbox_dhvgt004"] = 7
cut_dict["inbox_dslltm30"] = 8
cut_dict["reference"] = 9

selregion = sys.argv[3]

norm = ""
pf1 = ""
pf2 = ""
donorm = False

if len(sys.argv) > 4 :
    if sys.argv[4] == "norm" :
        donorm = True
        norm = sys.argv[4]
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
infilename = rollstring+"Rolls"+variabCap+(norm.capitalize())+"Fit.root"

rfile_in = ROOT.TFile.Open(infilename)
fittree = rfile_in.Get("T")

rpcfit = RpcRollEffFit()
fittree.SetBranchAddress("rpcfit", AddressOf(rpcfit, "rollidtree"))

nroll_sel = 0

roll_sel = []


selfilename = rollstring+"_"+selregion+".txt"
f = open(selfilename, 'w')

############################################################################

print rollstring, selregion, cut_dict[selregion]

for i in xrange(fittree.GetEntries()):
    fittree.GetEntry(i)
    if cut_dict[selregion] == 0 : # this is the inbox cut
        if (rpcfit.chi2 < chi2cut) and (math.fabs(rpcfit.p1)<0.02) and (rpcfit.refp>0.9)  :
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 1 : # this is the outofbox cut
        if (rpcfit.chi2 < chi2cut) and ( (math.fabs(rpcfit.p1)>0.02) or (rpcfit.refp<0.9) ) :
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 2 : # this is the chi2 > chi2cut but < 10000
        if (rpcfit.chi2 > chi2cut) and (rpcfit.chi2 < 10000):
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 3 : # this is the chi2 > 10000
        if (rpcfit.chi2 > 10000) :
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 4 : # this is the outofboxcut and only negative slopes
        if (rpcfit.chi2 < chi2cut) and ( (math.fabs(rpcfit.p1)>0.02) or (rpcfit.refp<0.9) ) and (rpcfit.p1 < - 0.02) :
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 5 : # this is the inboxcut but large negative deltahv
        if (rpcfit.chi2 < chi2cut) and (math.fabs(rpcfit.p1)<0.02) and (rpcfit.refp>0.9) and (rpcfit.deltahv < -0.5 ):
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 6 : # this is the outlier positive p1 cut 
        if (rpcfit.chi2 < chi2cut) and (rpcfit.p1>0.1):
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 7 : # this is the outlier region for deltahv within loose inbox
        if (rpcfit.chi2 < chi2cut) and (math.fabs(rpcfit.p1)<0.05) and (rpcfit.refp>0.9) and (rpcfit.deltahv > -0.5) and (rpcfit.deltahv > 0.04 ) :
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 8 : # this is the outlier region for deltaslope within loose inbox
        if (rpcfit.chi2 < chi2cut) and (math.fabs(rpcfit.p1)<0.05) and (rpcfit.refp>0.9) and (rpcfit.deltahv > -0.5) and (rpcfit.deltaslope < -30 ) :
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( rb3RollsDict[ str( rpcfit.rollidtree ) ] + "\t" + str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )
    elif cut_dict[selregion] == 9 : # this is the cut defining reference rolls
        if (rpcfit.chi2 < chi2cut) and (rpcfit.p1 > -0.02) :
            print rpcfit.rollidtree, rb3RollsDict[ str( rpcfit.rollidtree ) ]
            f.write( str(rpcfit.rollidtree) +"\n" )
            nroll_sel = nroll_sel + 1
            roll_sel.append( (rb3RollsDict[ str( rpcfit.rollidtree ) ] ).replace('+','p').replace('-','m')  )

f.close()

outplotdir = rollstring+"plots_"+selregion

if not (os.path.exists(outplotdir)) :
    os.mkdir(outplotdir)

for rollname in roll_sel :
    myh = ROOT.TGraphErrors()
    rfile_in.GetObject("eff_"+rollname,myh)
    myh.Draw("AP")
    c1.SaveAs(outplotdir+"/"+rollname+".png")


print "DONE"
print "selected", nroll_sel, "rolls"
print "writing results in", outplotdir






