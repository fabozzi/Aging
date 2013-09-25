#!/usr/bin/python

import sys
#sys.argv.append('-b')
import os, commands
import math
import ROOT

from ROOT import *
gROOT.ProcessLine(
    "struct RpcEffStruct {\
    double run_number;\
    double raw_id;\
    double eff_seg;\
    double eff_seg_error;\
    double n_extrap;\
    double cluster_size;\
    double clus_size_bin01;\
    double clus_size_bin02;\
    double clus_size_bin03;\
    double clus_size_bin04;\
    double clus_size_bin05;\
    double clus_size_bin06;\
    double clus_size_bin07;\
    double clus_size_bin08;\
    double clus_size_bin09;\
    double clus_size_bin10;\
} ;");

gROOT.ProcessLine(
    "struct RpcRollEffFit {\
    Int_t rollid;\
    Float_t p0;\
    Float_t p1;\
    Float_t chi2;\
    Float_t ndof;\
} ;");


ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()        # don't pop up canvases



# fill the list of bad chambers from input blacklist
badChambers = []
blackListFileName = "blacklist.txt"
blackListFile = open(blackListFileName, "r")
for chamber in blackListFile:
    chamber = (chamber.rstrip().split(' '))[1]
#    print chamber
    badChambers.append(chamber)


print "------------------------------------------------"



barrelRollsDict = {}
endcapRollsDict = {}
rb3RollsDict = {}



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
        else :
            endcapRollsDict[ chambID[0] ] = chambID[1]
#    else :
#        print " BAD ROLL FOUND !!!!!!!!!!!!!!"

#print barrelRollsDict
#print len(barrelRollsDict)
#print "------------------------------------------------"
#print endcapRollsDict
#print len(endcapRollsDict)
#print "------------------------------------------------"
#print rb3RollsDict
#print len(rb3RollsDict)
#print "------------------------------------------------"

rb3ids = rb3RollsDict.keys()
barrelids = barrelRollsDict.keys()

# NOW WE HAVE THE LISTS OF BARREL / ENDCAP / RB3_ROLLS (reference)

# GET THE RUN LIST
# open the official input run list from a file and make a vector from it
rpcRuns = []
runListFile = open("runslist.txt", "r")
for eachRun in runListFile:
    rpcRuns.append(int(eachRun.rstrip()))
ntotruns = len(rpcRuns) 
############################################################

# define dictionaries to store (for each roll) list of efficiencies and errors by run

barrelRollsEffDict = {}
barrelRollsErrDict = {}

for y in barrelids :
    barrelRollsEffDict[ y ] = []
    barrelRollsErrDict[ y ] = []
    for j in range(0, ntotruns):
        barrelRollsEffDict[ y ].append(0.0)
        barrelRollsErrDict[ y ].append(0.0)


############################################################

# define vectors to store reference (RB3) efficiencies and errors by run

rb3EffInRun = []
rb3EffErrInRun = []
rb3GoodRollsInRun = []

for y in range(0, ntotruns):
    rb3EffInRun.append(0.0)
    rb3EffErrInRun.append(0.0)
    rb3GoodRollsInRun.append(0)


firstRun = 190679
lastRun = 194076
runNumIndexFirst = rpcRuns.index( firstRun )
runNumIndexLast = rpcRuns.index( lastRun )

# open DB Tree file
rfile_in = ROOT.TFile.Open("TestFile.root")
dbTree = rfile_in.Get("DBTree")
# get the branch variables
rpcefficiency = RpcEffStruct()
dbTree.SetBranchAddress("rpcefficiency", AddressOf(rpcefficiency, "run_number") )

# loop over DB entries
for i in xrange(dbTree.GetEntries()):
    dbTree.GetEntry(i)
    temprunnum = int(rpcefficiency.run_number)
    if temprunnum in rpcRuns :
        temprunnumindex = rpcRuns.index( temprunnum )
        if temprunnumindex in range(runNumIndexFirst, runNumIndexLast+1) :
            rollid = str(int(rpcefficiency.raw_id))
            if rollid in barrelids :
                rolleff = rpcefficiency.eff_seg
                rollefferr = rpcefficiency.eff_seg_error
                if rollefferr > 0 :
                    (barrelRollsEffDict[ rollid ])[temprunnumindex] = rolleff
                    (barrelRollsErrDict[ rollid ])[temprunnumindex] = rollefferr
                    if rollid in rb3ids :
#                rolleff = rpcefficiency.eff_seg
#                rollefferr = rpcefficiency.eff_seg_error
#                if rollefferr > 0 :
#                    print "GOOD RB3 ROLL FOUND !!!", rollid
#                    print "RUN ", temprunnum, " with index ", temprunnumindex
#                    print "EFF = ", rolleff, " ERR =  ", rollefferr
#                    print "------------------------------------------------"
                        rb3EffInRun[temprunnumindex] = rb3EffInRun[temprunnumindex] + rolleff
#                        rb3EffErrInRun[temprunnumindex] = rb3EffErrInRun[temprunnumindex] + rollefferr * rollefferr
                        rb3EffErrInRun[temprunnumindex] = rb3EffErrInRun[temprunnumindex] + rollefferr 
                        rb3GoodRollsInRun[temprunnumindex] = rb3GoodRollsInRun[temprunnumindex] + 1 


for myRunInd in range(runNumIndexFirst, runNumIndexLast+1) :
    rb3EffInRun[myRunInd] = rb3EffInRun[myRunInd] / float( rb3GoodRollsInRun[myRunInd] )
#    rb3EffErrInRun[myRunInd] = math.sqrt( rb3EffErrInRun[myRunInd] ) / float( rb3GoodRollsInRun[myRunInd] )
    rb3EffErrInRun[myRunInd] = rb3EffErrInRun[myRunInd] / float( rb3GoodRollsInRun[myRunInd] )


############ WRITE REFERENCE EFFICIENCIES INTO ASCII FILE #################
#f = open('refeff_rb3.txt', 'w')
#for myRunInd in range(0, ntotruns) :
#    f.write( str(rpcRuns[myRunInd]) +"\t"+ str(rb3EffInRun[myRunInd])+"\t"+str(rb3EffErrInRun[myRunInd])+"\n" )
#f.close()
############################################################################

#for myids in barrelids :
#    print "EFF. for roll ID ", myids
#    print barrelRollsEffDict[ myids ]
#    print "ERR. for roll ID ", myids
#    print barrelRollsErrDict[ myids ]
#    print "-----------------------------------------"


# MAKE NORMALIZED HISTORY EFF PLOTS FOR EACH BARREL ROLL #######

outfilename = "effhistoryroll_barrel_tree.root"

rfile_out = ROOT.TFile.Open(outfilename,"RECREATE")

histodir = rfile_out.mkdir("rolleffhistory")
histodir.cd()

for myids in barrelids :
    h_name = "eff_"+myids
    h_effruns = ROOT.TH1F(h_name,h_name,ntotruns, 0.5, ntotruns+0.5) 

    rollEffInRun = barrelRollsEffDict[myids]
    rollEffErrInRun = barrelRollsErrDict[myids]
    for m in range(0, ntotruns) :
        eff_norm = 0
        erreff_norm = 0
        if (rb3EffInRun[m] * rollEffInRun[m] * rollEffErrInRun[m]) != 0 :
            eff_norm = rollEffInRun[m] / rb3EffInRun[m]
            relerr1 = rollEffErrInRun[m] / rollEffInRun[m]
            relerr2 = rb3EffErrInRun[m] / rb3EffInRun[m]
            erreff_norm = relerr1 + relerr2
            erreff_norm = erreff_norm * eff_norm 
        h_effruns.SetBinContent( m+1, eff_norm )
        h_effruns.SetBinError( m+1, erreff_norm )
        if m%50 == 0 :
            h_effruns.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))
    h_effruns.Write()

rfile_out.cd()


rpcfit = RpcRollEffFit()

mytree = ROOT.TTree('T', 'T')
mytree.Branch('rpcfit',rpcfit,'rollid/I:p0/F:p1/F:chi2/F:ndof/F') 


for myids in barrelids :
    h_name = "rolleffhistory/eff_"+myids
    myh = ROOT.TH1F()
    rfile_out.GetObject(h_name,myh)
    mchi2 = 0
    mndf = 0
    mp0 = 0
    mp1 = 0

    if myh.Integral() > 0 :
        fitstat = myh.Fit("pol1","S","",runNumIndexFirst+1,runNumIndexLast+1)
        mchi2 = fitstat.Chi2()
        mndf = fitstat.Ndf()
        mp0 = fitstat.Parameter(0)
        mp1 = fitstat.Parameter(1)
        
        print mchi2
        print mndf
        print mp0
        print mp1
        
    rpcfit.rollid = int(myids)
    rpcfit.p0 = mp0
    rpcfit.p1 = mp1
    rpcfit.chi2 = mchi2
    rpcfit.ndof = mndf
    mytree.Fill()

mytree.Print()
mytree.Write()

rfile_out.Close()   

print runNumIndexFirst, runNumIndexLast








