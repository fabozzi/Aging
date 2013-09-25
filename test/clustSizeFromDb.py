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

#gROOT.ProcessLine(
#    "struct RpcRollEffFit {\
#    Int_t rollid;\
#    Float_t p0;\
#    Float_t p1;\
#    Float_t chi2;\
#    Float_t ndof;\
#} ;");


ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()        # don't pop up canvases


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
#    else :
#        print " BAD ROLL FOUND !!!!!!!!!!!!!!"

rb3ids = rb3RollsDict.keys()
barrelids = barrelRollsDict.keys()
endcapids = endcapRollsDict.keys()

# NOW WE HAVE THE LISTS OF BARREL / ENDCAP / RB3_ROLLS (reference)

# GET THE RUN LIST
# open the official input run list from a file and make a vector from it
rpcRuns = []
runListFile = open("runslist.txt", "r")
for eachRun in runListFile:
    rpcRuns.append(int(eachRun.rstrip()))
ntotruns = len(rpcRuns) 
############################################################

# define dictionaries to store (for each roll) list of CS and CS errors by run

barrelRollsCSDict = {}
barrelRollsCSErrDict = {}
endcapRollsCSDict = {}
endcapRollsCSErrDict = {}

for y in barrelids :
    barrelRollsCSDict[ y ] = []
    barrelRollsCSErrDict[ y ] = []
    for j in range(0, ntotruns):
        barrelRollsCSDict[ y ].append(0.0)
        barrelRollsCSErrDict[ y ].append(0.0)

for k in endcapids :
    endcapRollsCSDict[ k ] = []
    endcapRollsCSErrDict[ k ] = []
    for j in range(0, ntotruns):
        endcapRollsCSDict[ k ].append(0.0)
        endcapRollsCSErrDict[ k ].append(0.0)


############################################################

# define vectors to store reference (RB3) CS and errors by run

rb3CSInRun = []
rb3CSErrInRun = []
rb3GoodRollsInRun = []

# FOR EACH WHEEL, define a vector to store reference CS values by run
rb3CSInRun_ByW = {}
rb3CSErrInRun_ByW = {}
rb3GoodRollsInRun_ByW = {}
for w in wheels_names :
    rb3CSInRun_ByW[w] = []
    rb3CSErrInRun_ByW[w] = []
    rb3GoodRollsInRun_ByW[w] = []


for y in range(0, ntotruns):
    rb3CSInRun.append(0.0)
    rb3CSErrInRun.append(0.0)
    rb3GoodRollsInRun.append(0)
    for w in wheels_names :
        (rb3CSInRun_ByW[w]).append(0.0)
        (rb3CSErrInRun_ByW[w]).append(0.0)
        (rb3GoodRollsInRun_ByW[w]).append(0)


firstRun = 190679
lastRun = 209151
runNumIndexFirst = rpcRuns.index( firstRun )
runNumIndexLast = rpcRuns.index( lastRun )


print "LOOPING OVER TREE"

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
            isBarrelRoll = (rollid in barrelids) 
            isEndcapRoll = (rollid in endcapids)
            if (isBarrelRoll or isEndcapRoll) :
# get average cluster size
                rollCS = rpcefficiency.cluster_size
# get efficiency
                rolleff = (rpcefficiency.eff_seg) / 100.0
# get efficiency error
                rollefferr = (rpcefficiency.eff_seg_error) / 100.0
# get square root of number of entries
                sqNentries = 0.0
                if rolleff * rollefferr > 0 :
                    sqNentries = (math.sqrt( rolleff * (1.0-rolleff) )) / rollefferr

# compute cluster size probabilities
                cs_prob = []
                cs_prob.append(rpcefficiency.clus_size_bin01 - rpcefficiency.clus_size_bin02)
                cs_prob.append(rpcefficiency.clus_size_bin02 - rpcefficiency.clus_size_bin03)
                cs_prob.append(rpcefficiency.clus_size_bin03 - rpcefficiency.clus_size_bin04)
                cs_prob.append(rpcefficiency.clus_size_bin04 - rpcefficiency.clus_size_bin05)
                cs_prob.append(rpcefficiency.clus_size_bin05 - rpcefficiency.clus_size_bin06)
                cs_prob.append(rpcefficiency.clus_size_bin06 - rpcefficiency.clus_size_bin07)
                cs_prob.append(rpcefficiency.clus_size_bin07 - rpcefficiency.clus_size_bin08)
                cs_prob.append(rpcefficiency.clus_size_bin08 - rpcefficiency.clus_size_bin09)
                cs_prob.append(rpcefficiency.clus_size_bin09 - rpcefficiency.clus_size_bin10)
                cs_prob.append(rpcefficiency.clus_size_bin10)

#               rollCSerr = rpcefficiency.cluster_size_error
#                rollCSerr = 1.0
                if rollCS * sqNentries > 0 :
#                    print "--------------------------------"
#                    print "MEAN CS = ", rollCS
#                    averollCS = 0

# compute error on the average CS
                    rms2CS = 0
                    sigmaSqCS = 0
                    sigmaCS = 0
                    for nstrips in range (1,11) :
#                        averollCS = averollCS + nstrips * cs_prob[nstrips-1]
                        rms2CS = rms2CS + nstrips * nstrips * cs_prob[nstrips-1]
#                        print rms2CS
                    sigmaSqCS = rms2CS - rollCS * rollCS
                    if sigmaSqCS > 0 :
                        sigmaCS = (math.sqrt(sigmaSqCS)) / sqNentries


                    if rollCS < 1 :
                        print "--------------------------------"
                        print "MEAN CS = ", rollCS
                        print cs_prob
                        print "rms CS = ", rms2CS
                        print "error CS = ", sigmaCS


# fill the CS disctionaries according to Barrel or Endcap roll
                    if isBarrelRoll :
                        (barrelRollsCSDict[ rollid ])[temprunnumindex] = rollCS
                        (barrelRollsCSErrDict[ rollid ])[temprunnumindex] = sigmaCS
                        if rollid in rb3ids :
                            rb3CSInRun[temprunnumindex] = rb3CSInRun[temprunnumindex] + rollCS
#                        rb3CSErrInRun[temprunnumindex] = rb3CSErrInRun[temprunnumindex] + rollCSerr * rollCSerr
                            rb3CSErrInRun[temprunnumindex] = rb3CSErrInRun[temprunnumindex] + sigmaCS 
                            rb3GoodRollsInRun[temprunnumindex] = rb3GoodRollsInRun[temprunnumindex] + 1 
                            wh = ((rb3RollsDict[rollid]).split('_'))[0]
                            (rb3CSInRun_ByW[wh])[temprunnumindex] = (rb3CSInRun_ByW[wh])[temprunnumindex] + rollCS
                            (rb3CSErrInRun_ByW[wh])[temprunnumindex] = (rb3CSErrInRun_ByW[wh])[temprunnumindex] + sigmaCS 
                            (rb3GoodRollsInRun_ByW[wh])[temprunnumindex] = (rb3GoodRollsInRun_ByW[wh])[temprunnumindex] + 1 
                    else:
                        (endcapRollsCSDict[ rollid ])[temprunnumindex] = rollCS
                        (endcapRollsCSErrDict[ rollid ])[temprunnumindex] = sigmaCS
                        

print "COMPUTING REFRENCE VALUES"

# compute average CS and error for RB3 rolls (reference value) 
for myRunInd in range(runNumIndexFirst, runNumIndexLast+1) :
    if (rb3GoodRollsInRun[myRunInd] * rb3CSInRun[myRunInd] * rb3CSErrInRun[myRunInd] ) != 0 :
        rb3CSInRun[myRunInd] = rb3CSInRun[myRunInd] / float( rb3GoodRollsInRun[myRunInd] )
#    rb3CSErrInRun[myRunInd] = math.sqrt( rb3CSErrInRun[myRunInd] ) / float( rb3GoodRollsInRun[myRunInd] )
        rb3CSErrInRun[myRunInd] = rb3CSErrInRun[myRunInd] / float( rb3GoodRollsInRun[myRunInd] )
    for wh in wheels_names :
        if ( (rb3GoodRollsInRun_ByW[wh])[myRunInd] * (rb3CSInRun_ByW[wh])[myRunInd] * (rb3CSErrInRun_ByW[wh])[myRunInd] ) != 0 :
            (rb3CSInRun_ByW[wh])[myRunInd] = (rb3CSInRun_ByW[wh])[myRunInd] / float( (rb3GoodRollsInRun_ByW[wh])[myRunInd] )
            (rb3CSErrInRun_ByW[wh])[myRunInd] = (rb3CSErrInRun_ByW[wh])[myRunInd] / float( (rb3GoodRollsInRun_ByW[wh])[myRunInd] )
        

############ WRITE REFERENCE CS INTO ASCII FILE #################
#f = open('refCS_rb3.txt', 'w')
#for myRunInd in range(0, ntotruns) :
#    f.write( str(rpcRuns[myRunInd]) +"\t"+ str(rb3CSInRun[myRunInd])+"\t"+str(rb3CSErrInRun[myRunInd])+"\n" )
#f.close()
############################################################################

#for myids in barrelids :
#    print "CS. for roll ID ", myids
#    print barrelRollsCSDict[ myids ]
#    print "ERR. for roll ID ", myids
#    print barrelRollsCSErrDict[ myids ]
#    print "-----------------------------------------"


# MAKE NORMALIZED HISTORY CS PLOTS FOR EACH ROLL #######

print "DUMPING HISTORY PLOTS"

outfilename = "CS_history_new.root"

rfile_out = ROOT.TFile.Open(outfilename,"RECREATE")

#histodir3 = rfile_out.mkdir("CShistory_endcap_normAllRB3")
#histodir4 = rfile_out.mkdir("CShistory_endcap")


print "WRITE HISTORY PLOTS NORMALIZED TO RB3"

histodir1 = rfile_out.mkdir("CShistory_normAllRB3")
histodir1.cd()

for myids in barrelids :
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')

    h1_name = "CS_"+tmpstring+"_RB3Norm"
    h1_CSruns = ROOT.TH1F(h1_name,h1_name,ntotruns, 0.5, ntotruns+0.5) 

    rollCSInRun = barrelRollsCSDict[myids]
    rollCSErrInRun = barrelRollsCSErrDict[myids]
    for m in range(0, ntotruns) :
        CS_norm = 0
        errCS_norm = 0
        if (rb3CSInRun[m] * rollCSInRun[m] * rollCSErrInRun[m]) != 0 :
            CS_norm = rollCSInRun[m] / rb3CSInRun[m]
            relerr1 = rollCSErrInRun[m] / rollCSInRun[m]
            relerr2 = rb3CSErrInRun[m] / rb3CSInRun[m]
            errCS_norm = relerr1 + relerr2
            errCS_norm = errCS_norm * CS_norm 
        h1_CSruns.SetBinContent( m+1, CS_norm )
        h1_CSruns.SetBinError( m+1, errCS_norm )
        if m%50 == 0 :
            h1_CSruns.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))
    h1_CSruns.Write()

for myids in endcapids :
    tmpstring = (endcapRollsDict[myids]).replace('+','p').replace('-','m')

    h3_name = "CS_"+tmpstring+"_RB3Norm"
    h3_CSruns = ROOT.TH1F(h3_name,h3_name,ntotruns, 0.5, ntotruns+0.5) 

    rollCSInRun = endcapRollsCSDict[myids]
    rollCSErrInRun = endcapRollsCSErrDict[myids]
    for m in range(0, ntotruns) :
        CS_norm = 0
        errCS_norm = 0
        if (rb3CSInRun[m] * rollCSInRun[m] * rollCSErrInRun[m]) != 0 :
            CS_norm = rollCSInRun[m] / rb3CSInRun[m]
            relerr1 = rollCSErrInRun[m] / rollCSInRun[m]
            relerr2 = rb3CSErrInRun[m] / rb3CSInRun[m]
            errCS_norm = relerr1 + relerr2
            errCS_norm = errCS_norm * CS_norm 
        h3_CSruns.SetBinContent( m+1, CS_norm )
        h3_CSruns.SetBinError( m+1, errCS_norm )
        if m%50 == 0 :
            h3_CSruns.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))
    h3_CSruns.Write()


rfile_out.cd()

print "WRITE HISTORY PLOTS"

histodir2 = rfile_out.mkdir("CShistory")
histodir2.cd()

for myids in barrelids :
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')

    h2_name = "CS_"+tmpstring
    h2_CS = ROOT.TH1F(h2_name,h2_name,ntotruns, 0.5, ntotruns+0.5) 

    rollCSInRun = barrelRollsCSDict[myids]
    rollCSErrInRun = barrelRollsCSErrDict[myids]
    for m in range(0, ntotruns) :
        h2_CS.SetBinContent( m+1, rollCSInRun[m] )
        h2_CS.SetBinError( m+1, rollCSErrInRun[m] )
        if m%50 == 0 :
            h2_CS.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))
    h2_CS.Write()

for myids in endcapids :
    tmpstring = (endcapRollsDict[myids]).replace('+','p').replace('-','m')

    h4_name = "CS_"+tmpstring
    h4_CS = ROOT.TH1F(h4_name,h4_name,ntotruns, 0.5, ntotruns+0.5) 

    rollCSInRun = endcapRollsCSDict[myids]
    rollCSErrInRun = endcapRollsCSErrDict[myids]
    for m in range(0, ntotruns) :
        h4_CS.SetBinContent( m+1, rollCSInRun[m] )
        h4_CS.SetBinError( m+1, rollCSErrInRun[m] )
        if m%50 == 0 :
            h4_CS.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))
    h4_CS.Write()


rfile_out.cd()


print "WRITE HISTORY PLOTS NORMALIZED TO RB3 in the wheel"

histodir5 = rfile_out.mkdir("CShistory_RB3Norm_InWh")
histodir5.cd()

for myids in barrelids :
    wh = ((barrelRollsDict[myids]).split('_'))[0]
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')

    h5_name = "CS_"+tmpstring+"_RB3Norm_InWh"
    h5_CS = ROOT.TH1F(h5_name,h5_name, ntotruns, 0.5, ntotruns+0.5) 

    rollCSInRun = barrelRollsCSDict[myids]
    rollCSErrInRun = barrelRollsCSErrDict[myids]
    for m in range(0, ntotruns) :
        CS_norm = 0
        errCS_norm = 0
        if ((rb3CSInRun_ByW[wh])[m] * rollCSInRun[m] * rollCSErrInRun[m]) != 0 :
            CS_norm = rollCSInRun[m] / (rb3CSInRun_ByW[wh])[m]
            relerr1 = rollCSErrInRun[m] / rollCSInRun[m]
            relerr2 = (rb3CSErrInRun_ByW[wh])[m] / (rb3CSInRun_ByW[wh])[m]
            errCS_norm = relerr1 + relerr2
            errCS_norm = errCS_norm * CS_norm 
        h5_CS.SetBinContent( m+1, CS_norm )
        h5_CS.SetBinError( m+1, errCS_norm )
        if m%50 == 0 :
            h5_CS.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))
    h5_CS.Write()


rfile_out.cd()

rfile_out.Close()   

print "DONE"
print runNumIndexFirst, runNumIndexLast







