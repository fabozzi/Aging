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

#for i in range(0, ntotruns):
#    print rpcRuns[i], rpcRunsTime[i]

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

# define vectors to store reference (RE1Ring3) CS and errors by run
re1Ring3CSInRun = []
re1Ring3CSErrInRun = []
re1Ring3GoodRollsInRun = []

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
    re1Ring3CSInRun.append(0.0)
    re1Ring3CSErrInRun.append(0.0)
    re1Ring3GoodRollsInRun.append(0.0)
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


# fill the CS dictionaries according to Barrel or Endcap roll
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
                        if rollid in re1Ring3ids :
                            re1Ring3CSInRun[temprunnumindex] = re1Ring3CSInRun[temprunnumindex] + rollCS
                            re1Ring3CSErrInRun[temprunnumindex] = re1Ring3CSErrInRun[temprunnumindex] + sigmaCS 
                            re1Ring3GoodRollsInRun[temprunnumindex] = re1Ring3GoodRollsInRun[temprunnumindex] + 1 
                        


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
    if (re1Ring3GoodRollsInRun[myRunInd] * re1Ring3CSInRun[myRunInd] * re1Ring3CSErrInRun[myRunInd] ) != 0 :
        re1Ring3CSInRun[myRunInd] = re1Ring3CSInRun[myRunInd] / float( re1Ring3GoodRollsInRun[myRunInd] )
        re1Ring3CSErrInRun[myRunInd] = re1Ring3CSErrInRun[myRunInd] / float( re1Ring3GoodRollsInRun[myRunInd] )


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

outfilename = "CS_history_Graph_new.root"

rfile_out = ROOT.TFile.Open(outfilename,"RECREATE")

print "WRITE HISTORY PLOTS NORMALIZED TO RB3"

vx = ROOT.TVectorF(ntotruns)
vy = ROOT.TVectorF(ntotruns)
vxerr = ROOT.TVectorF(ntotruns)
vyerr = ROOT.TVectorF(ntotruns)

for i in range(0, ntotruns):
    vx[i] = rpcRunsTime[i]
    vxerr[i] = 0.0


histodir1 = rfile_out.mkdir("CShistory_Norm")
histodir1.cd()

for myids in barrelids :
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')

    h1_name = "CS_"+tmpstring+"_Norm"
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
            vy[m] = CS_norm
            vyerr[m] = errCS_norm

    h1_CSruns = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h1_CSruns.SetNameTitle(h1_name,h1_name)

    h1_CSruns.Write()


for myids in endcapids :
    tmpstring = (endcapRollsDict[myids]).replace('+','p').replace('-','m')

    h3_name = "CS_"+tmpstring+"_Norm"
    rollCSInRun = endcapRollsCSDict[myids]
    rollCSErrInRun = endcapRollsCSErrDict[myids]

    for m in range(0, ntotruns) :
        CS_norm = 0
        errCS_norm = 0
        if (re1Ring3CSInRun[m] * rollCSInRun[m] * rollCSErrInRun[m]) != 0 :
            CS_norm = rollCSInRun[m] / re1Ring3CSInRun[m]
            relerr1 = rollCSErrInRun[m] / rollCSInRun[m]
            relerr2 = re1Ring3CSErrInRun[m] / re1Ring3CSInRun[m]
            errCS_norm = relerr1 + relerr2
            errCS_norm = errCS_norm * CS_norm 
            vy[m] = CS_norm
            vyerr[m] = errCS_norm

    h3_CSruns = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h3_CSruns.SetNameTitle(h3_name,h3_name)

    h3_CSruns.Write()

rfile_out.cd()


print "WRITE HISTORY PLOTS"

histodir2 = rfile_out.mkdir("CShistory")
histodir2.cd()

for myids in barrelids :
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')

    h2_name = "CS_"+tmpstring
    rollCSInRun = barrelRollsCSDict[myids]
    rollCSErrInRun = barrelRollsCSErrDict[myids]

    for m in range(0, ntotruns):
        vy[m] = rollCSInRun[m]
        vyerr[m] = rollCSErrInRun[m]

    h2_CS = ROOT.TGraphErrors(vx, vy, vxerr, vyerr) 
    h2_CS.SetNameTitle(h2_name,h2_name)

    h2_CS.Write()

for myids in endcapids :
    tmpstring = (endcapRollsDict[myids]).replace('+','p').replace('-','m')

    h4_name = "CS_"+tmpstring

    rollCSInRun = endcapRollsCSDict[myids]
    rollCSErrInRun = endcapRollsCSErrDict[myids]
    for m in range(0, ntotruns) :
        vy[m] = rollCSInRun[m]
        vyerr[m] = rollCSErrInRun[m]

    h4_CS = ROOT.TGraphErrors(vx, vy, vxerr, vyerr) 
    h4_CS.SetNameTitle(h4_name,h4_name)

    h4_CS.Write()


rfile_out.cd()

print "WRITE HISTORY PLOTS NORMALIZED TO RB3 in the wheel"

histodir5 = rfile_out.mkdir("CShistory_RB3Norm_InWh")
histodir5.cd()

for myids in barrelids :
    wh = ((barrelRollsDict[myids]).split('_'))[0]
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')

    h5_name = "CS_"+tmpstring+"_RB3Norm_InWh"
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
            vy[m] = CS_norm
            vyerr[m] = errCS_norm

    h5_CSruns = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h5_CSruns.SetNameTitle(h5_name,h5_name)

    h5_CSruns.Write()

rfile_out.cd()


rfile_out.Close()   

print "DONE"
print runNumIndexFirst, runNumIndexLast







