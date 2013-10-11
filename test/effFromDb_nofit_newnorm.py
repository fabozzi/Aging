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


print "------------------------------------------------"

# fill the list of rb3 rolls used as reference
rb3refRolls = []
rb3refRollsFile = open("rb3refrolls.txt", "r")
for roll in rb3refRollsFile:
    rb3refRolls.append(roll.rstrip())
#print rb3refRolls


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
re1Ring3ids = re1Ring3RollsDict.keys()
barrelids = barrelRollsDict.keys()
endcapids = endcapRollsDict.keys()

#barrelnames = barrelRollsDict.values()
#print 'barrel ids names:  ', barrelRollsDict.values()

# NOW WE HAVE THE LISTS OF BARREL / ENDCAP / RB3_ROLLS (reference)

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

# define dictionaries to store (for each roll) list of efficiencies and errors by run

barrelRollsEffDict = {}
barrelRollsEffErrDict = {}
endcapRollsEffDict = {}
endcapRollsEffErrDict = {}

for y in barrelids :
    barrelRollsEffDict[ y ] = []
    barrelRollsEffErrDict[ y ] = []
    for j in range(0, ntotruns):
        barrelRollsEffDict[ y ].append(0.0)
        barrelRollsEffErrDict[ y ].append(0.0)


for k in endcapids :
    endcapRollsEffDict[ k ] = []
    endcapRollsEffErrDict[ k ] = []
    for j in range(0, ntotruns):
        endcapRollsEffDict[ k ].append(0.0)
        endcapRollsEffErrDict[ k ].append(0.0)


############################################################

# define vectors to store reference (RB3) efficiencies and errors by run

rb3EffInRun = []
rb3EffErrInRun = []
rb3GoodRollsInRun = []

# define vectors to store reference (RE1Ring3) Eff and errors by run
re1Ring3EffInRun = []
re1Ring3EffErrInRun = []
re1Ring3GoodRollsInRun = []

# FOR EACH WHEEL, define a vector to store reference Eff values by run
rb3EffInRun_ByW = {}
rb3EffErrInRun_ByW = {}
rb3GoodRollsInRun_ByW = {}
for w in wheels_names :
    rb3EffInRun_ByW[w] = []
    rb3EffErrInRun_ByW[w] = []
    rb3GoodRollsInRun_ByW[w] = []

for y in range(0, ntotruns):
    rb3EffInRun.append(0.0)
    rb3EffErrInRun.append(0.0)
    rb3GoodRollsInRun.append(0)
    re1Ring3EffInRun.append(0.0)
    re1Ring3EffErrInRun.append(0.0)
    re1Ring3GoodRollsInRun.append(0.0)
    for w in wheels_names :
        (rb3EffInRun_ByW[w]).append(0.0)
        (rb3EffErrInRun_ByW[w]).append(0.0)
        (rb3GoodRollsInRun_ByW[w]).append(0)
        
                                                        


firstRun = 190679
#firstRun = 206513
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
# get efficiency
                rolleff = (rpcefficiency.eff_seg) / 100.0
# get efficiency error
                rollefferr = (rpcefficiency.eff_seg_error) / 100.0
# get square root of number of entries
                sqNentries = 0.0
                if rolleff * rollefferr > 0 :
                    sqNentries = (math.sqrt( rolleff * (1.0-rolleff) )) / rollefferr
                    if isBarrelRoll :
                        (barrelRollsEffDict[ rollid ])[temprunnumindex] = rolleff
                        (barrelRollsEffErrDict[ rollid ])[temprunnumindex] = rollefferr
                        if rollid in rb3refRolls :
                            rb3EffInRun[temprunnumindex] = rb3EffInRun[temprunnumindex] + rolleff
                            rb3EffErrInRun[temprunnumindex] = rb3EffErrInRun[temprunnumindex] + rollefferr 
                            rb3GoodRollsInRun[temprunnumindex] = rb3GoodRollsInRun[temprunnumindex] + 1 
                            wh = ((rb3RollsDict[rollid]).split('_'))[0]
                            (rb3EffInRun_ByW[wh])[temprunnumindex] = (rb3EffInRun_ByW[wh])[temprunnumindex] + rolleff
                            (rb3EffErrInRun_ByW[wh])[temprunnumindex] = (rb3EffErrInRun_ByW[wh])[temprunnumindex] + rollefferr
                            (rb3GoodRollsInRun_ByW[wh])[temprunnumindex] = (rb3GoodRollsInRun_ByW[wh])[temprunnumindex] + 1
                    else:
                        (endcapRollsEffDict[ rollid ])[temprunnumindex] = rolleff
                        (endcapRollsEffErrDict[ rollid ])[temprunnumindex] = rollefferr
                        if rollid in re1Ring3ids :
                            re1Ring3EffInRun[temprunnumindex] = re1Ring3EffInRun[temprunnumindex] + rolleff
                            re1Ring3EffErrInRun[temprunnumindex] = re1Ring3EffErrInRun[temprunnumindex] + rollefferr
                            re1Ring3GoodRollsInRun[temprunnumindex] = re1Ring3GoodRollsInRun[temprunnumindex] + 1



print "COMPUTING REFRENCE VALUES"

for myRunInd in range(runNumIndexFirst, runNumIndexLast+1) :
    if float( rb3GoodRollsInRun[myRunInd] * rb3EffInRun[myRunInd] * rb3EffErrInRun[myRunInd]) == 0: continue
    rb3EffInRun[myRunInd] = rb3EffInRun[myRunInd] / float( rb3GoodRollsInRun[myRunInd] )
#    rb3EffErrInRun[myRunInd] = math.sqrt( rb3EffErrInRun[myRunInd] ) / float( rb3GoodRollsInRun[myRunInd] )
    rb3EffErrInRun[myRunInd] = rb3EffErrInRun[myRunInd] / float( rb3GoodRollsInRun[myRunInd] )
    for wh in wheels_names :
        if ( (rb3GoodRollsInRun_ByW[wh])[myRunInd] * (rb3EffInRun_ByW[wh])[myRunInd] * (rb3EffErrInRun_ByW[wh])[myRunInd]) == 0 : continue
        (rb3EffInRun_ByW[wh])[myRunInd] = (rb3EffInRun_ByW[wh])[myRunInd] / float( (rb3GoodRollsInRun_ByW[wh])[myRunInd] )
        (rb3EffErrInRun_ByW[wh])[myRunInd] = (rb3EffErrInRun_ByW[wh])[myRunInd] / float( (rb3GoodRollsInRun_ByW[wh])[myRunInd] )
    if (re1Ring3GoodRollsInRun[myRunInd] * re1Ring3EffInRun[myRunInd] * re1Ring3EffErrInRun[myRunInd] ) != 0 :
        re1Ring3EffInRun[myRunInd] = re1Ring3EffInRun[myRunInd] / float( re1Ring3GoodRollsInRun[myRunInd] )
        re1Ring3EffErrInRun[myRunInd] = re1Ring3EffErrInRun[myRunInd] / float( re1Ring3GoodRollsInRun[myRunInd] )
        
                                                        

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

print "DUMPING HISTORY PLOTS"

outfilename = "eff_history_Graph_nofit_newnorm.root"

rfile_out = ROOT.TFile.Open(outfilename,"RECREATE")

print "WRITE HISTORY PLOTS NORMALIZED"


vx = ROOT.TVectorF(ntotruns)
vy = ROOT.TVectorF(ntotruns)
vxerr = ROOT.TVectorF(ntotruns)
vyerr = ROOT.TVectorF(ntotruns)
vgr = ROOT.TVectorF(ntotruns)
vgrerr = ROOT.TVectorF(ntotruns)

for i in range(0, ntotruns):
    vx[i] = rpcRunsTime[i]
    vxerr[i] = 0.0
    vgrerr[i] = 0.0
            

histodir1 = rfile_out.mkdir("effhistory_norm")
histodir1.cd()

for myids in barrelids :
    tmpstring = barrelRollsDict[myids]
    #    print 'string bef:  ', tmpstring
    tmpstring = tmpstring.replace('+','p').replace('-','m')
    #    print 'string aft:  ', tmpstring
    
    h1_name = "eff_"+tmpstring+"_norm"
    #    h1_effruns = ROOT.TH1F(h1_name,h1_name,ntotruns, 0.5, ntotruns+0.5) 

    rollEffInRun = barrelRollsEffDict[myids]
    rollEffErrInRun = barrelRollsEffErrDict[myids]

    for m in range(0, ntotruns) :
        eff_norm = 0
        erreff_norm = 0
        if (rb3EffInRun[m] * rollEffInRun[m] * rollEffErrInRun[m]) != 0 :
            eff_norm = rollEffInRun[m] / rb3EffInRun[m]
            relerr1 = rollEffErrInRun[m] / rollEffInRun[m]
            relerr2 = rb3EffErrInRun[m] / rb3EffInRun[m]
            erreff_norm = relerr1 + relerr2
            erreff_norm = erreff_norm * eff_norm 
            vy[m] = eff_norm
            vyerr[m] = erreff_norm
    h1_effruns = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h1_effruns.SetNameTitle(h1_name,h1_name)
#    h_effruns.SetBinContent( m+1, eff_norm )
#    h_effruns.SetBinError( m+1, erreff_norm )
#        if m%50 == 0 :
#            h_effruns.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))

    h1_effruns.Write()

for myids in endcapids :
    tmpstring = (endcapRollsDict[myids]).replace('+','p').replace('-','m')
    
    h3_name = "eff_"+tmpstring+"_norm"
    rollEffInRun = endcapRollsEffDict[myids]
    rollEffErrInRun = endcapRollsEffErrDict[myids]
    
    for m in range(0, ntotruns) :
        eff_norm = 0
        erreff_norm = 0
        if (re1Ring3EffInRun[m] * rollEffInRun[m] * rollEffErrInRun[m]) != 0 :
            eff_norm = rollEffInRun[m] / re1Ring3EffInRun[m]
            relerr1 = rollEffErrInRun[m] / rollEffInRun[m]
            relerr2 = re1Ring3EffErrInRun[m] / re1Ring3EffInRun[m]
            erreff_norm = relerr1 + relerr2
            erreff_norm = erreff_norm * eff_norm
            vy[m] = eff_norm
            vyerr[m] = erreff_norm
            
    h3_effruns = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h3_effruns.SetNameTitle(h3_name,h3_name)
            
    h3_effruns.Write()
    
print "WRITE HISTORY PLOTS FOR REFERENCES"

for m in range(0, ntotruns) :
    if rb3EffInRun[m] != 0 :
        vy[m] = rb3EffInRun[m]
        vyerr[m] = rb3EffErrInRun[m]
        vgr[m] = rb3GoodRollsInRun[m]
hreference_effruns = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
hreference_ngood = ROOT.TGraphErrors(vx, vgr, vxerr, vgrerr)
hreference_effruns.SetNameTitle("rb3ref_eff","rb3ref_eff")
hreference_ngood.SetNameTitle("rb3ref_ngood","rb3ref_ngood")
hreference_effruns.Write()
hreference_ngood.Write()

for m in range(0, ntotruns) :
    if re1Ring3EffInRun[m] != 0 :
        vy[m] = re1Ring3EffInRun[m]
        vyerr[m] = re1Ring3EffErrInRun[m]
        vgr[m] = re1Ring3GoodRollsInRun[m]
hreference1_effruns = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
hreference1_ngood = ROOT.TGraphErrors(vx, vgr, vxerr, vgrerr)
hreference1_effruns.SetNameTitle("re1ring3ref_eff","re1ring3ref_eff")
hreference1_ngood.SetNameTitle("re1ring3ref_ngood","re1ring3ref_ngood")
hreference1_effruns.Write()
hreference1_ngood.Write()

rfile_out.cd()


print "WRITE HISTORY PLOTS"

histodir2 = rfile_out.mkdir("effhistory")
histodir2.cd()

for myids in barrelids :
    tmpstring = barrelRollsDict[myids]
    #    print 'string bef:  ', tmpstring
    tmpstring = tmpstring.replace('+','p').replace('-','m')
    #    print 'string aft:  ', tmpstring
    
    
    h2_name = "eff_"+tmpstring
    #    h2_effruns = ROOT.TH1F(h_name,h_name,ntotruns, 0.5, ntotruns+0.5)
    
    rollEffInRun = barrelRollsEffDict[myids]
    rollEffErrInRun = barrelRollsEffErrDict[myids]
    for m in range(0, ntotruns) :
        #        h_effruns.SetBinContent( m+1, rollEffInRun[m] )
        #        h_effruns.SetBinError( m+1, rollEffErrInRun[m] )
        #        if m%50 == 0 :
        #            h_effruns.GetXaxis().SetBinLabel(m+1,str(rpcRuns[m]))
        vy[m] = rollEffInRun[m]
        vyerr[m] = rollEffErrInRun[m]
                        
    h2_eff = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h2_eff.SetNameTitle(h2_name,h2_name)
    h2_eff.Write()
                

for myids in endcapids :
    tmpstring = (endcapRollsDict[myids]).replace('+','p').replace('-','m')
    
    h4_name = "eff_"+tmpstring
    
    rollEffInRun = endcapRollsEffDict[myids]
    rollEffErrInRun = endcapRollsEffErrDict[myids]
    for m in range(0, ntotruns) :
        vy[m] = rollEffInRun[m]
        vyerr[m] = rollEffErrInRun[m]
        
    h4_eff = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h4_eff.SetNameTitle(h4_name,h4_name)
    h4_eff.Write()

                                                    
rfile_out.cd()

print "WRITE HISTORY PLOTS NORMALIZED TO RB3 in the wheel"

histodir5 = rfile_out.mkdir("effhistory_RB3norm_InWh")
histodir5.cd()

for myids in barrelids :
    wh = ((barrelRollsDict[myids]).split('_'))[0]
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')
    
    h5_name = "eff_"+tmpstring+"_RB3norm_InWh"
    rollEffInRun = barrelRollsEffDict[myids]
    rollEffErrInRun = barrelRollsEffErrDict[myids]

    for m in range(0, ntotruns) :
        eff_norm = 0
        erreff_norm = 0
        if ((rb3EffInRun_ByW[wh])[m] * rollEffInRun[m] * rollEffErrInRun[m]) != 0 :
            eff_norm = rollEffInRun[m] / (rb3EffInRun_ByW[wh])[m]
            relerr1 = rollEffErrInRun[m] / rollEffInRun[m]
            relerr2 = (rb3EffErrInRun_ByW[wh])[m] / (rb3EffInRun_ByW[wh])[m]
            erreff_norm = relerr1 + relerr2
            erreff_norm = erreff_norm * eff_norm
            vy[m] = eff_norm
            vyerr[m] = erreff_norm
            
    h5_effruns = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h5_effruns.SetNameTitle(h5_name,h5_name)
            
    h5_effruns.Write()
    

rfile_out.cd()
rfile_out.Close()


print "DONE"

print runNumIndexFirst, runNumIndexLast








