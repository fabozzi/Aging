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
    Float_t rollidtree;\
    Float_t p0;\
    Float_t p0err;\
    Float_t p1;\
    Float_t p1err;\
    Float_t chi2;\
    Float_t ndof;\
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

runListFile = open("runslist_small.txt", "r")
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
                        if rollid in rb3ids :
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

outfilename = "eff_history_Graph.root"

rfile_out = ROOT.TFile.Open(outfilename,"RECREATE")

print "WRITE HISTORY PLOTS NORMALIZED TO RB3"

rollidtree = ROOT.TVectorF(1200)
mchi2 = ROOT.TVectorF(1200)
mndf = ROOT.TVectorF(1200)
mp0 = ROOT.TVectorF(1200)
mp1 = ROOT.TVectorF(1200)
mp0err = ROOT.TVectorF(1200)
mp1err = ROOT.TVectorF(1200)
j = 0

vx = ROOT.TVectorF(ntotruns)
vy = ROOT.TVectorF(ntotruns)
vxerr = ROOT.TVectorF(ntotruns)
vyerr = ROOT.TVectorF(ntotruns)

for i in range(0, ntotruns):
    vx[i] = rpcRunsTime[i]
    vxerr[i] = 0.0
            
rfile_out.cd()
histodir1 = rfile_out.mkdir("Effhistory_normAllRB3")
histodir1.cd()

for myids in barrelids :

    histodir1.cd()
    tmpstring = barrelRollsDict[myids]
    #    print 'string bef:  ', tmpstring
    tmpstring = tmpstring.replace('+','p').replace('-','m')
    #    print 'string aft:  ', tmpstring
    
    
    h1_name = "eff_"+tmpstring+"_RB3Norm"
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

    tf1 = ROOT.TF1("tf1","[0]+x*[1]")
    fitstat1 = h1_effruns.Fit(tf1,"S","",2012.0,2012.95816911)
    rollidtree[j] = float(myids)
    mchi2[j] = fitstat1.Chi2()
    mndf[j] = fitstat1.Ndf()
    mp0[j] = fitstat1.Parameter(0) 
    mp0err[j] = fitstat1.ParError(0)
    mp1[j] = fitstat1.Parameter(1)
    mp1err[j] = fitstat1.ParError(1)
    
    #    p = ROOT.TPaveStats(0.1,0.8,0.3,0.95)
    #    p = h1_effruns.GetListOfFunctions().FindObject("stats")
    #    p.Write()
    h1_effruns.Write()
    j=j+1


#histodirTree1 = rfile_out.mkdir("FitParameters_normAllRB3")
histodir1.cd()

efffit_barrel = RpcRollEffFit()

mytree = ROOT.TTree('FitParameters_normAllRB3', 'FitParameters_normAllRB3')
mytree.Branch('efffit_barrel',efffit_barrel,'rollidtree/I:p0/F:p0err/F:p1/F:p1err/F:chi2/F:ndof/F')

for i in range(0,j):
    efffit_barrel.rollidtree = rollidtree[i]
    efffit_barrel.p0 = mp0[i] + 2012*mp1[i]
    efffit_barrel.p0err = math.sqrt(math.pow(mp0err[i],2) + math.pow(2012,2) * math.pow(mp1err[i],2))
    efffit_barrel.p1 = mp1[i]
    efffit_barrel.p1err = mp1err[i]
    efffit_barrel.chi2 = mchi2[i]
    efffit_barrel.ndof = mndf[i]

    mytree.Fill()
        
mytree.Print()
mytree.Write()

rfile_out.cd()
histodir1.cd()
j = 0
    

for myids in endcapids :
    tmpstring = (endcapRollsDict[myids]).replace('+','p').replace('-','m')
    
    h3_name = "eff_"+tmpstring+"_RE1RING3Norm"
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
    tf3 = ROOT.TF1("tf3","[0]+x*[1]")
    fitstat3 = h3_effruns.Fit(tf3,"S","",2012.0,2012.95816911)
    rollidtree[j] = float(myids)
    mchi2[j] = fitstat3.Chi2()
    mndf[j] = fitstat3.Ndf()
    mp0[j] = fitstat3.Parameter(0)
    mp0err[j] = fitstat3.ParError(0)
    mp1[j] = fitstat3.Parameter(1)
    mp1err[j] = fitstat3.ParError(1)
                            
    h3_effruns.Write()
    j=j+1
       
rfile_out.cd()

histodir1.cd()

efffit_endcap = RpcRollEffFit()

mytree = ROOT.TTree('FitParameters_normAllre1Ring3', 'FitParameters_normAllre1Ring3')
mytree.Branch('efffit_endcap',efffit_endcap,'rollidtree/I:p0/F:p0err/F:p1/F:p1err/F:chi2/F:ndof/F')

for i in range(0,j):
    efffit_endcap.rollidtree = rollidtree[i]
    efffit_endcap.p0 = mp0[i] + 2012*mp1[i]
    efffit_endcap.p0err = math.sqrt(math.pow(mp0err[i],2) + math.pow(2012,2) * math.pow(mp1err[i],2))
    efffit_endcap.p1 = mp1[i]
    efffit_endcap.p1err = mp1err[i]
    efffit_endcap.chi2 = mchi2[i]
    efffit_endcap.ndof = mndf[i]
    
    mytree.Fill()

mytree.Print()
mytree.Write()

                                    
print "WRITE HISTORY PLOTS"

rfile_out.cd()
histodir2 = rfile_out.mkdir("Effhistory")
histodir2.cd()
j = 0


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
    tf2 = ROOT.TF1("tf2","[0]+x*[1]")
    fitstat2 = h2_eff.Fit(tf2,"S","",2012.0,2012.95816911)
    rollidtree[j] = float(myids)
    mchi2[j] = fitstat2.Chi2()
    mndf[j] = fitstat2.Ndf()
    mp0[j] = fitstat2.Parameter(0)
    mp0err[j] = fitstat2.ParError(0)
    mp1[j] = fitstat2.Parameter(1)
    mp1err[j] = fitstat2.ParError(1)
    
    h2_eff.Write()
    j=j+1
    

##rfile_out.cd()
efffit_barrel = RpcRollEffFit()

mytree = ROOT.TTree('FitParameters_barrel', 'FitParameters_barrel')
mytree.Branch('efffit_barrel',efffit_barrel,'rollidtree/I:p0/F:p0err/F:p1/F:p1err/F:chi2/F:ndof/F')

for i in range(0,j):
    efffit_barrel.rollidtree = rollidtree[i]
    efffit_barrel.p0 = mp0[i] + 2012*mp1[i]
    efffit_barrel.p0err = math.sqrt(math.pow(mp0err[i],2) + math.pow(2012,2) * math.pow(mp1err[i],2))
    efffit_barrel.p1 = mp1[i]
    efffit_barrel.p1err = mp1err[i]
    efffit_barrel.chi2 = mchi2[i]
    efffit_barrel.ndof = mndf[i]
    
    mytree.Fill()
    
mytree.Print()
mytree.Write()
    
j = 0
rfile_out.cd()
histodir2.cd()

for myids in endcapids :
    tmpstring = (endcapRollsDict[myids]).replace('+','p').replace('-','m')
    
    h4_name = "Eff_"+tmpstring
    
    rollEffInRun = endcapRollsEffDict[myids]
    rollEffErrInRun = endcapRollsEffErrDict[myids]
    for m in range(0, ntotruns) :
        vy[m] = rollEffInRun[m]
        vyerr[m] = rollEffErrInRun[m]
        
    h4_eff = ROOT.TGraphErrors(vx, vy, vxerr, vyerr)
    h4_eff.SetNameTitle(h4_name,h4_name)
    tf4 = ROOT.TF1("tf4","[0]+x*[1]")
    fitstat4 = h4_eff.Fit(tf4,"S","",2012.0,2012.95816911)
    rollidtree[j] = float(myids)
    mchi2[j] = fitstat4.Chi2()
    mndf[j] = fitstat4.Ndf()
    mp0[j] = fitstat4.Parameter(0)
    mp0err[j] = fitstat4.ParError(0)
    mp1[j] = fitstat4.Parameter(1)
    mp1err[j] = fitstat4.ParError(1)
    
    h4_eff.Write()
    j=j+1


histodir2.cd()
efffit_endcap = RpcRollEffFit()

mytree = ROOT.TTree('FitParameters_endcap', 'FitParameters_encdap')
mytree.Branch('efffit_endcap',efffit_endcap,'rollidtree/I:p0/F:p0err/F:p1/F:p1err/F:chi2/F:ndof/F')

for i in range(0,j):
    efffit_endcap.rollidtree = rollidtree[i]
    efffit_endcap.p0 = mp0[i] + 2012*mp1[i]
    efffit_endcap.p0err = math.sqrt(math.pow(mp0err[i],2) + math.pow(2012,2) * math.pow(mp1err[i],2))
    efffit_endcap.p1 = mp1[i]
    efffit_endcap.p1err = mp1err[i]
    efffit_endcap.chi2 = mchi2[i]
    efffit_endcap.ndof = mndf[i]
    
    mytree.Fill()
    
mytree.Print()
mytree.Write()
                                    
print "WRITE HISTORY PLOTS NORMALIZED TO RB3 in the wheel"


histodir5 = rfile_out.mkdir("Effhistory_RB3Norm_InWh")
rfile_out.cd()
histodir5.cd()
j = 0


for myids in barrelids :
    wh = ((barrelRollsDict[myids]).split('_'))[0]
    tmpstring = (barrelRollsDict[myids]).replace('+','p').replace('-','m')
    
    h5_name = "Eff_"+tmpstring+"_RB3Norm_InWh"
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
    tf5 = ROOT.TF1("tf5","[0]+x*[1]")
    fitstat5 = h5_effruns.Fit(tf5,"S","",2012.0,2012.95816911)
    rollidtree[j] = float(myids)
    mchi2[j] = fitstat5.Chi2()
    mndf[j] = fitstat5.Ndf()
    mp0[j] = fitstat5.Parameter(0)
    mp0err[j] = fitstat5.ParError(0)
    mp1[j] = fitstat5.Parameter(1)
    mp1err[j] = fitstat5.ParError(1)
    
    h5_effruns.Write()
    j=j+1
    
histodir5.cd()
efffit_barrel = RpcRollEffFit()

mytree = ROOT.TTree('FitParameters_barrel', 'FitParameters_barrel')
mytree.Branch('efffit_barrel',efffit_barrel,'rollidtree/I:p0/F:p0err/F:p1/F:p1err/F:chi2/F:ndof/F')

for i in range(0,j):
    efffit_barrel.rollidtree = rollidtree[i]
    efffit_barrel.p0 = mp0[i] + 2012*mp1[i]
    efffit_barrel.p0err = math.sqrt(math.pow(mp0err[i],2) + math.pow(2012,2) * math.pow(mp1err[i],2))
    efffit_barrel.p1 = mp1[i]
    efffit_barrel.p1err = mp1err[i]
    efffit_barrel.chi2 = mchi2[i]
    efffit_barrel.ndof = mndf[i]
    
    mytree.Fill()
    
mytree.Print()
mytree.Write()

rfile_out.Close()


print "DONE"

print runNumIndexFirst, runNumIndexLast








