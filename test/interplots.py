#!/usr/bin/python

import sys
#sys.argv.append('-b')
import os, commands
import math
import ROOT

from ROOT import *

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gStyle.SetOptFit(1)
ROOT.gStyle.SetOptStat(1)


c1 = ROOT.TCanvas()

rfile_in = ROOT.TFile.Open("RB4RollsEffNormFit.root")

mytree = rfile_in.Get("T")

mytree.Print()

mytree.Draw("chi2")
c1.SetLogy()
c1.SaveAs("rb4norm_chi2all.png")
mytree.Draw("chi2/ndof")
c1.SaveAs("rb4norm_chi2rall.png")
mytree.Draw("chi2/ndof","chi2/ndof<1000")
c1.SaveAs("rb4norm_chi2rlt1k.png")
mytree.Draw("chi2/ndof","chi2/ndof<100")
c1.SaveAs("rb4norm_chi2rlt100.png")

c1.SetLogy(0)

mytree.Draw("refp","chi2/ndof<10")
c1.SaveAs("rb4norm_refp_chi2rlt10.png")
mytree.Draw("p1","chi2/ndof<10")
c1.SaveAs("rb4norm_p1_chi2rlt10.png")

mytree.Draw("refp:p1","chi2/ndof<1000","box")
c1.SaveAs("rb4norm_refpvsp1_chi2rlt1k.png")
mytree.Draw("refp:p1","chi2/ndof<100","box")
c1.SaveAs("rb4norm_refpvsp1_chi2rlt100.png")
mytree.Draw("refp:p1","chi2/ndof<10","box")
c1.SaveAs("rb4norm_refpvsp1_chi2rlt10.png")
mytree.Draw("refp:p1","chi2/ndof<2","box")
c1.SaveAs("rb4norm_refpvsp1_chi2rlt2.png")

os.system("mv rb4norm_sel rb4norm_sel_OLD")

os.system("mkdir rb4norm_sel")

os.system("mv rb4norm*.png rb4norm_sel/")

