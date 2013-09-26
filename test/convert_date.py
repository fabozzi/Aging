#!/usr/bin/python

import sys
#sys.argv.append('-b')
import os, commands
import math
import ROOT

from datetime import datetime as dt
import time

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

runtimeDict = {}

runtimeFile = open("rpcrunlist_date.txt", "r")
runepochFile = open("rpcrun_epoch.txt", "w")
for r_entry in runtimeFile:
#    print r_entry.rstrip()
     ent = r_entry.rstrip().split('\t')
     rundate = ent[1].split('/')
     runday = int(rundate[0])
     runmonth = int(rundate[1])
     runyear = int(rundate[2]) + 2000
     runtime = ent[2].split(':')
     runh = int(runtime[0])
     runm = int(runtime[1])
     runs = int(runtime[2])
     print runday, runmonth, runyear, runh, runm, runs
     run_dt = dt(year=runyear, month=runmonth, day=runday, hour=runh, minute=runm, second=runs)
     print toYearFraction(run_dt)
     runepochFile.write( ent[0] +"\t"+ str(toYearFraction(run_dt)) +"\n" )
runepochFile.close()

