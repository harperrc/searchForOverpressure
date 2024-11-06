#!/usr/bin/python3

import math
import numpy as np
import os
from scipy import optimize as opt
import sys
from nuclib import *

R2D = 180.0 / math.pi 

def DPbyhb(x,gr,yld,dpd):
   dp = dynamicPressureFt(gr,x,yld)

   return (dp - dpd)

def getBnds(gr,yld,dpd):
   yg = []

   oldans = 1.0

   hb     = 0.0

   while (hb <= hbmax):
      dp = dynamicPressureFt(gr,hb,yld)
      ans = dp - dpd

      if (np.sign(oldans) != np.sign(ans)):
         yg.append(hb)

      oldans = ans
      
      hb = hb + 5.0 * stp

   if (len(yg) == 0):
      return (None)

   if (len(yg) == 1):
      tm = yg[0]
      yg = []
      yg = [0.0,tm]

   if (yg[0] == yg[1]):
      return (None)

   return (yg)

yld = 1.0

if (len(sys.argv) < 3):
   print('useage %s dp oname' % (sys.argv[0]))
   sys.exit()

dpd   = float(sys.argv[1])
aname = sys.argv[2]

oname = '%s.dp' % (aname)

if (dpd > 99.0):
   grmax = 450
   hbmax = 450
   stp   = 0.25
elif (dpd > 19.0):
   grmax = 650
   hbmax = 650
   stp   = 0.25
elif (dpd > 7.9):
   grmax = 950
   hbmax = 750
   stp   = 0.25
elif (dpd > 2.9):
   grmax = 1500
   hbmax =  950
   stp   = 0.50
elif (dpd > 0.9):
   grmax = 2400
   hbmax = 1250
   stp   = 1.00
elif (dpd > 0.39):
   grmax = 3250
   hbmax = 1500
   stp   = 1.00
else:
   print('no max for %10.3f' % (dpd))

gr = 0.00

x = []
y = []
a = []

maxgr = 0.0
maxhb = 0.0
gndb  = 0.0

while (gr <= grmax):
   yg = getBnds(gr,yld,dpd)

   if (yg == None):
      break

   i0 = 0
   i1 = 1

   nyg = len(yg) - 1

   for i in range(0,nyg):
      h0 = yg[i0]
      h1 = yg[i1]

      ans = opt.brentq(DPbyhb,h0,h1,args=(gr,yld,dpd))

      if (ans <= 1.0):
         gndb = max(gndb,gr)

      if (gr > maxgr):
         maxgr = gr
         maxhb = ans

      x.append(gr)
      y.append(ans)
      a.append(math.atan2(ans,gr) * R2D)

      i0 = i1
      i1 = i1 + 1

   gr = gr + stp

#  sort based on clockwise angle from (0,0)

ndx = np.argsort(a)

fo = open(oname,'w')

for (i,id) in enumerate(ndx):
   fo.write('%15.5f %15.5f\n' % (x[id],y[id]))
fo.close()

print('%15.5f %15.5f %15.5f %15.5f' % (dpd,maxgr,maxhb,gndb))
