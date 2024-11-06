#!/usr/bin/python3

#  use brents method to find root looping over ground range 
#  for a given height of burst all for a given overpressure

#  do it for 1Kt then use cube root scaling for others

import math
import numpy as np
import os
from scipy import optimize as opt
import sys
from nuclib import *
from other import *

def quartic(x,y):

#  takes 5 UNIFORMLY SPACED POINTS IN X and fits a quartic
#  and returns the root,  the y data MUST cross the zero
#  axis (e.g. negative and positive values)
#
#  inputs:
#    x     uniformly space dependant values
#    y     evaluated function at x
#
#  outputs:
#    value  value of x where y ~ 0
#
#  from Vallado Fundamentals of Astrodynamics and Applications

   accuracy = 1.0e-5
   value    = 0.0

#  this maps the user values to 0 to 4 space

   x1 = 0.0
   y1 = x[0]
   x2 = 4.0
   y2 = x[4]

   slope = (y2 - y1) / (x2 - x1)
   yint  = 0.50 * (y2 + y1) - 0.50 * slope * (x2 + x1)

#  solve the quartic

   c5 = y[0]
   c4 = (-50.00 * y[0] + 96.00 * y[1]   -72.00 * y[2] + \
               32.00 * y[3]   -6.00 * y[4]) / 24.00
   c3 = ( 35.00 * y[0] -104.00 * y[1] + 114.00 * y[2] - \
               56.00 * y[3] + 11.00 * y[4]) / 24.00
   c2 = (-10.00 * y[0] + 36.00 * y[1]   -48.00 * y[2] + \
               28.00 * y[3]   -6.00 * y[4]) / 24.00
   c1 = (         y[0]   -4.00 * y[1] +   6.00 * y[2] - \
               4.00 * y[3] +         y[4]) / 24.00

   q  = x1
   q0 = x2

#  iterate to solve the quartic
#  normally gets to accuracy within 3 loops

   n   = 0
   f   = 0.0
   fp  = 0.0
   fdp = 0.0

   while (abs(q-q0) > accuracy):
      q0 = q
      f   = c5 + q * (c4 + q * (c3 + q * (c2 + c1 * q)))
      fp  = c4 + q * (2.0 * c3 + q * (3.0 * c2 + 4.0 * c1 * q))
      fdp = 2.0 * c3 + q * (6.0 * c2 + 12.0 * c1 * q)

      n = n + 1
   
      if (n > 128):
         print('failed to find solution in quartic')
         value = -1.0e20

         return value

      q = q0 - f / (fp - f * fdp / (2.0 * fp))

#  convert the 0-4 value back to the user space

   value = yint + q * slope

   return value

def OPbygr(x,hob,yld,ovpd):
   o = idealOverpressureFt(x,hob,yld)
   return (o - ovpd)

def findMaxHB(yld,ovpd):

#  for 0 ground range find the maximum hob to 
#  acheive the desired overpressure

#  takes somewhat small steps then uses a
#  quartic to find the zero crossing

   hbmin =     0.000
   hbmax = 10000.000
   stp   = max(0.0001,1.0 / (2.0 * ovpd))

   npt   = int((hbmax - hbmin) / stp) + 1

   x  = [0.0] * npt
   y  = [0.0] * npt

   xx = [0.0] * 5
   yy = [0.0] * 5
   
   gr = 0.0
   
   nneg = 0
   
   hb   = hbmin
   
   for i in range(0,npt):
      ovp = OPbygr(gr,hb,1.0,0.0)
   
      dif  = ovp - ovpd

      if (dif < 0.0):
         nneg = nneg + 1
   
      x[i] = hb
      y[i] = dif
   
      if (nneg == 3):
         nfn = i
   
         if (nfn > 5):
            xx = x[nfn-5:nfn]
            yy = y[nfn-5:nfn]
   
            h  = quartic(xx,yy)
   
            return (h,stp)
   
      hb = hb + stp

   print('failed to find solution')
   sys.exit()

yld = 1.0

#  get inputs from user
#  ovp    overpressure in psi
#  aname  name for output files

ovp   = float(sys.argv[1])
aname = sys.argv[2]

#  find max hob at zero ground range

(hb,stp) = findMaxHB(yld,ovp)
hbmax    = math.ceil(hb) + 500.0
grmax    = 10000.0

pname = '%s.plt' % (aname)
oname = '%s.out' % (aname)

maxgr = 0.0
maxhb = 0.0

nf    = 0
hob   = 0.0

x = []
y = []

xmin = 1.0e20
xmax = -xmin
ymin = 1.0e20
ymax = -ymin

hobmax = 10000.0

while (hob <= hobmax):
   try:
      ans = opt.brentq(OPbygr,0,grmax,args=(hob,yld,ovp))

      nf  = nf + 1

      x.append(ans)
      y.append(hob)

      xmin = min(xmin,ans)
      xmax = max(xmax,ans)
      ymin = min(ymin,hob)
      ymax = max(ymax,hob)

      if (ans > maxgr):
         maxgr = ans
         maxhb = hob

   except ValueError as ve:

      xmin = math.floor(xmin)
      xmax = math.ceil(xmax)
      ymin = math.floor(ymin)
      ymax = math.ceil(ymax)

      xlab = looseLabel(xmin,xmax,ntick=11)
      ylab = looseLabel(ymin,ymax,ntick=11)

      print('%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f ' % \
            (ovp,maxgr,maxhb,
             xlab['min_val'],xlab['max_val'],
             ylab['min_val'],ylab['max_val']))

      if (x[-1] != 0.0):
         x.append(0.0)
         y.append(hb)

      x.reverse()
      y.reverse()

      fo = open(oname,'w')
      fo.write('2\n')
      fo.write('%d\n' % (len(x)))
   
      fp = open(pname,'w')
   
      for i in range(0,nf):
         fo.write('%15.5f %15.5f\n' % (x[i],y[i]))
         fp.write('%15.5f %15.5f\n' % (x[i],y[i]))

      fo.close()

      sys.exit()

   hob = hob + stp
