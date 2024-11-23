#!/usr/bin/python3

import math
import numpy as np
import sys

m2ft = 100.0 / 30.48

def psiToGetPd_P(pd,vn):

# pd 0->1
# vn p type vn (adjusted)

   psi = math.exp(pd / 1.12820) * 0.719784 * 1.20**vn

   return (psi)

def psiToGetPd_Q(pd,vn):

# pd 0->1
# vn q type vn (adjusted)

   psi = math.exp(pd / 0.32200) * 6.10868e-3 * 1.44**vn

   return (psi)

def adjustPVN(yld,vn,k):

#  yld Kt
#  vn  vulnerability number
#  k   adjustment factor

   if (k > 0):
      vnn   = vn
      part1 = 0.10 * k
      part2 = (20.0 / yld)**0.33333333333333
      b     = -part1 * part2
      c     = part1 - 1.0
      z     = 0.50 * (-b + math.sqrt(b * b - 4.0 * c))
      vnn   = vnn + math.log(z * z) / 0.182321556
   else:
      vnn = vn

   return (vnn)

def adjustQVN(yld,vn,k):

#  yld Kt
#  vn  vulnerability number
#  k   adjustment factor

   if (k > 0):
      vnn   = vn
      part1 = 0.10 * k
      part2 = (20.0 / yld)**0.33333333333333333
      a3    = 1.0
      a2    = 0.0
      a1    = -part1 * part2
      a0    = part1 - 1.0
      p     = a2 / a3
      q     = a1 / a3
      r     = a0 / a3

      ylo   = 0.0
      yhi   = 1000.0
    
      z     = rootc(p,q,r,ylo,yhi)
      vnn   = vnn + math.log(z * z * z) / 0.364643
   else:
      vnn   = vn

   return (vnn)

def DPCoeff(M,y):

#  M Xq / x
#  y kft/yld**0.33333333

#  equations from PSR report 1419_3 for peak horizontal dynamic pressure
#  pg 121

   L = math.log10(M)
   A = -236.1 + (17.72 * M**0.593) / (1.0 + 10.4 * M**3.124)
   B = 12.27 - (21.69 * M**2.24) / (1.0 + 6.976 * M**0.484)
   C = 20.26 + (14.7 * M**2) / (1.0 + 0.08747 * M**3.05)
   D = -1.137 - (0.5606 * M**0.895) / (1.0 + 3.406 * M**7.48)
   E = 1.731 + (10.84 * M**1.12) / (1.0 + 12.26 * M**0.0014)
   F = 2.84 + (0.855 * M**0.9) / (1.0 + 1.05 * M**2.84)
   G = 50.0 - (1843.0 * y**2.153) / (1.0 + 3.95 * y**5.08)
   H = 0.294 + (71.56 * y**8.7) / (1.0 + 115.3 * y**6.186)
   I = abs(-3.324 + (987.5 * y**4.77) / (1.0 + 211.8 * y**5.166))
   J = 1.955 + (169.7 * y**9.317) / (1.0 + 97.36 * y**6.513)
   K = 8.123e-6 + (0.001613 * y**6.428) / (1.0 + 60.26 * y**7.358)

   return (A,B,C,D,E,F,G,H,I,J,K,L)

def dynamicPressureFt(grft_,hobft_,yld_):

#  grkm    ground range from burst (Ft)
#  hobkm   height of burst (Ft)
#  yld     yield (Kt)

# equations from PSR report 1419_3 for peak horizontal dynamic pressure

   yld    = yld_
   yld13  = yld**0.33333333333333

#  convert to Kft

   grkft  = grft_ / 1000.0
   hobkft = hobft_ / 1000.0

   x      = max(1.0e-8,grkft / yld13)
   y      = max(1.0e-8,hobkft / yld13)
   r      = math.sqrt(x * x + y * y)

   Xq     = (63.5 * y**7.26) / (1.0 + 67.11 * y**4.746) + 0.6953 * y**0.808

   M = max(1.0e-6,Xq / x)

   (A,B,C,D,E,F,G,H,I,J,K,L) = DPCoeff(M,y)
 
   J = min(150,max(-150,J))
   E = min(150,max(-150,E))

   Qs = (A * r**D) / (1.0 + B * r**E) + C / r**F

   if (x >= Xq):
      return (Qs)

# if in mach region then evaluate above term first at M = 1

   (A1,B1,C1,D1,E1,F1,G1,H1,I1,J1,K1,L1) = DPCoeff(1.0,y)

   rp = math.sqrt(Xq * Xq + y * y)
   Qm = (A1 * rp**D1) / (1.0 + B1 * rp**E1) + C1 / rp**F1

# apply other terms

   t1 = (G * L**I) / (1.0 + 649.0 * L**I)
   t2 = (4.01 * L**J) / (1.0 + H * L**J)
   t3 = 7.67e-6 * (1.0 / (K + L**3.22) - 1.0 / K)
   arg = t1 - t2 + t3

   Qs = Qm * math.exp(arg)

   return (Qs)

def dynamicPressureKm(grkm_,hobkm_,yld_):
   grft_  = grkm_ * 3280.8
   hobft_ = hobkm_ * 3280.8

   dp = dynamicPressureFt(grft_,hobft_,yld_)

   return (dp)

def idealOverpressureFt(grft,hobft,yld):

#  PSR report 1419-3 (better at high overpressure than previous)

#  grft   ground range from burst to target (ft)
#  hobft  height of burs (ft)
#  yld    yield of burst (Kt)

   y13 = yld**0.3333333333333333

   x     = max(1.0e-9,grft / y13)
   y     = max(1.0e-9,hobft / y13)

   capr  = math.sqrt(x * x + y * y)
   capr2 = capr * capr

   r     = capr / 1000.0
   oz    = x / y
   z     = y / x

   y2  = y * y
   y3  = y * y2
   y4  = y * y3
   y5  = y * y4
   y6  = y * y5
   y7  = y * y6

   z2  = z * z
   z3  = z * z2
   z4  = z * z3
   z5  = z * z4
   z6  = z * z5
   z7  = z * z6
   z8  = z * z7
   z9  = z * z8
   z17 = z9 * z8
   z18 = z9 * z9

   a  = 1.22 - 3.908 * z2 / (1.0 + 810.2 * z5)

   b  = 2.321 + 6.195 * z18 / (1.0 + 1.113 * z18) - \
        0.03831 * z17 / (1.0 + 0.02415 * z17) + \
        0.6692 / (1.0 + 4164.0 * z8)

   bb = 0.0629 * oz**8.34 / \
        (1.0 + 0.00509 * oz**13.05) * 0.05 * y /  \
        (1.0 + 2.56e-8 * y5)

   c  = 4.153 - 1.149 * z18 / (1.0 + 1.641 * z18) - \
        1.10 / (1.0 + 2.771 * z**2.50)

   d  = -4.166 + 25.76 * z**1.75 / \
         (1.0 + 1.382 * z18) + 8.257 * z /  \
         (1.0 + 3.219 * z)

   e  = 1.00 - 0.004642 * z18 / (1.0 + 0.003886 * z18)

   f  = 0.6069 + 2.879 * z**9.250 / \
        (1.0 + 2.359 * z**14.5) -  \
        17.15 * z2 / (1.0 + 71.66 * z3)

   g  = 1.83 + 5.361 * z2 / (1.0 + 0.3139 * z6)

   h = -(0.2905 + 64.67 * z5) / (1.0 + 441.5 * z5) - \
        1.389 * z / (1.0 + 49.03 * z5) + \
        (8.808 * z**1.5) / (1.0 + 154.5 * z**3.5) + \
        (1.094 * capr2) / ((0.7813e9 - 1.234e5 * capr + \
        1201.0 * capr**1.5 + capr2) * \
        (1.0 + 2.0 * y))
   
   p  = 1.8008e-7 * y4 / (1.0 + 0.0002863 * y4) - \
         2.121 * y2 / (794300.0 + y**4.3)

   q = 5.18 + 8.864 * y**3.5 / (3.788e6 + y4)

   ppeak = 10.47 / r**a + (b - bb) / r**c + d * e / (1.0 + f * r**g) + h + p / (r**q)

   return (ppeak)

def idealOverpressureKm(grkm,hobkm,yld):

#   this function computes the overpressure over an ideal surface
#   based on speicher and brodes work (pacific serria research)  
#   airblast from nuclear burst analytic approximations (PSR Report 1419-3)
#
#   inputs
#     grkm   ground range from target to burst (Km)
#     hobkm  height of burst (Km)
#     yld    yield of weapon (kt)
#
#   outputs   overpressure  - overpressure (psi)


#  convert to feet

   grft   = gfkm * m2ft * 1000.0
   hobftt = hobkm * m2ft * 1000.0

   ovp  = idealOverpressureFt(grft,hobft,yld)

   return (ovp)

def timar(grft,hobft,yld):

# PSR report 1419-3

#  input
#     grft    ground range from burst (ft)
#     hobft   height of burst (ft)
#     yld     yield of detonation (Kt)

#  output
#     tau     time of arrival (sec)

   y13 = yld**0.3333333333

#  compute scaled gr,height

   x = grft / y13
   y = hobft / y13

#  handle 0's

   x = max(x,1.0e-9)
   y = max(y,1.0e-9)

   capr = math.sqrt(x * x + y * y)

   r  = capr / 1000.0
   z  = y / x

   if (z > 100.0):
      z = 100.0

   r2 = r * r
   r3 = r2 * r

#  compute time of shock arrival

   xm  = 170.00 * y / (1.0 + 60.0 * y**0.25) + 2.89 * (y / 100.0)**2.5

   top = (0.543 - 21.8 * r + 386.0 * r2 + 2383.00 * r3) * r**8

   bot = 2.99e-14 - 1.91e-10 * r2 + 1.032e-6 * r**4 - \
          4.43e-6 * r**6 + (1.028 + 2.087 * r +  \
          2.69 * r2) * r**8

   u   = top / bot

#  N.B. tau in milliseconds so multiply by 0.001 when returning

   tau = u

   if (x < xm):
      return (tau * 0.001 * y13)

#  if in mach region

   top = (1.086 - 34.605 * r + 486.30 * r2 + 2383.0 * r3) * r**8
  
   bot = 3.0137e-13 - 1.2128e-9 * r2 + 4.128e-6 * r**4 - \
         1.116e-5 * r**6 + (1.632 + 2.629 * r + \
         2.69 * r2) * r**8

   w   = top / bot

   tau = u * xm / x + w * (1.0 - xm / x)

#  N.B. tau in milliseconds so multiply by 0.001 when returning

   return (tau * 0.001 * y13)

def rootc(p,q,r,ylo,yhi):
#
#  solve cubic in form x^3 + p * x^2 + q * x + r = 0

   third   = 0.33333333333
   pi2thd  = 2.0 * 3.141592654 / 3.0

   povr3   = p * third
   a       = q - p * povr3
   b       = p * p * p / 13.5 - povr3 * q + r
   alpha   = -b / 2.0
   gama    = a * third
   alphasq = alpha * alpha
   gamcub  = gama * gama * gama

   det     = alphasq + gamcub

   if (det == 0):
      root = -(alpha**third) - povr3
      if (ylo < root and yhi > root):
         return (root)

      root = -2.0 * root - p
      return (root)

   elif (det > 0.0):
      sdet = math.sqrt(det)
      arg1 = alpha + sdet
      arg2 = alpha - sdet

      s1 = 1.0
      if (arg1 != 0.0):
         s1 = arg1 / abs(arg1)

      s2 = 1.0
      if (arg2 != 0.0):
         s2 = arg2 / abs(arg2)

      root = s1 * (abs(arg1)**third) + s2 * (abs(arg2)**third) - povr3
      return (root)

   else:
      xm = 2.0 * math.sqrt(-gama)
      theta = 3.0 * b / a / xm
      theta = math.acos(theta) * third

      for i in range(0,3):
         root = xm * math.cos(theta) - povr3
         if (ylo < root and yhi > root):
            return (root)

         theta = theta + pi2thd

   return (0.0)

def pdam(grx,hobx,vn,yld):
#
#  return probability of damage for P type targets
#  grx    ground range (Km)
#  hobx   height of burst (Km)
#  vn     vulnerability number (if < 57 is vn if > is 100 * vn + k)
#  yld    weapon yield (kt)

   if (vn < 0.0):
      return (0.0)

   third  = 0.33333333333

   gr     = max(grx,1.0e-6)
   hob    = max(hobx,1.0e-6)

#  scale

   gr     = gr / yld**third
   hob    = hob / yld**third

   psi    = idealOverpressureKm(gr,hob,yld)
   psi    = max(0.0,psi)

   if (vn > 57.0):
      k     = int(vn) % 100
      vnn   = (vn - k) / 100.0
      part1 = 0.10 * k
      part2 = (20.0 / yld)**third
      b     = -part1 * part2
      c     = part1 - 1.0
      z     = 0.50 * (-b + math.sqrt(b * b - 4.0 * c))
      vnn   = vnn + math.log(z * z) / 0.182322
   else:
      vnn   = vn

   pd = 1.12820 * math.log(psi / (0.719784 * (1.20)**vnn));
   pd = min(1.0,max(0.0,pd))

   return (pd,psi)

def qdam(grx,hobx,vn,yld):
#
#  return probability of damage for Q type targets
#  grx    ground range (Km)
#  hobx   height of burst (Km)
#  vn     vulnerability number (if < 35 is vn if > is 100 * vn + k)
#  yld    weapon yield (kt)

   if (vn < 0.0):
      return (0.0)

   third  = 0.33333333333

   gr     = max(grx,1.0e-6)
   hob    = max(hobx,1.0e-6)

#  scale

   gr     = gr / yld**third
   hob    = hob / yld**third

   dp     = dynamicPressureKm(grx,hobx,yld)

   if (vn > 35):
      k     = int(vn) % 100
      vnn   = (vn - k) / 100.0
      part1 = 0.10 * k
      part2 = (20.0 / yld)**third
      a3    = 1.0
      a2    = 0.0
      a1    = -part1 * part2
      a0    = part1 - 1.0
      p     = a2 / a3
      q     = a1 / a3
      r     = a0 / a3

      ylo   = 0.0
      yhi   = 1000.0
    
      z     = rootc(p,q,r,ylo,yhi)
      vnn   = vnn + math.log(z * z * z) / 0.364643
   else:
      vnn   = vn

   pd = 0.3220 * math.log(dp / (6.10867e-3 * (1.44)**vnn))
   pd = min(1.0,max(0.0,pd))

   return (pd,dp)

def rangeFromYieldHOBOVP(overPressure,heightOfBurst,yieldOfBurst):
#
#  computes ground range given overpressure, height of burst and yield
#
#  ADA194873 "Direct Determination of Range From Current Nuclear Overpressure Equations"
#            Roger S Wolczek, March 1988
#            Air Force Institute Of Technology Air University
#
#  inputs:
#     overPressure (PSI)
#     heightOfBurst (Km)
#     yieldOfBurst  (Kt)

   MMC_metersToFeet = 3.2808398950131233595

   w    = yieldOfBurst
   h    = heightOfBurst * 1000.0 * MMC_metersToFeet
   ps   = overPressure

   m = w**0.33333333
   y = (h/m)/1000.0
   p = ps / 1000.0

   psq  = p * p
   pcub = psq * p
   pqt  = psq * psq
   p5   = pqt * p
   py   = p * y
   p2y  = psq * y
   p3y  = pcub * y
   p4y  = pqt * y
   p5y  = p5 * y
   ysq  = y * y
   ycub = ysq * y
   yqt  = ysq * ysq
   y5   = yqt * y
   pysq = psq * ysq
   pyc  = pcub * ycub
   pyqt = pqt * yqt
   p5y5 = p5 * y5
   y2p  = ysq * p
   y3p  = ycub * p
   y4p  = yqt * p
   y5p  = y5 * p
   p3y2 = pcub * ysq
   p4y2 = pqt * ysq
   p5y2 = p5 * ysq
   p4y3 = pqt * ycub
   p5y3 = p5 * ycub
   p5y4 = p5 * yqt
   y3p2 = ycub * psq
   y4p2 = yqt * psq
   y5p2 = y5 * psq
   y4p3 = yqt * pcub
   y5p3 = y5 * pcub
   y5p4 = y5 * pqt

   f = 1000.0 * m / 3280.8

   if (p > 100.0):
      return -1.0

   if (p < 0.001):
      return -1.0

   if (p > 10.0):
      x=0.09125435 - 0.00259593 * p - 0.011649 * py + 0.0000538051 * psq + \
        0.0001111895 * p2y + 1.55066730 * y2p - 0.0156667 * pysq - \
        5.2376e-07 * pcub + 0.00004916175 * p3y2 + 97.31778176 * ycub -\
        25.706 * y3p + 1.86726e-09 * pqt - 776.975 * yqt + \
        80.37406298 * y4p + 0.10802327 * y3p2 + 0.01220085 * y4p3 - \
        0.000106692 * pyqt

      return x * f

   if (p > 1.0):
      x = 0.20154762 - 0.0621137  * p + 0.26548571 * y + 0.01362242 *  psq - \
          7.29235 * ysq + 0.07036665 * pysq - 0.00141394 * pcub + \
          39.26591463 * ycub + 47.69269458 * y3p + 1.41551464 * pyc + \
          0.00005402714 * pqt - 0.0671128 * p4y3 - 314.529 * y4p - \
          12.1885 * y3p2 + 64.01100943 * y4p2 - 8.28738 * y4p3 + \
          0.42330598 * pyqt

      return x * f

   if (p > 0.10):
      x = 0.49576495 - 1.76963 * p + 0.2873705 * y + 4.11651645 * psq - \
          7.91446 * ysq + 30.88152847 * y2p - 157.536 * pysq - \
          4.38405 * pcub + 237.78055884 * p3y2 + 30.88701918 * ycub - \
          90.0484 * y3p - 1383.45 * pyc + 1.69350254 * pqt - \
          109.732 * p4y2 + 631.54888112 * p4y3 - 29.9186 * yqt + \
          876.59060473 * y3p2 - 1444.15 * y4p2 + 2183.11803 * y4p3 - \
          935.61 * pyqt

      return x * f

   if (p > 0.01):
      x = 1.62275752 - 83.2005 * p + 2705.11147 * psq + 6.12568509 * ysq - \
          546.803 * y2p + 10261.68969 * pysq - 45666.2 * pcub  - \
          56907.7 * p3y2 - 14.3825 * ycub + 1654.94706 * y3p + \
          375478.28531 * pqt + 15.10099190 * yqt - 1985.89 * y4p  - \
          22491.0 * y3p2 + 18051.29615 * y4p2 + 4798899.81 * pyqt - \
          1186924.0 * p5 + 8777484.91 * p5y3 - 48910776.0 * p5y4 - \
          5.23814 * y5 + 660.83710834 * y5p - 5735942.0 * y5p4 + \
          49158205.9 * p5y5 

      return x * f

   if (p > 0.001):
      x = 8.23926119 - 5306.71 * p + 6.05714948 * y - 5471.77 * py + \
          1828923.98 * psq + 1686307.36 * p2y - 2.89962 * ysq  + \
          4167.17167 * y2p - 352120.0 * pysq - 31.9297925e7 * pcub - \
          27.9705884e7 * p3y + 0.87118870 * ycub - 1187.34 * y3p  + \
          162.39402197e6 * pyc + 269.81456091e8 * pqt  + \
          240.17907175e8 * p4y - 103.60283247e8 * p4y3  - \
          219.221 * y4p - 855214.0 * y3p2 + 501833.53514 * y4p2  - \
          47.312304e6 * y4p3 - 8.75407e11 * p5 - 8.08007e11 * p5y + \
          3361.23611805e8 * p5y3 - 0.0104459 * y5 + 34.21894224 * y5p  - \
          31509.7 * y5p2 - 11250907.0 * y5p3 + 203.1223752e7 * y5p4  - \
          6059.3995968e7 * p5y5 

      return x * f

   return -1.0

def pt(gr,hob,yld,tasa):

#  input
#     gr    ground range from burst (ft)
#     hob   height of burst (ft)
#     yld   yield of detonation (Kt)
#     tasa  time after shock arrival (msec)

#  output
#     deltap overpressure (PSI)
#     tau    time of shock arrival at gr (msec)
#     capd   duration of positive phase (msec)

   y13 = yld**0.3333333333

#  compute scaled gr,height

   x     = gr / y13
   y     = hob / y13
   sigma = tasa / y13

#  handle 0's

   x = max(x,1.0e-9)
   y = max(y,1.0e-9)

   capr = math.sqrt(x * x + y * y)

   r  = capr / 1000.0
   z  = y / x

#  compute peak overpressure, change from fortran to pass acutal gr,hob & yld

   deltps = idealOverpressureFt(gr,hob,yld)
   
   if (z > 100.0):
      z = 100.0

   r2 = r * r
   r3 = r2 * r

#  compute time of shock arrival

   xm = 170.0 * y / (1.0 + 60.0 * y**0.25) + 2.89 * (y / 100.0)**2.5

   top = (0.543 - 21.8 * r + 386.0 * r2 + 2383.0 * r3) * r**8

   bot = 2.99e-14 - 1.91e-10 * r2 + 1.032e-6 * r**4 -  \
         4.43e-6 * r**6 + (1.028 + 2.087 * r +  \
         2.69 * r2) * r**8

   u   = top / bot

   tau = u

   if (x >= xm):

#  if in mach region

      top = (1.086 - 34.605 * r + 486.30 * r2 + 2383.0 * r3) * r**8

      bot = 3.0137e-13 - 1.2128e-9 * r2 + 4.128e-6 * r**4 - \
            1.116e-5 * r**6 + (1.632 + 2.629 * r + \
            2.69 * r2) * r**8

      w   = top / bot
      tau = u * xm / x + w * (1.0 - xm / x)

#  handle extremely small tau which can cause extreme underflow
#  when raised to the 8th power.... micro milli seconds seems short enough

   tau = max(tau,1.0e-6)

   sigma = sigma + tau

#  compute some things to save time

   y2 = y * y

   r2 = r * r
   r3 = r * r2
   r4 = r * r3
   r5 = r * r4
   r6 = r * r5
   r7 = r * r6
   r8 = r * r7

   tau2  = tau * tau

   s2 = 1.0 - 15.18 * ((y / 100.0)**3.5) /  \
        (1.0 + 15.18 * ((y / 100.0)**3.5)) -  \
        (0.02441 * ((y / 1.0e6)**2) /  \
        (1.0 + 9000.0 * ((y / 100.0))**7)) *  \
        (1.0e10 / (0.441 + (x / 100.00)**10))

#  duration of positive phase in milliseconds pg 75 (95/210)

   capd = ((1640700.0 + 24629.0 * tau + 416.15 * tau2) /  \
          (10800.0 + 619.76 * tau + tau2)) * (0.40 +  \
           0.001204 * (tau**1.5) / (1.0 +  \
           0.001559 * tau**1.5) +  \
           (0.6126 + 0.5486 * (tau**0.25) / \
          (1.0 + 0.003570 * tau**1.5) -  \
          3.47 * (tau**0.637) / \
          (1.0 + 5.696 * tau**0.645)) * s2)

   s = 1.0 - 1100.0 * ((y / 100.0)**7) /  \
       (1.0 + 1100.0 * ((y / 100.0)**7)) -  \
       (2.441e-14 * y2/ (1.0 + 9000.0 *  \
       ((y / 100.0)**7))) * (1.0e10 / \
       (0.441 + ((x / 100.0)**10)))

   f2 = (0.445 - 5.44 * (r**1.02) /  \
        (1.0 + 100000.0 * (r**5.84)) +  \
        7.571 * (z**7.15) / (1.0 - 5.135 * (z**12.9)) -  \
        8.07 * (z**7.31) / (1.0 + 5.583 * (z**12.23))) *  \
        (0.435 * ((y / 10.0)**1.26) / ( 1.0 +  \
        0.03096 * ((y / 10.0)**3.12))) *  \
        (1.0 - 0.000019 * (tau**8) / (1.0 + \
         0.000019 * (tau**8)))

   f  = (0.01477 * (tau**0.75) / (1.0 + 0.005836 * tau) +  \
         7.402e-5 * (tau**2.5) /  \
        (1.0 + 1.4290e-8 * (tau**4.75)) - 0.216) * s +  \
         0.7076 - 3.077e-5 * (tau**3) /  \
        (1.0 + 4.3670e-5 * (tau**3.00)) + f2 -  \
        (0.452 - 9.94e-7 * (x**4.130) / \
        (1.0 + 2.1868e-6 * (x**4.130))) *  \
        (1.0 - 1.5397e-4 * (y**4.3) / \
        (1.0 + 1.5397e-4 * (y**4.30)))

   g  = 10.0 + (77.58 - 64.99 * (tau**0.125) /  \
                 (1.0 + 0.04348 * (tau**5.0))) * s

   h  = 3.003 + 0.05601 * tau / (1.0 + 1.473e-9 * (tau**5)) +  \
        (0.01769 * tau /  \
        (1.0 + 3.207e-10 * (tau**4.25)) -  \
         0.03209 * (tau**1.25) / (1.0 + \
        9.914e-8 * (tau**4)) - 1.6) * s -  \
        0.1966 * (tau**1.22) / (1.0 + \
        0.767 * (tau**1.22))

   b  = (f * ((tau / sigma)**g) + (1.0 - f) *  \
        ((tau / sigma)**h)) *  \
        (1.0 - (sigma - tau) / capd)

   if ((x < xm) or (y > 380.0)):
      deltap = deltps * b

      return (deltap,tau,capd,xm)

   xe = 3.039 * y / (1.0 + 0.0067 * y)
   ak = abs((x- xm) / (xe - xm))

   if (ak > 50.0):
      ak = 50.0

   d2 = 2.99 + 31240.0 * ((y / 100.0)**9.86) /  \
        (1.0 + 15530.0 * ((y / 100.0)**9.87))

   d  = 0.23 + 0.583 * y2 / (26667.0 + y2) +  \
        0.27 * ak + (0.50 - 0.583 * y2 / \
        (26667.0 + y2)) * (ak**d2)

   a  = (d - 1.0) * (1.0 - (ak**20) / (1.0 + (ak**20)))

#   if (abs(x-xm) < 1.0e-10):
#      aj = 0.0
#   else:
#      aj = 11860.0 * (sigma - tau) / (y * ((x - xm)**1.25))
   aj = 11860.0 * (sigma - tau) / (y * ((x - xm)**1.25))

   if (aj > 200.0):
      aj = 200.0

   v  = 1.0 + (0.003744 * ((y / 10.0)**5.185) / \
        (1.0 + 0.004684 * ((y / 10.0)**4.189)) +  \
        0.004755 * ((y / 10.0)**8.049) /  \
        (1.0 + 0.003444 * ((y / 10.0)**7.497)) - \
        0.04852 * ((y / 10.0)**3.423) /  \
        (1.0 + 0.03038 * ((y / 10.0)**2.538))) * \
        (aj**3) / (6.13 + (aj**3)) * (1.0 /  \
        (1.0 + 9.23 * (ak**2)))

   c3 = 1.0 + (1.094 * (ak**0.738) /  \
        (1.0 + 3.687 * (ak**2.63)) *  \
        (1.0 - 83.01 * ((y / 100.0)**6.5) /  \
        (1.0 + 172.3 * ((y / 100.0)**6.04)) - 0.15)) * \
        (1.0 / (1.0 + 0.5089 * (ak**13)))

   c2 = 23000.0 * ((y / 100.0)**9) /  \
          (1.0 + 23000.0 * ((y / 100.0)**9))

   temp = (x / 100.0)**4

   c  = (1.04 - 0.02409 * temp / (1.0 + 0.02317 * temp)) *  \
        (aj**7) / ((1.0 + a) * (1.0 + 0.923 * (aj**8.5))) *  \
        (c2 + (1.0 - c2) * (1.0 - 0.09 * (ak**2.5) /  \
                             (1.0 + 0.09 * (ak**2.5)))) * \
        c3 * (1.0 - (((sigma - tau) / capd)**8))

   deltap = deltps * (1.0 + a) * (b * v + c)

   return (deltap,tau,capd,xm)

