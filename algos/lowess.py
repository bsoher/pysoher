# -*- coding: utf-8 -*-
"""Smooth 1D data array using a Lowess algorithm  

The Lowess algorithm computes the smooth of a scatterplot of Y against X 
using robust locally weighted regression. Fitted values, YS, are computed 
at each of the values of the horizontal axis in X..

Functions:
    lowess (ndarray, float, float): ndarray
        Simplified, Pythonic wrapper around py_lowess
    py_lowess (...): 
        The original ported Fotran LOWESS funtion
    py_lowest (...): 
        Helper function for py_lowess. From original Fortran LOWEST funtion.

Example:

    yout = lowess(x, frac=0.1, delta=3.2)
        
        
Addional Info and Licenses: 
        
    This is our version of the Lowess code originally written in Fortran by 
    Cleveland, Grosse, and Shyu (see copyright below). It was ported from Fortran
    to Python by Brian J Soher 2012, 2019.

    Original Fortran Documentation and License:

        * The authors of this software are Cleveland, Grosse, and Shyu.
        * Copyright (c) 1989, 1992 by AT&T.
        * Permission to use, copy, modify, and distribute this software for any
        * purpose without fee is hereby granted, provided that this entire notice
        * is included in all copies of any software which is or includes a copy
        * or modification of this software and in all copies of the supporting
        * documentation for such software.
        * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
        * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
        * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
        * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 
        research!alice!wsc Mon Dec 30 16:55 EST 1985
        W. S. Cleveland
        ATT Bell Laboratories
        Murray Hill NJ 07974
 
        outline of this file:
        lines 1-72   introduction
            73-177   documentation for lowess
           178-238   ratfor version of lowess
           239-301   documentation for lowest
           302-350   ratfor version of lowest
           351-end   test driver and fortran version of lowess and lowest

        a multivariate version is available by "send dloess from a"
              
              COMPUTER PROGRAMS FOR LOCALLY WEIGHTED REGRESSION
             
             This package consists  of  two  FORTRAN  programs  for
        smoothing    scatterplots   by   robust   locally   weighted
        regression, or lowess.   The  principal  routine  is  LOWESS
        which   computes   the  smoothed  values  using  the  method
        described in The Elements of Graphing Data, by William S.
        Cleveland    (Wadsworth,    555 Morego   Street,   Monterey,
        California 93940).
             
             LOWESS calls a support routine, LOWEST, the code for
        which is included. LOWESS also calls a routine  SORT,  which
        the user must provide.
             
             To reduce the computations, LOWESS  requires  that  the
        arrays  X  and  Y,  which  are  the  horizontal and vertical
        coordinates, respectively, of the scatterplot, be such  that
        X  is  sorted  from  smallest  to  largest.   The  user must
        therefore use another sort routine which will sort X  and  Y
        according  to X.
             To summarize the scatterplot, YS,  the  fitted  values,
        should  be  plotted  against X.   No  graphics  routines are
        available in the package and must be supplied by the user.
             
             The FORTRAN code for the routines LOWESS and LOWEST has
        been   generated   from   higher   level   RATFOR   programs
        (B. W. Kernighan, ``RATFOR:  A Preprocessor for  a  Rational
        Fortran,''  Software Practice and Experience, Vol. 5 (1975),
        which are also included.
             
             The following are data and output from LOWESS that  can
        be  used  to check your implementation of the routines.  The
        notation (10)v means 10 values of v.

        
        X values:
          1  2  3  4  5  (10)6  8  10  12  14  50
        
        Y values:
           18  2  15  6  10  4  16  11  7  3  14  17  20  12  9  13  1  8  5  19

        YS values with F = .25, NSTEPS = 0, DELTA = 0.0
         13.659  11.145  8.701  9.722  10.000  (10)11.300  13.000  6.440  5.596
           5.456  18.998
        
        YS values with F = .25, NSTEPS = 0 ,  DELTA = 3.0
          13.659  12.347  11.034  9.722  10.511  (10)11.300  13.000  6.440  5.596
            5.456  18.998
        
        YS values with F = .25, NSTEPS = 2, DELTA = 0.0
          14.811  12.115  8.984  9.676  10.000  (10)11.346  13.000  6.734  5.744
            5.415  18.998

                                   
                                   LOWESS
        
        Calling sequence
        
        CALL LOWESS(X,Y,N,F,NSTEPS,DELTA,YS,RW,RES)
        
        Purpose
        
        LOWESS computes the smooth of a scatterplot of Y  against  X
        using  robust  locally  weighted regression.  Fitted values,
        YS, are computed at each of the  values  of  the  horizontal
        axis in X.
        
        Argument description
              
              X = Input; abscissas of the points on the
                  scatterplot; the values in X must be ordered
                  from smallest to largest.
              Y = Input; ordinates of the points on the
                  scatterplot.
              N = Input; dimension of X,Y,YS,RW, and RES.
              F = Input; specifies the amount of smoothing; F is
                  the fraction of points used to compute each
                  fitted value; as F increases the smoothed values
                  become smoother; choosing F in the range .2 to
                  idea which value to use, try F = .5.
         NSTEPS = Input; the number of iterations in the robust
                  fit; if NSTEPS = 0, the nonrobust fit is
                  returned; setting NSTEPS equal to 2 should serve
                  most purposes.
          DELTA = input; nonnegative parameter which may be used
                  to save computations; if N is less than 100, set
                  DELTA equal to 0.0; if N is greater than 100 you
                  should find out how DELTA works by reading the
                  additional instructions section.
             YS = Output; fitted values; YS(I) is the fitted value
                  at X(I); to summarize the scatterplot, YS(I)
                  should be plotted against X(I).
             RW = Output; robustness weights; RW(I) is the weight
                  given to the point (X(I),Y(I)); if NSTEPS = 0,
                  RW is not used.
            RES = Output; residuals; RES(I) = Y(I)-YS(I).

        
        Other programs called
               
               LOWEST
               SSORT
        
        Additional instructions
        
        DELTA can be used to save computations.   Very  roughly  the
        algorithm  is  this:   on the initial fit and on each of the
        NSTEPS iterations locally weighted regression fitted  values
        are computed at points in X which are spaced, roughly, DELTA
        apart; then the fitted values at the  remaining  points  are
        computed  using  linear  interpolation.   The  first locally
        weighted regression (l.w.r.) computation is carried  out  at
        X(1)  and  the  last  is  carried  out at X(N).  Suppose the
        l.w.r. computation is carried out at  X(I).   If  X(I+1)  is
        greater  than  or  equal  to  X(I)+DELTA,  the  next  l.w.r.
        computation is carried out at X(I+1).   If  X(I+1)  is  less
        than X(I)+DELTA, the next l.w.r.  computation is carried out
        at the largest X(J) which is greater than or equal  to  X(I)
        but  is not greater than X(I)+DELTA.  Then the fitted values
        for X(K) between X(I)  and  X(J),  if  there  are  any,  are
        computed  by  linear  interpolation  of the fitted values at
        X(I) and X(J).  If N is less than 100 then DELTA can be  set
        to  0.0  since  the  computation time will not be too great.
        For larger N it is typically not necessary to carry out  the
        l.w.r.  computation for all points, so that much computation
        time can be saved by taking DELTA to be  greater  than  0.0.
        If  DELTA =  Range  (X)/k  then,  if  the  values  in X were
        uniformly  scattered  over  the  range,  the   full   l.w.r.
        computation  would be carried out at approximately k points.
        Taking k to be 50 often works well.
        
        Method
        
        The fitted values are computed by using the nearest neighbor
        routine  and  robust locally weighted regression of degree 1
        with the tricube weight function.  A few additional features
        have  been  added.  Suppose r is FN truncated to an integer.
        Let  h  be  the  distance  to  the  r-th  nearest   neighbor
        from X(I).   All  points within h of X(I) are used.  Thus if
        the r-th nearest neighbor is exactly the  same  distance  as
        other  points,  more  than r points can possibly be used for
        the smooth at  X(I).   There  are  two  cases  where  robust
        locally  weighted regression of degree 0 is actually used at
        X(I).  One case occurs when  h  is  0.0.   The  second  case
        occurs  when  the  weighted  standard error of the X(I) with
        respect to the weights w(j) is  less  than  .001  times  the
        range  of the X(I), where w(j) is the weight assigned to the
        j-th point of X (the tricube  weight  times  the  robustness
        weight)  divided by the sum of all of the weights.  Finally,
        if the w(j) are all zero for the smooth at X(I), the  fitted
        value is taken to be Y(I).


                                   
                                   LOWEST
        
        Calling sequence
        
        CALL LOWEST(X,Y,N,XS,YS,NLEFT,NRIGHT,W,USERW,RW,OK)
        
        Purpose
        
        LOWEST is a support routine for LOWESS and  ordinarily  will
        not  be  called  by  the  user.   The  fitted  value, YS, is
        computed  at  the  value,  XS,  of  the   horizontal   axis.
        Robustness  weights,  RW,  can  be employed in computing the
        fit.
        
        Argument description
              
              X = Input; abscissas of the points on the
                  scatterplot; the values in X must be ordered
                  from smallest to largest.
              Y = Input; ordinates of the points on the
                  scatterplot.
              N = Input; dimension of X,Y,W, and RW.
             XS = Input; value of the horizontal axis at which the
                  smooth is computed.
             YS = Output; fitted value at XS.
          NLEFT = Input; index of the first point which should be
                  considered in computing the fitted value.
         NRIGHT = Input; index of the last point which should be
                  considered in computing the fitted value.
              W = Output; W(I) is the weight for Y(I) used in the
                  expression for YS, which is the sum from
                  I = NLEFT to NRIGHT of W(I)*Y(I); W(I) is
                  defined only at locations NLEFT to NRIGHT.
          USERW = Input; logical variable; if USERW is .TRUE., a
                  robust fit is carried out using the weights in
                  RW; if USERW is .FALSE., the values in RW are
                  not used.
             RW = Input; robustness weights.
             OK = Output; logical variable; if the weights for the
                  smooth are all 0.0, the fitted value, YS, is not
                  computed and OK is set equal to .FALSE.; if the
                  fitted value is computed OK is set equal to

        Method
        
        The smooth at XS is computed using (robust) locally weighted
        regression of degree 1.  The tricube weight function is used
        with h equal to the maximum of XS-X(NLEFT) and X(NRIGHT)-XS.
        Two  cases  where  the  program  reverts to locally weighted
        regression of degree 0 are described  in  the  documentation
        for LOWESS.

          test driver for lowess
          for expected output, see introduction
              real x(20), y(20), ys(20), rw(20), res(20)
              data x /1,2,3,4,5,10*6,8,10,12,14,50/
              data y /18,2,15,6,10,4,16,11,7,3,14,17,20,12,9,13,1,8,5,19/
              call lowess(x,y,20,.25,0,0.,ys,rw,res)
              write(6,*) ys
              call lowess(x,y,20,.25,0,3.,ys,rw,res)
              write(6,*) ys
              call lowess(x,y,20,.25,2,0.,ys,rw,res)
              write(6,*) ys
              end
  
"""

# Python modules
from __future__ import division

# 3rd party modules
import numpy

def lowess(y, frac = 0.1, delta = 3.2, nsteps=1):
    """
    Smooth 1D data array using a Lowess algorithm 

    The Lowess algorithm computes the smooth of a scatterplot of y against x 
    using robust locally weighted regression. Fitted values, ys, are computed 
    at each of the values of the horizontal axis in x.
    
    This wrapper assumes an even distribution of x values. Use py_lowess for
    more granular control of Lowess parameters.
    
    Args:
        y (ndarray): Y-values of data points on even distribution plot.
        frac (float, optional): default = 0.1, Specifies the amount of 
            smoothing; F is the fraction of points used to compute each fitted
            value; as F increases the smoothed values become smoother; 
            choosing F in the range 0.2 to idea which value to use, try 
            F = 0.5.
        delta (float, optional): default = 3.2, nonnegative parameter which 
            may be used to save computations; if n is less than 100, set delta
            equal to 0.0; if n is greater than 100 you should find out how 
            DELTA works by reading the additional instructions section for the
            py_lowess function.
        nsteps (int, optional): default = 1, the number of iterations in the 
            robust fit; if NSTEPS = 0 the nonrobust fit is returned; setting 
            NSTEPS equal to 1 or 2 should serve most purposes.
             
    Returns:
        yout (ndarray): fitted/smoothed value at each x location

    Example:
        yout = lowess(y, frac=0.3, delta=4.6, nsteps=3)    
    
    
    """
    ny    = numpy.size(y)
    yout  = numpy.zeros(ny, float)
    ysend = y
    
    # frac (# pts to smooth over), and delta control degree of smoothing
    if (frac >= 1.0) or (frac <= 0.0):
        frac = 0.1
    if delta > ny : delta = ny
    if delta < 1.0: delta = 1.0
    
    x   = numpy.arange(ny, dtype=float)
    # weights are in here if you want to take a look
    rw  = numpy.zeros(ny, float)
    # residuals are in here if you want to take a look
    res = numpy.zeros(ny, float)
    
    yout, rw, res = \
        py_lowess(x, ysend, ny, frac, nsteps, delta, yout, rw, res)
    
    return yout


def py_lowest(x,y,n,xs,ys,nleft,nright,userw,rw, w):
    """
    Python port of original helper function lowest in Fortran

    py_lowest is a support routine for py_lowess and ordinarily will not 
    be called by the user.  The fitted value, YS, is computed at the 
    value, XS, of the horizontal axis. Robustness weights, RW, can be 
    employed in computing the fit.
    
    Args:
        x (ndarray): X-corrdinates of the points on the scatterplot; the values 
            in X must be ordered from smallest to largest.
        y (ndarray): Y-values of the points on the scatterplot.
        n (int): dimension of x, y, w, and rw
        xs (ndarray): value of the horizontal axis at which the smooth 
            is computed.
        ys (): output, Left here for backward compatibility
        nleft (int): index of the first point which should be considered in 
            computing the fitted value.
        nright (int): index of the last point which should be considered in 
            computing the fitted value.
        userw (int): flag, if userw is 1, a robust fit is carried out using 
            the weights in rw; if userw is 0, the values in rw are not used.
        rw (ndarray): robustness weights.
        w (): output, Left here for backward compatibility
             
    Returns:
        ys (ndarray): fitted value at xs
        w  (ndarray): W(I) is the weight for Y(I) used in the expression for 
            YS, which is the sum from I = NLEFT to NRIGHT of W(I)*Y(I); W(I) 
            is defined only at locations NLEFT to NRIGHT.
        ok (int): flag, if the weights for the smooth are all 0.0, the 
            fitted value, YS, is not computed and OK is set equal to 
            0; if the fitted value is computed OK is set equal to 1

    Example:
        ys, w, ok = py_lowest(x,y,n,xs,ys,nleft,nright,userw,rw,w)

    
    """
    xrange = x[n-1]-x[0]
    hleft = xs-x[nleft-1]
    hright= x[nright-1]-xs
    h     = hleft * (hleft > hright) + hright * (hleft <= hright)
    h9    = 0.999*h
    h1    = 0.001*h
    a     = 0.0        # sum of weights
    w     = w * 0
    
    for j in range(nleft,n+1):     # compute weights (pick up all ties on right)
        w[j-1] = 0.0
        r = abs(x[j-1]-xs)
        
        if r <= h9:    # small enough for non-zero weight
            if r > h1:
                w[j-1] = (1.0 - (float(r)/h)**3)**3
            else:
                w[j-1] = 1.0
            
            if userw: w[j-1] = rw[j-1]*w[j-1]
            a = a+w[j-1]
        else:
            if x[j-1] > xs: break   # get out at first zero wt on right

    
    nrt = j        # rightmost pt (may be greater than nright because of ties)
    if a <= 0.0:
        ok = 0
    else:  # weighted least squares
        ok = 1
        w[nleft-1:nrt] = w[nleft-1:nrt]/a   # make sum of w(j) == 1
        
        if h > 0.0:    # use linear fit
            a = sum(w[nleft-1:nrt]*x[nleft-1:nrt]) # weighted center of x values
            b = xs-a
            c = sum(w[nleft-1:nrt]*(x[nleft-1:nrt]-a)**2)
            
            if numpy.sqrt(c) > 0.001*xrange:
                # points are spread out enough to compute slope
                b = b/c
                w[nleft-1:nrt] = w[nleft-1:nrt]*(1.0+b*(x[nleft-1:nrt]-a))
        
        ys = sum(w[nleft-1:nrt]*y[nleft-1:nrt])
    
    return ys, w, ok


def py_lowess(x,y,n,f,nsteps,delta,ys,rw,res):
    """
    Python port of original Fortran lowess algorithm. 

    LOWESS computes the smooth of a scatterplot of y against x using 
    robust locally weighted regression. Fitted values, ys, are computed 
    at each of the values of the horizontal axis in x.
    
    Calling sequence:
        
        ys, rw, res = py_lowess(x,y,n,f,nsteps,delta,ys,rw,res)

    Args:
        x (ndarray): X-corrdinates of the points on the scatterplot; the values 
            in X must be ordered from smallest to largest.
        y (ndarray): Y-values of the points on the scatterplot.
        n (int): dimension of x, y, w, and rw
        f (float): Specifies the amount of smoothing; F is the fraction of 
            points used to compute each fitted value; as F increases the 
            smoothed values become smoother; choosing F in the range 0.2 to 
            idea which value to use, try F = 0.5.
        nsteps (int): the number of iterations in the robust fit; if NSTEPS = 0
            the nonrobust fit is returned; setting NSTEPS equal to 2 should 
            serve most purposes.
        delta (float): nonnegative parameter which may be used to save 
            computations; if n is less than 100, set delta equal to 0.0; if 
            n is greater than 100 you should find out how DELTA works by 
            reading the additional instructions section.
        
    Returns:
        ys (ndarray): fitted values; ys(i) is the fitted value at x(i); to 
            summarize the scatterplot, ys(i) should be plotted against x(i).
        rw (ndarray): robustness weights; rw(i) is the weight given to the 
            point (x(i),y(i)); if nsteps = 0, rw is not used.
        res (ndarray): residuals; RES(I) = Y(I)-YS(I).  

    Additional instructions

        DELTA can be used to save computations. Very roughly the algorithm is 
        this:  on the initial fit and on each of the NSTEPS iterations 
        locally weighted regression fitted values are computed at points in 
        X which are spaced, roughly, DELTA apart; then the fitted values at 
        the remaining points are computed using linear interpolation.  The 
        first locally weighted regression (l.w.r.) computation is carried out 
        at X(1) and the last is carried out at X(N). Suppose the l.w.r. 
        computation is carried out at X(I).  If X(I+1) is greater than or 
        equal to X(I)+DELTA, the next l.w.r. computation is carried out at 
        X(I+1).  If X(I+1) is less than X(I)+DELTA, the next l.w.r. 
        computation is carried out at the largest X(J) which is greater than 
        or equal to X(I) but is not greater than X(I)+DELTA. Then the fitted 
        values for X(K) between X(I) and X(J), if there are any, are computed 
        by linear interpolation of the fitted values at X(I) and X(J). If N 
        is less than 100 then DELTA can be set to 0.0 since the computation 
        time will not be too great. For larger N it is typically not necessary 
        to carry out the l.w.r. computation for all points, so that much 
        computation time can be saved by taking DELTA to be greater than 0.0.
        If DELTA = Range (X)/k then, if the values in X were uniformly 
        scattered over the range, the full  l.w.r. computation would be 
        carried out at approximately k points. Taking k to be 50 often works 
        well.
        
    """
    ok = 0
    
    w = numpy.zeros(n,float)  # create once for many uses in py_lowest
    
    if n < 2:
        ys[0] = y[0]
        return ys, rw, res
    
    # at least two, at most n points
    ns = int(f*n)
    ns = numpy.where(ns < n-1, ns, n-1)
    ns = numpy.where(ns > 2, ns, 2)
    
    for iter in range(1,nsteps+2): # robustness iterations
        nleft  = 1
        nright = ns
        last   = 0        # index of prev estimated point
        i      = 1        # index of current point
        
        while True:
            while nright < n:
                # move nleft, nright to right if radius decreases
                d1 = x[i-1]-x[nleft-1]
                d2 = x[nright]-x[i-1]
                
                # if d1<=d2 with x(nright+1)==x(nright), lowest fixes
                if (d1 <= d2): break
                
                # radius will not decrease by move right
                nleft  = nleft+1
                nright = nright+1
            
            xi  = x[i-1]
            ysi = ys[i-1]
            
            ysi, res, ok = py_lowest(x,y,n,xi,ysi,nleft,nright,(iter > 1),rw, w)
            
            x[i-1]  = xi
            ys[i-1] = ysi
            
            # fitted value at x(i)
            if not ok: ys[i-1] = y[i-1]
            
            # all weights zero - copy over value (all rw==0)
            if last < i-1: # skipped points -- interpolate
                denom = x[i-1]-x[last-1]    # non-zero - proof?
                alpha = (x[last:i-1]-x[last-1])/float(denom)
                ys[last:i-1] = alpha*ys[i-1]+(1.0-alpha)*ys[last-1]
            
            last = i                  # last point actually estimated
            cut  = x[last-1]+delta    # x coord of close points
            
            for i in range(last+1,n+1):      # find close points
                if x[i-1] > cut: break       # i one beyond last pt within cut
                if x[i-1] == x[last-1]:      # exact match in x
                    ys[i-1] = ys[last-1]
                    last = i
            
            if i == n:
                i += 1
            i = numpy.where(last+1 > i-1, last+1, i-1)
            
            # back 1 point so interpolation within delta, but always go forward
            if last >= n: break
        
        res = y-ys
        
        if iter > nsteps: break    # compute robustness weights except last time

        rw = abs(res)
        rw = sorted(rw)
        
        m1   = int(1+n/2)
        m2   = n-m1+1
        cmad = 3.0*(rw[m1-1]+rw[m2-1])      # 6 median abs resid
        c9   = 0.999*cmad
        c1   = 0.001*cmad
        
        for i in range(n):
            r = abs(res[i])
            
            if r <= c1:
                rw[i] = 1.0        # near 0, avoid underflow
            else:
                if r > c9:
                    rw[i] = 0.0    # near 1, avoid underflow
                else:
                    rw[i] = (1.0-(r/cmad)**2)**2
    
    return ys, rw, res


#--------------------------------------------------------------------

def _test():
    
    import numpy as np

    n = 64
    f = 0.1
    nsteps = 1
    delta = 3.2
    ys = numpy.zeros(n)
    rw = numpy.zeros(n)
    res = numpy.zeros(n)
    x = numpy.arange(n)

    y=[ -0.836854, -0.172280,   0.187117,   1.61544,  -0.176774,    0.653145, 
        -0.546364,  0.194146,   0.925709,   1.20432,   1.53055,    -1.35556,    
         0.0514889, 1.02018,   -1.22616,    0.708497,  0.871673,   -0.789721,  
         0.332079,  0.205603,  -0.169367,  -0.318417, -0.295643,    0.522291,
        -2.23105,   0.258274,  -0.0877757, -1.64685,   0.286812,    0.299986,  
         1.10391,   0.742706,  -0.157581,  -0.597687,  0.659809,   -0.0328137,  
         1.16512,  -1.04800,    0.817815,  -1.40729,   0.519207,   -0.733439,  
         0.325304, -0.0428672, -0.871454,   0.771570,  0.00988832, -0.894773,
        -0.649426, -0.00869429, 1.87727,   -2.47856,  -1.68368,    -0.764296,  
         0.145749,  0.221329,   2.39680,    1.91128,   1.69614,     0.808025,   
         2.78748,   0.0275070,  0.104824,   1.78011]

    y = numpy.array(y)
    
    yout1            = lowess(y,f,delta)
    yout2, rw2, res2 = py_lowess(x,y,n,f,nsteps,delta,ys,rw,res)

    print "Max abs difference (lowess - py_lowess) = " + str(max(abs(yout1 - yout2)))

    
    # this is lowess unit test result for y array above
    y0=[-0.77791489, -0.36606162,  0.04579164,  0.4576449,   0.35414759,  0.25065028,
         0.14715296,  0.40698497,  0.66681697,  0.92664898,  0.6450437,   0.36343843,
         0.08183315,  0.12254237,  0.16325159,  0.20396081,  0.14572313,  0.08748545,
         0.02924777, -0.02934921, -0.08794618, -0.14654316, -0.18142722, -0.21631129,
        -0.25119536, -0.25733781, -0.26348027, -0.26962273,  0.01774429,  0.3051113,
         0.59247832,  0.40541144,  0.21834456,  0.03127769,  0.08988462,  0.14849156,
         0.2070985,   0.07388618, -0.05932614, -0.19253847, -0.17348405, -0.15442962,
        -0.1353752,  -0.12315045, -0.11092569, -0.09870093, -0.18450075, -0.27030057,
        -0.35610039, -0.58174107, -0.80738175, -1.03302243, -0.70402216, -0.37502189,
        -0.04602162,  0.52206765,  1.09015692,  1.6582462,   1.45195496,  1.24566373,
         1.0393725,   1.03840063,  1.03742877,  1.0364569 ]
         
    y0 = numpy.array(y0)
    
    print "Max abs difference (lowess - ground_truth, >1e-3 is good) = " + str(max(abs(yout1 - y0)))

    

if __name__ == '__main__':
    
    _test()

