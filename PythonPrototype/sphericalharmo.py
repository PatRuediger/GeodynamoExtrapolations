#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 09:19:17 2014

@author: Patrick

Algortithms for data extraction and Extrapolation

"""

import datastructure
#from datastructure import dot
from datastructure import Point3D
from datastructure import AvsUcdAscii
from numpy import *
import math
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
import OpenGLVisFunc as Vis
import NeHeGL


#List which contains the evaluated datapoints from the SHA
plotPos_x=[]
plotPos_z=[]
plotPos_y=[]
plotValue_x=[]
plotValue_y=[]
plotValue_z=[]
plotColor=[]

#Gauss coefficients with max degree m
n=4
gIGRF=ndarray((6,6))
hIGRF=ndarray((6,6))
gRE=ndarray((6,6))
hRE=ndarray((6,6))
gICB=ndarray((6,6))
hICB=ndarray((6,6))
#earth radius
ar= 2.91
#inner core boundary radius
icb = 5.380000000000000E-001
# data Object
data = AvsUcdAscii()

# x,y,z should be list of number type objects
def drawStreamLine(x,y,z,color):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for x0,y0,z0,c0 in [(x,y,z,color)]:
        ax.scatter(x0, y0, z0, c=c0,cmap=mpl.cm.gray)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')



def drawVectorField(x,y,vx,vy):
    fig = plt.figure()
    ax= fig.add_subplot(111)
    ax.scatter(x,y)
    ax.quiver(x,y,vx,vy)

#computes a whole streamline
def streamline(start,t_max,step):
    """
    @Input:start postion (Point3D), max number of timesteps (Integer), stepsize (float)
    @Output:Returns a list of tuples, evaluated points from the rk4 method which can be plottet afterwards
    Last item in the list can be uses as SHA input
    """
    t=0
    sl=[]
    vl=[]
    nextPos = start
    nextVal = data.getValue(nextPos)
    while t < t_max:
        if data.getValue(start) is None:
            print("Leaving data domain ....... start Extrapolation")
            sl.append(rk4(nextPos,nextVal,evalSHA,t)[0])
            vl.append(rk4(nextPos,nextVal,evalSHA,t)[1])
        else:
            ## Inside data domain
            sl.append(rk4(nextPos,nextVal,evalf,t)[0])
            vl.append(rk4(nextPos,nextVal,evalf,t)[1])
        nextPos = sl[len(sl)-1]
        nextVal = vl[len(vl)-1]
        t+=step
    return sl

def evalf(x,dt):
    return data.getValue(x)

def evalSHA(x,dt):
    return sphericalHarmoAnalysis(x)

def getGaussCoef(radius):
    
    """IGRF only for testing pruposes"""
    if (radius>ar):
      #  print("IGRF used")
        return gIGRF,hIGRF
    elif (radius<=ar and (radius >(ar-icb)*0.8)):
     #   print("RE used")
        return gRE,hRE
    elif (radius <=(ar-icb)*0.8):
     #   print("ICB used")
        return gICB,hICB
    else:
        print("Error in radius")

def sphericalHarmoAnalysis(x):
    """Evaluates the spherical harmonics equation for the magnetic field, with degree n
    @Input: Point3D Position
    @Output: Point3D Magnetic Field
    """
    ##Coord Transformation
    #print("Pos Cartesian: "+ str(x._x)+","+str(x._y)+","+str(x._z))
    #print(x)
    v = toSpherical(x)
    #print(v)
    #print("Pos Spherical: "+ str(v._x)+","+str(v._y)+","+str(v._z))
    result=Point3D(0,0,0)

    """Get Gaus Coef respective to radius"""  
    g,h = getGaussCoef(v._z)
    
    for l in range(1,n):
        for m in range(l+1):
    #m = 1
    #l = 1
            result._x+= -((ar/v._z)**(l +2))*(-math.sin(v._x))*deltaSN(m,l,math.cos(v._x))*(g[l][m]*math.cos(m*v._y)+h[l][m]*math.sin(m*v._y))
            result._y+= -(ar/(v._z))**(l +2)*SN(m,l,cos(v._x))*(-g[l][m]*m*math.sin(m*v._y)+h[l][m]*m*math.cos(m*v._y))
            result._z+= (l +1)*((ar/v._z)**(l +2))*SN(m,l,cos(v._x))*(g[l][m]*math.cos(m*v._y)+h[l][m]*math.sin(m*v._y))

    """    for l in range(1,n+1):
        for m in range(l+1):
            result._x+=((a/v._z)**(l+1))*deltaSN(m,l,math.cos(v._x))*(g[l][m]*math.cos(m*v._y)+h[l][m]*math.sin(m*v._y))
            #print(a/v._z, (l+1), deltaSN(m,l,math.cos(v._x)), g[l][m], math.cos(m*v._y), h[l][m], math.sin(m*v._y))
            result._y+=((a/v._z)**(l+1))*SN(m,l,cos(v._x))*(-g[l][m]*m*math.sin(m*v._y)+h[l][m]*m*math.cos(m*v._y))
            result._z+=(l+1)*((a/v._z)**(l+2))*SN(m,l,cos(v._x))*(g[l][m]*math.cos(m*v._y)+h[l][m]*math.sin(m*v._y))
            #print(l,m,result)"""
#    print("sha,spherical:",result)
    vf = toCartesianVecfield(v,result)
#    print("sha,cartesian:",vf)
    return vf

#Schmidt-Normalized Legendrefunction
def SN(m,l,x):
    #print("compute Schmidt- Normalized Legendre " + str(l)+" , " + str(m)
    if(math.isnan(Legendre(m,l,x))): print("NaN error")
    if(Legendre(m,l,x)==0): print("is Zero error")
    return ((-1)**m)*math.sqrt(2.0*math.factorial(l-m)/math.factorial(l+m))*Legendre(m,l,x)
def deltaSN(m,l,x):
    #print("compute derivative Schmidt- Normalized Legendre " + str(l)+" , " + str(m))
    if(math.isnan(deltaLegendre(m,l,x))): print("NaN error")
    if(deltaLegendre(m,l,x)==0): print("is Zero error")
    return ((-1)**m)*math.sqrt(2.0*math.factorial(l-m)/math.factorial(l+m))*deltaLegendre(m,l,x)


#Legendre Functions till degree 4
def Legendre(m,l,x):
    if l==1 and m==0:
        return x
    elif l==1 and m==1:
  #      if x**2 >1.0:
  #          return 0.0
  #      else:
            return -math.sqrt(1.0-x**2.0)
    elif l==2 and m ==0:
        return 0.5*(3.0*x**2.0 -1.0)
    elif l==2 and m==1:
  #      if x**2 > 1.0:
  #          return 0.0
  #      else:
            return -3.0*x*math.sqrt(1.0-x**2.0)
    elif l==2 and m ==2:
        return 3.0*(x**2.0 -1.0)*(-1.0)
    elif l==3 and m ==0:
        return 0.5*(5.0*x**3.0 -3.0*x)
    elif l==3 and m==1:
   #     if x**2 >1.0:
   #         return 0.0
   #     else:
            return(-3.0/2.0)*math.sqrt(1.0-x**2.0)*(5.0*x**2.0 -1.0)
    elif l==3 and m ==2:
        return -15.0*x*(x**2.0 -1.0)
    elif l==3 and m ==3:
        return -15.0*(1.0-x**2.0)**(3.0/2.0)
    elif l==4 and m ==0:
        return (1.0/8.0)*(35.0*x**.04 -30.0*x**2.0 +3.0)
    elif l==4 and m==1:
   #     if x**2 >1.0:
   #         return 0.0
   #    else:
            return -(5.0/2.0)*math.sqrt(1.0-x**2.0)*(7.0**3.0 -3.0*x)
    elif l==4 and m==2:
        return -(15.0/2.0)*(x**2.0 -1.0)*(7.0*x**2.0 -1.0)
    elif l==4 and m==3:
        return 105.0*x*(1.0-x**2.0)**(3.0/2.0)
    elif l==4 and m==4:
        return 105*(x**2.0 -1.0)**2.0

#First derivative of Legende Functions
def deltaLegendre(m,l,x):
    if l==1 and m == 0:
        return 1.0
    elif l==1 and m==1:
#        if x**2 > 1.0:
#            return 0.0
#        else:
            return x/math.sqrt(1.0-x**2.0)
    elif l==2 and m ==0:
        return 3.0*x
    elif l==2 and m==1:
 #       if x**2 >1.0:
 #           return 0.0
 #      else:
            return (6.0*x**2.0 -3.0)/(math.sqrt(1.0-x**2.0))
    elif l==2 and m ==2:
        return -6.0*x
    elif l==3 and m ==0:
        return (3.0/2.0) *(5.0*x**2.0 -1.0)
    elif l==3 and m==1:
 #       if x**2 > 1.0:
 #           return 0.0
 #       else:
            return (3.0*x*(15.0*x*2.0 -11.0))/(2.0*math.sqrt(1.0-x**2.0))
    elif l==3 and m ==2:
        return 15.0 -45.0*x**2.0
    elif l==3 and m ==3:
 #       if x**2 > 1.0:
 #          return 0.0
 #       else:
            return 45.0*x*math.sqrt(1.0-x**2.0)
    elif l==4 and m ==0:
        return (5.0/2.0)*x*(7.0*x**-3.0)
    elif l==4 and m==1:
 #       if x**2 > 1.0:
 #           return 0.0
 #       else:
            return (5.0*(28.0*x**4.0 - 27.0*x**2.0 + 3.0))/(2*math.sqrt(1.0-x**2.0))
    elif l==4 and m==2:
        return -30.0*x*(7.0*x**2.0 -4.0)
    elif l==4 and m==3:
 #       if (x**2) > 1.0:
 #           return 0.0
 #       else:
            return 105.0*math.sqrt(1.0-x**2.0)*(4.0*x**2.0 -1.0)
    elif l==4 and m==4:
        return 420.0*x*(x**2.0 -1.0)

#Coordinate transformation Spherical (theta,phi,r)  to Cartesian (x,y,z)
def toCartesian(x):
    v=datastructure.Point3D(0,0,0)
    v._x = x._z*math.sin(x._x)*math.cos(x._y)
    v._y = x._z*math.sin(x._x)*math.sin(x._y)
    v._z = x._z*math.cos(x._x)
    return v
#Coordinate transformation Cartesian  to Spherical
def toSpherical(x):
    v=datastructure.Point3D(0,0,0)
    #theta
    v._x = math.acos(x._z/x._length())
    #phi
    v._y = math.atan2(x._y,x._x)
    #r
    v._z = x._length()
    return v


def toSphericalVecfield(x,v):
    """ Pos x has to be in spherical coordinates
    """
    vf=datastructure.Point3D(0,0,0)
    vf._z= v._x*math.sin(x._x)*math.cos(x._y) + v._y*math.sin(x._x)*math.sin(x._y) + v._z*math.cos(x._x)
    vf._x = v._x*math.cos(x._x)*math.cos(x._y) + v._y*math.cos(x._x)*math.sin(x._y) - v._y*math.sin(x._x)
    vf._y = -v._x*math.sin(x._y) + v._y*math.cos(x._y)
    return vf

def toCartesianVecfield(x,v):
    """ Pos x has to be in spherical coordinates
    """
    vf=datastructure.Point3D(0,0,0)
    vf._x = v._z*math.sin(x._x)*math.cos(x._y) + v._x*math.cos(x._x)*math.cos(x._y)- v._y*math.sin(x._y)
    vf._y = v._z*math.sin(x._x)*math.sin(x._y) + v._x*math.cos(x._x)*math.sin(x._y)+ v._y*math.cos(x._y)
    vf._z = v._z*math.cos(x._x) - v._x*math.sin(x._x)
    return vf

#Method for the parrallel SHA to check if the 2 evaluated points are in a certain range
def rangeCheck(v1,v2):
    return True
    
def adaptStep(x1,v1,x2,v2,dt):
    """
    Input Position, Vector, time(refering to v(t) as a vectorfield), current "time" t, an initial stepsize dt 
    Output boolean,adapted stepsize
    Domain specific knowledge regarding the Vectorfield can be added here, to speed up the estimation of the stepsize
    """
    #Error respective to cos(angle) with angle between the two steps 
    # 1%    = 0.9980267284282716
    # 0.1 % = 0.9999802608561371
    err =  0.9999802608561371
    dtn = dt    
    #decrease stepsize when error is too high 
    if (datastructure.dot(v1,v2)/(v1._length()*v2._length() ) ) < err: 
        #print(datastructure.dot(v1,v2)/(v1._length()*v2._length() ),"error too high")
        dtn/= 4.0
        return True, dtn
    #increase stepsize when error is very small
        """    elif (datastructure.dot(v1,v2)/(v1._length()*v2._length() ) ) > (err + (1.0-err)*0.9):
     #   print(datastructure.dot(v1,v2)/(v1._length()*v2._length() ),"error too low")
        dtn*=100.0        
        return True, dtn
        """        
    else :
    #    print(datastructure.dot(v1,v2)/(v1._length()*v2._length() ) ,"error 0.1%")
        return False, dtn

def rk4(x, v, a, t, dt):
    """Returns final (position, magnetic field) tuple after
    time dt has passed.

    x: initial position (Point3D)
    v: initial magnetic field (Point3D)
    a: evaluation fucntion a(x,t) (must be callable) should return a vector valued item (e.g. trilinear interpolation)
    dt: timestep (number)"""
    x0 = datastructure.Point3D(x._x,x._y,x._z)
    v0 = datastructure.Point3D(v._x,v._y,v._z)
    t1= t + dt
    
    """k1 should be equal to v0 """
    k1 = a(x0,t)
    k2 = a(x0.add(k1.mult(dt/2.0)),t + dt/2.0)
    k3 = a(x0.add(k2.mult(dt/2.0)),t + dt/2.0)
    k4 = a(x0.add(k3.mult(dt)),t + dt)
    
    k1 = k1.mult(1.0/6.0 * dt)    
    k2 = k2.mult(2.0/6.0 * dt)
    k3 = k3.mult(2.0/6.0 * dt)
    k4 = k4.mult(1.0/6.0 * dt)
    
    x1 = x0.add(k1)
    x1 = x1.add(k2)
    x1 = x1.add(k3)
    x1 = x1.add(k4)
    
    v1 = a(x1,t1)
    
    ## steer the gaps inbetween here
    dt2 = dt/4.0
    
    t2 = t+dt2
    k1 = a(x0,t)
    k2 = a(x0.add(k1.mult(dt2/2.0)),t + dt2/2.0)
    k3 = a(x0.add(k2.mult(dt2/2.0)),t + dt2/2.0)
    k4 = a(x0.add(k3.mult(dt2)),t + dt2)
    
    k1 = k1.mult(1.0/6.0 * dt2)    
    k2 = k2.mult(2.0/6.0 * dt2)
    k3 = k3.mult(2.0/6.0 * dt2)
    k4 = k4.mult(1.0/6.0 * dt2)
    
    x2 = x0.add(k1)
    x2 = x2.add(k2)
    x2 = x2.add(k3)
    x2 = x2.add(k4)
    
    v2 = a(x2,t2)
    
    #adapt stepsize
    needAdapt = adaptStep(x1,v1,x2,v2,dt)[0]
    while needAdapt:
        needAdapt, dt = adaptStep(x1,v1,x2,v2,dt)
        t1 = t+dt
        x1,v1,t1,x2,v2,t2 = rk4(x0,v0,a,t1,dt)
        

   # print(x1,v1,t1)
    return x1, v1, t1, x2, v2, t2

def eulerForward(x,v,a,step=1):
    """Returns final (position, magnetic field,stepsize) triple after
    time dt has passed.
    x: initial position (Point3D)
    v: initial magnetic field (Point3D)
    a: evaluation fucntion a(x,v,dt) (must be callable) should return a number like object (e.g. trilinear interpolation)
    step: initial timestep, only usefull, if it takes very long to find the first step (number)"""
    #error treshold for adaptive stepsize
    err = 0.99999999
    x0 = datastructure.Point3D(x._x,x._y,x._z)
    v0 = datastructure.Point3D(v._x,v._y,v._z)

    dt = step

    x1=x0.add(v0.mult(dt))
    x2=x0.add(v0.mult(dt/2))
    v1=a(x1,dt)
    v2=a(x2,dt/2)
        
    print(x1,v1)
    return Point3D(x1._x,x1._y,x1._z),Point3D(v1._x,v1._y,v1._z),dt
    
def testSHA(x,y,z,mx,my,mz):
    loadGaussCoefIGRF("../GausCoef.txt")
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")

    """for theta in range(10,2*3141,300):
        for phi in range(10,3141,300):
            tpr = datastructure.Point3D(theta/1000.0, phi/1000.0, 292/100.0)
            xyz = toCartesian(tpr)
            plotPos_x.append(xyz._x)
            plotPos_y.append(xyz._y)
            plotPos_z.append(xyz._z)
            v = evalSHA(xyz, None)
            plotValue_x.append(v._x)
            plotValue_y.append(v._y)
            plotValue_z.append(v._z)
            vf = toSphericalVecfield(tpr,v)
            plotColor.append(1.0)
          #  print(tpr,vf,v)"""
    
    sl=[]
    vl=[]
    tmax = 1000.0
    #initial stepsize, optional
    step = 0.8e-2
    max_steps = 4000
    t=step
    tpr0 = datastructure.Point3D(0.75*math.pi, 0.5*math.pi,2.9)
    xyz0 = toCartesian(tpr0)
    v0 = evalSHA(xyz0, None)
    nextVal=v0
    nextPos=xyz0
    plotPos_x.append(xyz0._x)
    plotPos_y.append(xyz0._y)
    plotPos_z.append(xyz0._z)
    print(xyz0,v0)
    i= 0

        
    while (max_steps>i) and (t<tmax):
       # print(((t-step)*100.0)/(tmax-step),"% finished ..... ")
        if(toSpherical(nextPos)._z<icb):
            print(toSpherical(nextPos))
            print("Inner Core reached")
            break
       # print(toSpherical(nextPos));
        xf, vf ,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,evalSHA,t,step)
      #  print(xf,vf,t2)
        #xfs=toSpherical(xf)
        #vfs=toSphericalVecfield(xfs,vf)
        plotPos_x.append(xf._x)
        plotPos_y.append(xf._y)
        plotPos_z.append(xf._z)
        plotValue_x.append(vf._x)
        plotValue_y.append(vf._y)
        plotValue_z.append(vf._z)
        plotColor.append(0.5)
        sl.append(xf)
        vl.append(vf)
        nextPos = xf
        nextVal = vf
        t=t2
        print("step #",i)
        i+=1     
    
    #vf = toSphericalVecfield(tpr,v)
    #print(tpr,vf,v)
    Vis.setStreamLine(sl,vl)
    NeHeGL.main()
    #drawStreamLine(plotPos_x,plotPos_y,plotPos_z,plotColor)
    return

def testSHAwithVecfield():
    loadGaussCoefIGRF("../GausCoef.txt")
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    
    """define Vectorfield Positions in spherical coordinates"""
    xf = []
    for theta in range(10,2*3141,315):
        for phi in range(10,3141,158):
            tpr = datastructure.Point3D(theta/1000.0, phi/1000.0, 2.0)
            xyz = toCartesian(tpr)
            xf.append(xyz)
          #print(tpr._x, tpr._y, tpr._z)
          
    #calculate VF Values
    vf = []
    for x in xf:
        v = evalSHA(x,1.0)
        print("Pos:",x,"Value:",v)
        vf.append(v)
 #   vf = []    
   # tpr=datastructure.Point3D(3.141,0.1,2.92)
   # xyz= toCartesian(tpr)
   # xf.append(xyz)
  #  vf.append(evalSHA(tpr,1.0))
 #   print(xf)
  #  print(vf)
  #  value_n = datastructure.Point3D((vf[0]._x/vf[0]._length),(vf[0]._y/vf[0]._length),(vf[0]._z/vf[0]._length))
  #  print(vf[0]._x/vf[0]._length())
    Vis.setVectorfield(xf,vf)
    NeHeGL.main()
    return

def loadGaussCoefIGRF(filename):
    with open(filename,'r') as file:
        count = 0
        for line in file:
            if(count<2):
                ##skip first 2 lines
                doNothing=0
                count +=1
            else:
                entries = line.split()
                if entries[0] is 'g':
                    gIGRF[int(entries[1])][int(entries[2])] =float(entries[3])
                    count +=1
                elif entries[0] is 'h':
                    hIGRF[int(entries[1])][int(entries[2])] =float(entries[3])
                    count +=1
                else :
                    print("Error: Invalid File Format")

def loadGaussCoefSimu(filenameRE,filenameICB):
    with open(filenameRE,'r') as file:
        count = 0
        for line in file :
            ##skip as many lines, as needed to get to you disered timestep
            if(count<3):
                ##skip first 3 lines
                doNothing=0
                count +=1
            ## point to the line with the disered timestep    
            elif (count==3):
                entries = line.split()
                ##we only use Coefs until the degree of 5
                m = 0
                row = 2                
                for k in range(0,6):                    
                    for i in range (0,m+1):
                        #print("g",n,i,entries[row])
                        gRE[m][i]=entries[row]
                        row +=1                        
                    m +=1
                    if(k<5):
                        for q in range(m,0,-1):
                            #print("h",n,q,entries[row])
                            hRE[m][q]=entries[row]
                            row +=1
                count +=1
            else:
                #print("GausCoeffs for Mantel Area read")
                break
                
    with open(filenameICB,'r') as file:
        count = 0
        for line in file :
            ##skip as many lines, as needed to get to you disered timestep
            if(count<3):
                ##skip first 3 lines
                doNothing=0
                count +=1
            ## point to the line with the disered timestep    
            elif (count==3):
                entries = line.split()
                ##we only use Coefs until the degree of 5
                m = 0
                row = 2
                for k in range(1,6):
                    for i in range (0,m+1):
                        #print("g",n,i,entries[row])
                        gICB[m][i]=entries[row]
                        row +=1
                    m +=1
                    if(k<5):
                        for q in range(m,0,-1):
                            #print("h",n,q,entries[row])
                            hICB[m][q]=entries[row]
                            row +=1
                count +=1
            else:
                #print("GausCoeffs for ICB Area read")
                break
                

def main():
    """only dummy values here """
    #testSHA(1.0879, 0.0, -1.0879, 0,0,0)
    testSHAwithVecfield()
    #plt.show()

if __name__ == "__main__":
    main()

