#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 09:19:17 2014

@author: Patrick

Algortithms for data extraction and Extrapolation

"""

import datastructure
from datastructure import Point3D
from datastructure import AvsUcdAscii
from numpy import *
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib as mpl
#from mayavi import mlab



# List which contains the streamline datapoints of the forward/backward rungeKutta
fw_points=[]
bw_points=[]

#List which contains the evaluated datapoints from the SHA
sha_points=[]
plotPos_x=[]
plotPos_z=[]
plotPos_y=[]
plotValue_x=[]
plotValue_y=[]
plotValue_z=[]
plotColor=[]

#Gauss coefficients with max degree m
n=4
g=ndarray((6,6))
h=ndarray((6,6))
#earth radius
#a = 6371000
ar= 2.91
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

def evalf(x,v,dt):
    return data.getValue(x)

def evalSHA(x,v,dt):
    return sphericalHarmoAnalysis(x)

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
    #print("compute Schmidt- Normalized Legendre " + str(l)+" , " + str(m))
    return ((-1)**m)*math.sqrt(2.0*math.factorial(l-m)/math.factorial(l+m))*Legendre(m,l,x)
def deltaSN(m,l,x):
    #print("compute derivative Schmidt- Normalized Legendre " + str(l)+" , " + str(m))
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

def rk4(x, v, a, dt):
    """Returns final (position, magnetic field) tuple after
    time dt has passed.

    x: initial position (Point3D)
    v: initial magnetic field (Point3D)
    a: evaluation fucntion a(x,v,dt) (must be callable) should return a number like object (e.g. trilinear interpolation)
    dt: timestep (number)"""
    x1 = datastructure.Point3D(x._x,x._y,x._z)
    v1 = datastructure.Point3D(v._x,v._y,v._z)
    a1 = a(x1, v1, 0)

    x2 = x.add(v1.mult(0.5*dt))
    v2 = v.add(a1.mult(0.5*dt))
    a2 = a(x2, v2, dt/2.0)

    x3 = x.add(v2.mult(0.5*dt))
    v3 = v.add(a2.mult(0.5*dt))
    a3 = a(x3, v3, dt/2.0)

    x4 = x.add(v3.mult(dt))
    v4 = v.add(a3.mult(dt))
    a4 = a(x4, v4, dt)

    print(v1, v2, v3, v4)
    print(a1, a2, a3, a4)
    xf = x
    xf=xf.add(v1.mult(dt/6))
    xf=xf.add(v2.mult(dt/6*2))
    xf=xf.add(v3.mult(dt/6*2))
    xf=xf.add(v4)
#    xf = x + (dt/6.0)*(v1 + 2*v2 + 2*v3 + v4)
    vf = v.add(a1.mult(dt/6))
    vf=vf.add(a2.mult(dt/6*2))
    vf=vf.add(a3.mult(dt/6*2))
    vf=vf.add(a4)
#    vf = v + (dt/6.0)*(a1 + 2*a2 + 2*a3 + a4)

    print("Pos: " + str(xf._x)+ ","+ str(xf._y)+ "," + str(xf._z))
    print("Val: " + str(vf._x)+ ","+ str(vf._y)+ "," + str(vf._z))
    return xf, vf

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
    v1=a(x1,v0,dt)
    v2=a(x2,v0,dt/2)
#    print("x0,v0:" , x0,v0)
#    print("x1,v1:",x1,v1)
#    print("x2,v2:",x2,v2)
    #adaptive stepsize

#    print("Angle between:",datastructure.dot(v1,v2)/(v1._length() * v2._length()))
    """
    while (datastructure.dot(v1,v2)/(v1._length() * v2._length())) < err :
        dt/=2
        x1=x0.add(v0.mult(dt))
        x2=x0.add(v0.mult(dt/2))
        v1=a(x1,v0,dt)
        v2=a(x2,v0,dt/2)
        """
        
    print(x1,v1)
    return Point3D(x1._x,x1._y,x1._z),Point3D(v1._x,v1._y,v1._z),dt
    
def testSHA(x,y,z,mx,my,mz):
    loadGaussCoef("../GausCoef.txt")

    for theta in range(10,2*3141,300):
        for phi in range(10,3141,300):
            tpr = datastructure.Point3D(theta/1000.0, phi/1000.0, 292/100.0)
            xyz = toCartesian(tpr)
            plotPos_x.append(xyz._x)
            plotPos_y.append(xyz._y)
            plotPos_z.append(xyz._z)
            v = evalSHA(xyz, None, None)
            plotValue_x.append(v._x)
            plotValue_y.append(v._y)
            plotValue_z.append(v._z)
            vf = toSphericalVecfield(tpr,v)
            plotColor.append(1.0)
          #  print(tpr,vf,v)
            
    #drawVectorField(plotPos_x,plotPos_y,plotValue_x,plotValue_y)
    #drawStreamLine(plotPos_x,plotPos_y,plotPos_z,plotColor)        
    #return
    
    sl=[]
    vl=[]
    tmax = 6.63e-03
    #initial stepsize, optional
    step = 5.0e-06
    max_steps = 10000
    t=step
    tpr0 = datastructure.Point3D(0.75*math.pi, 0.5*math.pi, 292/100.0)
    xyz0 = toCartesian(tpr0)
    v0 = evalSHA(xyz0, None, None)
    nextVal=v0
    nextPos=xyz0
    plotPos_x.append(xyz0._x)
    plotPos_y.append(xyz0._y)
    plotPos_z.append(xyz0._z)
    print(xyz0,v0)
    i= 0
    while (max_steps>i) and (t<tmax):
        print(((t-step)*100.0)/(tmax-step),"% finished ..... ")
        xf, vf ,t2= eulerForward(nextPos,nextVal,evalSHA,step)
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
        t+=t2
        print("step #",i)
        i+=1
    #vf = toSphericalVecfield(tpr,v)
    #print(tpr,vf,v)  
    drawStreamLine(plotPos_x,plotPos_y,plotPos_z,plotColor)
    return


def loadGaussCoef(filename):
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
                    g[int(entries[1])][int(entries[2])] =float(entries[3])
                    count +=1
                elif entries[0] is 'h':
                    h[int(entries[1])][int(entries[2])] =float(entries[3])
                    count +=1
                else :
                    print("Error: Invalid File Format")


def main():
    testSHA(1.0879, 0.0, -1.0879, 0,0,0)
    plt.show()

if __name__ == "__main__":
    main()

