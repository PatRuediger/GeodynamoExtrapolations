#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 09:19:17 2014

@author: Patrick

Algortithms for data extraction and Extrapolation

"""

import datastructure2
#from datastructure import dot
from datastructure2 import Point3D
from datastructure2 import dot
#from datastructure2 import AvsUcdAscii
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
gRE=ndarray((96,96))
hRE=ndarray((96,96))
gICB=ndarray((96,96))
hICB=ndarray((96,96))
#earth radius
ar= 2.91
#inner core boundary radius
icb = 5.380000000000000E-001
ocb = 1.537983852128
# data Object
#data = AvsUcdAscii()
g_data = None
isOuterCore = False
   # isInnerCore = False
isMantle = False



def setData(data):
    global g_data
    g_data = data

def evalSHA(x,dt):
    global isOuterCore
    global isMantle            
    if ((toSpherical(x)._z)<=ocb) and ((toSpherical(x)._z)>icb):
        result = g_data.getValueKDTree(x,dt) 
        #if not isOuterCore:print(">>>>>>>>>>>>>>>>>>Going from Mantle to OC: Pos:",toSpherical(x)," Value:",toSpherical(result))
        #print("OC")
        isOuterCore = True
        isMantle = False
        return result
    else:
        result = sphericalHarmoAnalysis(x)
    #if not isMantle:print(">>>>>>>>>>>>>>>>>>Going from OC to Mantle: Pos:",toSpherical(x)," Value:",toSpherical(result))
    #print("Mantle")
    isMantle = True
    isOuterCore = False         
    return result

def useIGRFonly():
    gRE = gIGRF
    gICB = gIGRF
    hRE = hIGRF
    hICB = hIGRF
    return
    
def getGaussCoef(radius):
    
    """IGRF only for testing pruposes"""
    if (radius>=ocb):
       # print("RE used",radius)
        return gRE,hRE
    elif (radius <=icb):
      #  print("ICB used",radius)
        return gICB,hICB
    else:
        print("Error in radius",radius)

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
    degree=95
    g,h = getGaussCoef(v._z)
    Lp = SN_Legendre(math.cos(v._x),degree)
    dLp = deltaSN_Legendre(Lp,degree)
    """for l in range(1,n):
        for m in range(l+1):
    #m = 1
    #l = 1
            result._x+= -((ar/v._z)**(l +2))*(-math.sin(v._x))*deltaSN(m,l,math.cos(v._x))*(g[l][m]*math.cos(m*v._y)+h[l][m]*math.sin(m*v._y))
            result._y+= -(ar/(v._z))**(l +2)*SN(m,l,cos(v._x))*(-g[l][m]*m*math.sin(m*v._y)+h[l][m]*m*math.cos(m*v._y))
            result._z+= (l +1)*((ar/v._z)**(l +2))*SN(m,l,cos(v._x))*(g[l][m]*math.cos(m*v._y)+h[l][m]*math.sin(m*v._y))
    """
    """Using implicit Legendre"""
    for l in range(1,degree):
        for m in range(l+1):
            result._x+= -((ar/v._z)**(l +2))*(-math.sin(v._x))*dLp[l][m]*(g[l][m]*math.cos(m*v._y)+h[l][m]*math.sin(m*v._y))
            result._y+= -(ar/(v._z))**(l +2)*Lp[l][m]*(-g[l][m]*m*math.sin(m*v._y)+h[l][m]*m*math.cos(m*v._y))
            result._z+= (l +1)*((ar/v._z)**(l +2))*Lp[l][m]*(g[l][m]*math.cos(m*v._y)+h[l][m]*math.sin(m*v._y))

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

def p(m,l,x):
    """Legendre Poynomial Basis"""
    if m==0 and l==0 :
        return 1.0
    elif m==0 and l==1:
        return cos(x)
    else:
        return p(0,l-1,x) * float(2*l-1)/float(l) * cos(x)- p(0,l-2,x) * float(l-1)/float(l)

def SN_Legendre(x,degree):
    """Schmid Normalized Legendrefunction
       returns a 2D Array with evaluated Legendre functions at position x
       use as L(m,l,x) = p_sn[m,l]""" 
    p_sn=[[0 for xl in range(degree)]*degree for xl in range(degree)]
    """Normalization"""
    df=[1.0 for xl in range(degree+1)]    
    for m in range(1,degree):   
        #df.insert(m,1.0)
        for k in range(1,m):
            df[m]=df[m] * float(2*k-1) / float(2*k)
        df[m+1] = sqrt( 2.0 * df[m] * float(2*m+1) ) * cos(x)
        df[m] = sqrt(2.0*df[m])
        
        if( m < degree-1):
            for l in range(m+2,degree):
                df[l]=( cos(x) * float(2*l-1) * df[l-1] - sqrt( float( (l-1)*(l-1) - m*m )) * df[l-2] ) / sqrt( float( l*l - m*m ))
        
        for l in range(m,degree):
            p_sn[m][l]=df[l] * sin(x)**m
    return p_sn

def deltaSN_Legendre(p_sn,degree):
   """Schmid Normalized Legendrefunction
   returns a 2D Array with evaluated Legendre functions at position x
   use as L(m,l,x) = p_sn[m,l]""" 
   dp_sn=[[0 for xl in range(degree)]*degree for xl in range(degree)]
   """Normalization"""
   dp_sn[0][0] = 0.0
   for l in range(1,degree):
       dp_sn[0][l] = - sqrt( float(l*(l+1)/2) ) * p_sn[1][l]
   dp_sn[1][1]=p_sn[0][1]

   if degree < 2: return dp_sn
   for l in range(2,degree):
       dp_sn[1][l]= 0.5 * ( sqrt( float( 2*l*(l+1) ) ) * p_sn[0][l] - sqrt( float((l-1)*(l+2)) ) * p_sn[2][l] )
       dp_sn[l][l]= 0.5 * sqrt(float(2*l))*p_sn[l-1][l]
   if degree <3: return dp_sn
   for l in range(3,degree):
       for m in range(2,l-1):
           dp_sn[m][l] = 0.5* ( sqrt( float( (l+m)*(l-m+1) ) ) *p_sn[m-1][l] - sqrt( float( (l-m)*(l+m+1) ) ) *p_sn[m+1][l] )
   return dp_sn
   
   
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
    v=Point3D(0,0,0)
    v._x = x._z*math.sin(x._x)*math.cos(x._y)
    v._y = x._z*math.sin(x._x)*math.sin(x._y)
    v._z = x._z*math.cos(x._x)
    return v
#Coordinate transformation Cartesian  to Spherical
def toSpherical(x):
    v=Point3D(0,0,0)
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
    vf=Point3D(0,0,0)
    vf._z= v._x*math.sin(x._x)*math.cos(x._y) + v._y*math.sin(x._x)*math.sin(x._y) + v._z*math.cos(x._x)
    vf._x = v._x*math.cos(x._x)*math.cos(x._y) + v._y*math.cos(x._x)*math.sin(x._y) - v._y*math.sin(x._x)
    vf._y = -v._x*math.sin(x._y) + v._y*math.cos(x._y)
    return vf

def toCartesianVecfield(x,v):
    """ Pos x has to be in spherical coordinates
    """
    vf=Point3D(0,0,0)
    vf._x = v._z*math.sin(x._x)*math.cos(x._y) + v._x*math.cos(x._x)*math.cos(x._y)- v._y*math.sin(x._y)
    vf._y = v._z*math.sin(x._x)*math.sin(x._y) + v._x*math.cos(x._x)*math.sin(x._y)+ v._y*math.cos(x._y)
    vf._z = v._z*math.cos(x._x) - v._x*math.sin(x._x)
    return vf


    
def adaptStep(v1,v2,dt):
    """
    Input Position, Vector, time(refering to v(t) as a vectorfield), current "time" t, an initial stepsize dt 
    Output boolean,adapted stepsize
    Domain specific knowledge regarding the Vectorfield can be added here, to speed up the estimation of the stepsize
    """
    #Error respective to cos(angle) with angle between the two steps 
    # 1%    = 0.9980267284282716
    # 0.1 % = 0.999980261
    err_up =  0.9999    
    err_down = 0.99999
    dtn = dt    
    #decrease stepsize when error is too high 
    if (datastructure2.dot(v1,v2)/( v1._length() * v2._length() ) ) < err_up: 
        #print(datastructure.dot(v1,v2)/(v1._length()*v2._length() ),"error too high")
        dtn/= 5.0
        #print("stepsize decreased",dtn)
        return True, dtn
    #increase stepsize when error is very small
    elif (datastructure2.dot(v1,v2)/(v1._length()*v2._length() ) ) > err_down:
     #   print(datastructure.dot(v1,v2)/(v1._length()*v2._length() ),"error too low")
        dtn*=2.0
        #print("stepsize increased",dtn)        
        return True, dtn
    else:
        return False, dtn

def rk4_backwards(x,v,a,t,dt):
    """Returns final (position, magnetic field) tuple after
    time dt has passed. In backwards order

    x: initial position (Point3D)
    v: initial magnetic field (Point3D)
    a: evaluation fucntion a(x,t) (must be callable) should return a vector valued item (e.g. trilinear interpolation)
    dt: timestep (number)"""
    return

def rk4(x, v, a, t, dt,direction="forward"):
    """Returns final (position, magnetic field) tuple after
    time dt has passed.
    x: initial position (Point3D)
    v: initial magnetic field (Point3D)
    a: evaluation fucntion a(x,t) (must be callable) should return a vector valued item (e.g. trilinear interpolation)
    dt: timestep (number)"""
    x0 = Point3D(x._x,x._y,x._z)
    v0 = Point3D(v._x,v._y,v._z)
    t1= t + dt
    
    """k1 should be equal to v0 """
    k1 = a(x0,t)
    k2 = a(x0.add(k1.mult(dt/2.0)),t + dt/2.0)
    k3 = a(x0.add(k2.mult(dt/2.0)),t + dt/2.0)
    k4 = a(x0.add(k3.mult(dt)),t + dt)
    
    if(direction == "backward"):
        """invert Vector"""
        k1 = k1.mult(-1.0)
        k2 = k2.mult(-1.0)
        k3 = k3.mult(-1.0)
        k4 = k4.mult(-1.0)
        #print("RK4 backwards used")
    
    k1 = k1.mult(1.0/6.0 * dt)    
    k2 = k2.mult(2.0/6.0 * dt)
    k3 = k3.mult(2.0/6.0 * dt)
    k4 = k4.mult(1.0/6.0 * dt)
    
    x1 = x0.add(k1)
    x1 = x1.add(k2)
    x1 = x1.add(k3)
    x1 = x1.add(k4)
    
    v1 = a(x1,t1)
    if(v1 ==None):
        print("v1 set back to v0' due to numerical error")
        print(k1,k2,k3,k4)
        v1=v0.mult(1.0-1e-12)
    ## steer the gaps inbetween here
    dt2 = dt/4.0
    
    t2 = t+dt2
    k1 = a(x0,t)
    k2 = a(x0.add(k1.mult(dt2/2.0)),t + dt2/2.0)
    k3 = a(x0.add(k2.mult(dt2/2.0)),t + dt2/2.0)
    k4 = a(x0.add(k3.mult(dt2)),t + dt2)
    
    if(direction == "backward"):
        """invert Vector"""
        k1 = k1.mult(-1.0)
        k2 = k2.mult(-1.0)
        k3 = k3.mult(-1.0)
        k4 = k4.mult(-1.0)
        #print("RK4 backwards used")
    
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
    #needAdapt = adaptStep(x1,v1,x2,v2,dt)[0]
    #while needAdapt:
        #needAdapt, dt = adaptStep(x1,v1,x2,v2,dt)
        #t1 = t+dt
        #x1,v1,t1,x2,v2,t2 = rk4(x0,v0,a,t1,dt)
        

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
    x0 = Point3D(x._x,x._y,x._z)
    v0 = Point3D(v._x,v._y,v._z)

    dt = step

    x1=x0.add(v0.mult(dt))
    x2=x0.add(v0.mult(dt/2))
    v1=a(x1,dt)
    v2=a(x2,dt/2)
        
    print(x1,v1)
    return Point3D(x1._x,x1._y,x1._z),Point3D(v1._x,v1._y,v1._z),dt
    

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
                for k in range(0,95):                    
                    for i in range (0,m+1):
                        #print("g",n,i,entries[row])
                        gRE[m][i]=entries[row]
                        row +=1                        
                    m +=1
                    if(k<95):
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
                for k in range(1,95):
                    for i in range (0,m+1):
                        #print("g",n,i,entries[row])
                        gICB[m][i]=entries[row]
                        row +=1
                    m +=1
                    if(k<95):
                        for q in range(m,0,-1):
                            #print("h",n,q,entries[row])
                            hICB[m][q]=entries[row]
                            row +=1
                count +=1
            else:
                #print("GausCoeffs for ICB Area read")
                break
    return gRE,hRE,gICB,hICB                

def main():
    """only dummy values here """
    #testSHA(1.0879, 0.0, -1.0879, 0,0,0)
    #testSHAwithVecfield()
    #plt.show()

if __name__ == "__main__":
    main()

