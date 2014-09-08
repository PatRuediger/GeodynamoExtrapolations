# -*- coding: utf-8 -*-
"""
Created on Sat Sep 06 11:22:36 2014

@author: Patrick
"""
from sphericalharmo import *

def dipol_vf():
    return    

def test_Dipol_VF():
    loadGaussCoefIGRF("../GausCoef.txt")
    useIGRFonly() #overrites all other Coefficients
    """define Vectorfield Positions in spherical coordinates"""
    xf = []
    for theta in range(10,2*3141,315):
        for phi in range(10,3141,158):
            for radi in range(292,1000,120):
                tpr = datastructure.Point3D(theta/1000.0, phi/1000.0, radi/100.0)
                xyz = toCartesian(tpr)
                xf.append(xyz)          
    """calculate VF Values"""
    vf = []
    for x in xf:
        v = evalSHA(x,1.0)
        print("Pos:",x,"Value:",v)
        vf.append(v)
    Vis.setVectorfield(xf,vf)
    return

def testRK_Dipol_SL(theta,phi,r,tmax=1000.0,t0=0.8e-3,max_steps=10000):
    loadGaussCoefIGRF("../GausCoef.txt")
    useIGRFonly() #overrites all other Coefficients
    sl=[]
    vl=[]
    t=t0
    tpr0 = datastructure.Point3D(theta, phi,r)
    xyz0 = toCartesian(tpr0)
    xyz0 = xyz0.mult(10.0)
    v0 = evalSHA(xyz0, None)
    sl.append(xyz0.mult(0.1))
    vl.append(v0)
    nextVal=v0
    nextPos=xyz0
    print("Start new Streamline with:", tpr0,toSphericalVecfield(tpr0,v0))
    #print(tpr0,toSphericalVecfield(tpr0,v0))
    i= 0 
    while (max_steps>i) and (t<tmax):
       # print(((t-t0)*100.0)/(tmax-t0),"% finished ..... ")
        if(((toSpherical(nextPos)._z*0.1)<2.91) and i>1):
            print("Mantle reached",toSpherical(nextPos))
            break
        #print("next pos spherical :",toSpherical(nextPos));
        xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,evalSHA,t,t0)
        """adapt stepsize"""
        needAdapt= adaptStep(xf,vf,xf2,vf2,t0)[0]
        while needAdapt:
            needAdapt, t0 = adaptStep(xf,vf,xf2,vf2,t0)
            xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,evalSHA,t,t0)
        nextPos = xf
        nextVal = vf
        t=t2
        xfscaled = xf.mult(0.1)
        sl.append(xfscaled)
        vl.append(vf)
        #print("step #",i,toSpherical(xf), vf)
        i+=1
    Vis.addStreamLine(sl,vl)
    return

def testSHA_SL(x,y,z,mx,my,mz):
    #loadGaussCoefIGRF("../GausCoef.txt")
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    
    sl=[]
    vl=[]
    tmax = 1000.0
    #initial stepsize, optional
    step = 0.8e-2
    max_steps = 4000
    t=step
    tpr0 = datastructure.Point3D(0.75*math.pi, 0.5*math.pi,2.91)
    xyz0 = toCartesian(tpr0)
    v0 = evalSHA(xyz0, None)
    sl.append(xyz0)
    vl.append(v0)
    nextVal=v0
    nextPos=xyz0
    i= 0        
    while (max_steps>i) and (t<tmax):
       # print(((t-step)*100.0)/(tmax-step),"% finished ..... ")
        if(toSpherical(nextPos)._z<icb):
            print(toSpherical(nextPos))
            print("Inner Core reached")
            break
        xf, vf ,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,evalSHA,t,step)
        #adapt stepsize
        needAdapt = adaptStep(xf, vf,xf2,vf2,tf2)[0]
        while needAdapt:
            needAdapt, dt = adaptStep(xf, vf,xf2,vf2,tf2)
            t1 = t+dt
            x1,v1,t1,x2,v2,t2 = rk4(x0,v0,a,t1,dt)

        nextPos = xf
        nextVal = vf
        t=t2
        xfscaled = xf.mult(0.1)
        sl.append(xfscaled)
        vl.append(vf)
        #print("step #",i,xf,vf)
        i+=1     
    
    #vf = toSphericalVecfield(tpr,v)
    #print(tpr,vf,v)
    Vis.addStreamLine(sl,vl)
    return

def testSHAwithVecfield():
    loadGaussCoefIGRF("../GausCoef.txt")
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    
    """define Vectorfield Positions in spherical coordinates"""
    xf = []
    for theta in range(10,2*3141,315):
        for phi in range(10,3141,158):
            tpr = datastructure.Point3D(theta/1000.0, phi/1000.0, (ar-icb)*0.79)
            xyz = toCartesian(tpr)
            xf.append(xyz)
          #print(tpr._x, tpr._y, tpr._z)
          
    #calculate VF Values
    vf = []
    for x in xf:
        v = evalSHA(x,1.0)
        print("Pos:",x,"Value:",v)
        vf.append(v)
    Vis.setVectorfield(xf,vf)
    return
    
def main():
    """only dummy values here """
    #testSHA_SL(1.0879, 0.0, -1.0879, 0,0,0)
    #testSHAwithVecfield()
    for phi in range(10,2*3141,500):
        testRK_Dipol_SL(2.8,phi/1000.0,3.0)
        testRK_Dipol_SL(2.6,phi/1000.0,3.0)
  #      testRK_Dipol_SL(0.4,0.8,3.0)
    #test_Dipol_VF()
    NeHeGL.main()

if __name__ == "__main__":
    main()
