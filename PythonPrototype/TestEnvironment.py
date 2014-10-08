# -*- coding: utf-8 -*-
"""
Created on Sat Sep 06 11:22:36 2014

@author: Patrick
"""
from sphericalharmo import *
import datastructure2 as ds

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
                tpr = ds.Point3D(theta/1000.0, phi/1000.0, radi/100.0)
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

def testRK_Dipol_SL(theta,phi,r,direction,tmax=1.0e+10,t0=0.8e-3,max_steps=10000):
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    #data = ds.AvsUcdAscii()
    #data.loadFile('E:/Uni/GeodynamicsProject/Datasets/out.1550.inp')
    sl=[]
    vl=[]
    t=t0
    tpr0 = ds.Point3D(theta, phi,r)
    xyz0 = toCartesian(tpr0)
    v0 = evalSHA(xyz0, None)
    sl.append(xyz0)
    vl.append(v0)
    nextVal=v0
    nextPos=xyz0
    print("Start new Streamline with:", tpr0,toSphericalVecfield(tpr0,v0),direction)
    #print(tpr0,toSphericalVecfield(tpr0,v0))
    i= 0 
    while (max_steps>i) and (t<tmax):
       # print(((t-t0)*100.0)/(tmax-t0),"% finished ..... ")
        if(((toSpherical(nextPos)._z)<1.538461538462E+0)):
            print("inner core reached",toSpherical(nextPos))
            #xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,data.getValue,t,t0,direction)
            break
        #print("next pos spherical :",toSpherical(nextPos));
        else:
            xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,evalSHA,t,t0,direction)
        """adapt stepsize"""
        needAdapt= adaptStep(xf,vf,xf2,vf2,t0)[0]
        while needAdapt:
            needAdapt, t0 = adaptStep(xf,vf,xf2,vf2,t0)
            if(((toSpherical(nextPos)._z)<1.538461538462E+0)):
                xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,data.getValue,t,t0,direction)
            #break
        #print("next pos spherical :",toSpherical(nextPos));
            else:
                xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,evalSHA,t,t0,direction)
        nextPos = xf
        nextVal = vf
        t=t2
        sl.append(xf)
        vl.append(vf)
        #print("step #",i,toSpherical(xf), vf,t)
        i+=1
    print("i",i,"sl[i-1]",toSpherical(sl[i-1]),"vl[i-1]",toSphericalVecfield(toSpherical(sl[-1]),vl[i-1]))
    Vis.addStreamLine(sl,vl)
    return
    
def testRK_Whole_SL(theta,phi,r,direction,tmax=1.0e+10,t0=0.8e-3,max_steps=40000):
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    data = ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    sl=[]
    vl=[]
    t=t0
    tpr0 = ds.Point3D(theta, phi,r)
    xyz0 = toCartesian(tpr0)
    v0 = evalSHA(xyz0, None)
    sl.append(xyz0)
    vl.append(v0)
    nextVal=v0
    nextPos=xyz0
    max_step=4.0e-5
    print("Start new Streamline with:", tpr0,toSphericalVecfield(tpr0,v0),direction)
    #print(tpr0,toSphericalVecfield(tpr0,v0))
    i= 0 
    while (max_steps>i) and (t<tmax):
        #print("step: ", i)
       # print(((t-t0)*100.0)/(tmax-t0),"% finished ..... ")
        xdt = nextPos.add(nextVal.mult(t0))
        #print("x + v*dt: ", toSpherical(xdt))
        if( ((toSpherical(nextPos)._z)<1.538461538462E+0) and ((toSpherical(nextPos)._z)>0.538461538462E+0)
            and ((toSpherical(xdt)._z)<1.538461538462E+0) and ((toSpherical(xdt)._z)>0.538461538462E+0)):
          #  print("inner core reached",toSpherical(nextPos))
            xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,data.getValue,t,t0,direction)
            dtmax=data._currentCell.gridSize()/2.0
           # print("dtmax: ",dtmax)
          #  print("outer core tracing :  xf[",i,"] ",toSpherical(xf));
        else:
            xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,evalSHA,t,t0,direction)
            dtmax=0.8e-3
        """adapt stepsize"""
        if(t0>dtmax): 
           # print(i,t0,dtmax)            
            t0=dtmax
        else:
            needAdapt= adaptStep(vf,vf2,t0)[0]
            while needAdapt:
                #print("Stepsize adapted")
                needAdapt, t0 = adaptStep(vf,vf2,t0)
                if(t0>dtmax): 
                  #  print(i,t0,dtmax)            
                    t0=dtmax
                    needAdapt = False                    
                if(((toSpherical(nextPos)._z)<1.538461538462E+0)):
                    xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,data.getValue,t,t0,direction)
                #break
            #print("next pos spherical :",toSpherical(nextPos));
                else:
                    xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,evalSHA,t,t0,direction)
        nextPos = xf
        nextVal = vf
        t=t2
        sl.append(xf)
        vl.append(vf)
        if(i%100)==0:print("step #",i,toSpherical(xf), vf,t)
        i+=1
    print("i",i,"sl[i-1]",toSpherical(sl[i-1]),"vl[i-1]",toSphericalVecfield(toSpherical(sl[-1]),vl[i-1]))
    Vis.addStreamLine(sl,vl)
    return
    
def project_on_boundary(x,v,r):
    """assume x,v are in cartesian"""
    xs = toSpherical(x)
    vs = toSphericalVecfield(xs,v)
    a = (xs._z - r)/(vs._z)        
    sp_t= xs._x + a*vs._x
    sp_p= xs._y + a*vs._y
    sp_r=r
    sp = ds.Point3D(sp_t,sp_p,sp_r)
    print("SP with Boundary in Spherical Coords: ", sp,vs)    
    return toCartesian(sp)

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
    tpr0 = ds.Point3D(0.75*math.pi, 0.5*math.pi,2.91)
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
            tpr = ds.Point3D(theta/1000.0, phi/1000.0, 3.0)
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

def testPerfectDipol(theta, phi,r,direction,tmax=1.0e+10,t0=0.8e-3,max_steps=10000):
    sl=[]
    vl=[]
    t=t0
    tpr0 = ds.Point3D(theta, phi,r)
    xyz0 = toCartesian(tpr0)
    #xyz0 = xyz0.mult(10.0)
    v0 = perfectDipol(xyz0, None)
    sl.append(xyz0)
    vl.append(v0)
    nextVal=v0
    nextPos=xyz0
    print("Start new Streamline with:", tpr0,toSphericalVecfield(tpr0,v0),direction)
    #print(tpr0,toSphericalVecfield(tpr0,v0))
    i= 0 
    while (max_steps>i) and (t<tmax):
        if((toSpherical(nextPos)._z)<2.85):
            print("inner core reached",toSpherical(nextPos))
            break
        xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,perfectDipol,t,t0,direction)
        """adapt stepsize"""
        needAdapt= adaptStep(xf,vf,xf2,vf2,t0)[0]
        while needAdapt:
            needAdapt, t0 = adaptStep(xf,vf,xf2,vf2,t0)
            xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,perfectDipol,t,t0,direction)
        nextPos = xf
        nextVal = vf
        t=t2
        sl.append(xf)
        vl.append(vf)
        #print("step #",i,toSpherical(xf), vf,t)
        i+=1
    print(len(sl))
    print("i",i,"sl[i]",toSpherical(sl[i]),"vl[i]",toSphericalVecfield(toSpherical(sl[i]),vl[i]))
    Vis.addStreamLine(sl,vl)
    return toSpherical(sl[i])
    
def perfectDipol(x,dt):
    xs = toSpherical(x)
    m_dipol= 1.0 ##magnetic dipol moment, only got the r component, as it is (0,0,1)
    m = m_dipol*x._z ## scalar produkt of dipol moment and position
    mu_dipol=1.0##magnetic field constant
    c = mu_dipol/(4*math.pi*(xs._z**3))
    force = 1.8
    x_dipol = c*(3.0*x._x*m) 
    y_dipol = c*(3.0*x._y*m)
    z_dipol = c*(3.0*x._z*m - force*xs._z**2)
    vs = ds.Point3D(x_dipol,y_dipol,z_dipol)
    return vs
    
def main():
    """Test of Extrapolation Method """
    #for phi in range(10,2*3141,500):
       #testRK_Dipol_SL(1.4,phi/1000.0,2.8,"forward")
       # testRK_Dipol_SL(0.3,phi/1000.0,2.3,"forward")
       # testRK_Dipol_SL(0.6,phi/1000.0,2.3,"forward")
    """Test of whole Streamline Vis"""    
    testRK_Whole_SL(0.3,0.1,0.9,"forward")
    """test of Integration Method"""
    #endpoint = testPerfectDipol(0.2,10/1000.0,2.9,"forward")
    #endpoint1 = testPerfectDipol(endpoint._x,endpoint._y,endpoint._z,"backward")
    #err=endpoint._x - endpoint1._x + endpoint._y - endpoint1._y +endpoint._z - endpoint1._z
    print("Error for RK4: ", err)
    """end of Test"""
    NeHeGL.main()

if __name__ == "__main__":
    main()
