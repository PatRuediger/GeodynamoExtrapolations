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

def testRK_Dipol_SL():
    loadGaussCoefIGRF("../GausCoef.txt")
    useIGRFonly() #overrites all other Coefficients
    sl=[]
    vl=[]
    tmax = 1000.0
    #initial stepsize, optional
    step = 0.8e-2
    max_steps = 100
    t=step
    tpr0 = datastructure.Point3D(0.5, 0.5,2.92)
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
        if(toSpherical(nextPos)._z<2.91):
            print(toSpherical(nextPos))
            print("Core reached")
            break
        #print(toSpherical(nextPos));
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
        print("step #",i,xf,vf)
        i+=1     
    
    #vf = toSphericalVecfield(tpr,v)
    #print(tpr,vf,v)
    Vis.setStreamLine(sl,vl)
    NeHeGL.main()
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
        print("step #",i,xf,vf)
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
    
def main():
    """only dummy values here """
    #testSHA_SL(1.0879, 0.0, -1.0879, 0,0,0)
    #testSHAwithVecfield()
    #testRK_Dipol_SL()
    test_Dipol_VF()

if __name__ == "__main__":
    main()
