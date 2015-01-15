# -*- coding: utf-8 -*-
"""
Created on Sat Sep 06 11:22:36 2014

@author: Patrick
"""
from sphericalharmo import *
import sphericalharmo as sph
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
    dtmax = 10000.0
    while (max_steps>i) and (t<tmax):
       # print(((t-t0)*100.0)/(tmax-t0),"% finished ..... ")
        if(((toSpherical(nextPos)._z)<1.6)):
            print("inner core reached",toSpherical(nextPos))
            #xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,data.getValue,t,t0,direction)
            break
        #print("next pos spherical :",toSpherical(nextPos))
        else:
            xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,evalSHA,t,t0,direction)
        """adapt stepsize"""
        needAdapt= adaptStep(vf,vf2,t0)[0]
        while needAdapt:
            needAdapt, t0 = adaptStep(vf,vf2,t0)
            #print("timestep after adapt",t0)
            if(t0>dtmax):         
                t0=dtmax
                needAdapt = False
        #print("next pos spherical :",toSpherical(nextPos));
            else:
                xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,evalSHA,t,t0,direction)
        nextPos = xf
        nextVal = vf
        t=t2
        sl.append(xf)
        vl.append(vf)
        #print("step #",i,toSpherical(xf), vf,t)
        if(i%10)==0:print("step #",i,toSpherical(xf), vf,t0)
        i+=1
    print("i",i,"sl[i-1]",toSpherical(sl[i-1]),"vl[i-1]",toSphericalVecfield(toSpherical(sl[-1]),vl[i-1]))
    Vis.addStreamLine(sl,vl)
    return

def getSmallestGridSize():
    data = ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    data.builtKDTree()
    l_grid=[]
    for ver in data._vertexList:
        xtupple = (ver._pos._x,ver._pos._y,ver._pos._z)
        d,ni = data._kdTree.query(xtupple,1) ## return the indices of the nearest neigbhour, d and ni are arrays
        l_grid.append(d)
    print(min(l_grid))
   
def test_OC_only(theta,phi,r,direction,data,tmax=1.0e+10,t0=1.0e-10,max_steps=5000):
    sph.setData(data)
    sl=[]
    vl=[]
    t=t0
    """dtmin is dependent on the OS (32 or 64 bit)"""
    dtmin=1.0e-16
    dtmax = 1.0e-10
    tpr0 = ds.Point3D(theta, phi,r)
    xyz0 = toCartesian(tpr0)
    v0 = evalSHA(xyz0, None)
    dtmax = data._CurrentMaxStep
    #print("dtmax_ first:",dtmax)
    sl.append(xyz0)
    vl.append(v0)
    nextVal=v0
    nextPos=xyz0
    print("Start new Streamline with:", tpr0,toSphericalVecfield(tpr0,v0),direction)
    #print(tpr0,toSphericalVecfield(tpr0,v0))
    i= 0
    outOfBounds = False
    while (max_steps>i) and (t<tmax) and  (not outOfBounds):
        xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,evalSHA,t,t0,direction)
        if(((toSpherical(xf)._z)>=1.537983852128)):
            outOfBounds = True
            break
       
        dtmax = data._CurrentMaxStep
        print("dtmax:",dtmax)
        if(t0<dtmin):           
            t0=dtmin
        if(t0>dtmax):           
            t0=dtmax
        needAdapt= adaptStepGearHairer(vf,vf2,t0)[0]
        print("Before while", needAdapt)
        while needAdapt:
            needAdapt, t0 = adaptStepGearHairer(vf,vf2,t0)
            print("inside while", needAdapt)
            if(t0>dtmax):         
                t0=dtmax
                needAdapt = False
                print("stepsize adapted 2") 
            if(t0<dtmin):           
                t0=dtmin
                needAdapt = False
                print("stepsize adapted 3")                    
            xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,evalSHA,t,t0,direction)
            dtmax = data._CurrentMaxStep
            if(((toSpherical(xf)._z)>=1.537983852128)):
                outOfBounds = True
                break
        nextPos = xf
        nextVal = vf
        t=t2
        sl.append(xf)
        vl.append(vf)
        """if(i%100)==0:"""
        print("--------- step #",i,toSpherical(xf), vf,t,"-------------")
        i+=1
    #print("i",i,"sl[i-1]",toSpherical(sl[i-1]),"vl[i-1]",toSphericalVecfield(toSpherical(sl[-1]),vl[i-1]))
    Vis.addStreamLine(sl,vl)
    return
    
def loadData(path):
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    data = ds.VTKData()
    data.loadFile(path)
    return data
    
def testRK_Whole_SL(theta,phi,r,direction,tmax=1.0e+10,t0=0.8e-4,max_steps=1000):
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    data = ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    data.builtKDTree()
    sph.setData(data)
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
        xf,vf,t2,xf2,vf2,tf2 = rk4(nextPos,nextVal,evalSHA,t,t0,direction)
        if(((toSpherical(nextPos)._z)<1.537983852128+0)):
            dtmax=data._currentCell.gridSize()/2.0
        else:
            dtmax=1.0
        """adapt stepsize"""
        if(t0>dtmax):           
            t0=dtmax
        needAdapt= adaptStep(vf,vf2,t0)[0]
        while needAdapt:
            needAdapt, t0 = adaptStep(vf,vf2,t0)
            if(t0>dtmax):         
                t0=dtmax
                needAdapt = False                    
            xf,vf,t2,xf2,vf2,tf2 = rk4(xf,vf,evalSHA,t,t0,direction)
        nextPos = xf
        nextVal = vf
        t=t2
        sl.append(xf)
        vl.append(vf)
        if(i%100)==0:print("step #",i,toSpherical(xf), vf,t)
        i+=1
    #print("i",i,"sl[i-1]",toSpherical(sl[i-1]),"vl[i-1]",toSphericalVecfield(toSpherical(sl[-1]),vl[i-1]))
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
    
def testBoundaryVecField():
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    data = ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    data.builtKDTree()
    sph.setData(data)
    
    xf_mantle=[]
    xf_oc = []
    for theta in range(10,2*3141,700):
        for phi in range(100,3141,350):
            tpr_m = ds.Point3D(theta/1000.0, phi/1000.0, 1.537983852128 + 1.0e-4)
            tpr_oc = ds.Point3D(theta/1000.0, phi/1000.0, 1.537983852128 - 1.0e-3)
            xyz_m = ds.toCartesian(tpr_m)
            xyz_oc = ds.toCartesian(tpr_oc)
            xf_mantle.append(xyz_m)
            xf_oc.append(xyz_oc)
    vf_mantle = []
    vf_oc = []
    i=0
    for x_oc in xf_oc:
        i+=1
        print(float(i)*100/len(xf_mantle))
        #print("OC>>>Pos:",toSpherical(x_oc) )
        print("input pos: ", ds.toSpherical(x_oc))
        v = data.getValueKDTree(x_oc,1.0)
        if((i*100.0/len(xf_oc))%10)==0: print("Evaluation reached: " +str(i*100.0/len(xf_oc)) +"%")
       # print("OC>>>Pos:",x_oc ,"Value:",v)
        vf_oc.append(v)
    j=0    
    for x_m in xf_mantle:
        j+=1
        print(float(j)*100/len(xf_mantle))
        #print("MANTEL>>>Pos:",toSpherical(x_m) )
        v = evalSHA(x_m,1.0)
        #print("MANTEL>>>Pos:",x_m ,"Value:",v)
        if((j*100.0/len(xf_mantle))%10)==0: 
            print("Evaluation reached: " +str(j*100.0/len(xf_mantle)) +"%")
        vf_mantle.append(v)
        
   # xf = xf_mantle + xf_oc
   # vf = vf_mantle + vf_oc    
    Vis.setVectorfield(xf_mantle,vf_mantle)  #purple
    Vis.setVectorfield2(xf_oc,vf_oc) # light blue
    return    

def testBoundaryVecField2():
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    data = ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    data.builtKDTree()
    sph.setData(data)
    vf_mantle = []
    vf_oc = []
    xf_mantle=[]
    xf_oc = []
    for theta in range(10,2*3141,315):
        for phi in range(10,3141,158):
            ##choose verteces out of the data set -> avoid interpolation
            tpr_oc = ds.Point3D(theta/1000.0, phi/1000.0, 1.537983852128 - 1.0e-4)
            xyz_oc = toCartesian(tpr_oc)
            xyz_tripple = (xyz_oc._x,xyz_oc._y,xyz_oc._z)
            di,oc_index = data._kdTree.query(xyz_tripple)
            xyz_oc = data._vertexList[oc_index]._pos
            vf_oc.append(data._vertexList[oc_index]._mag)
            xf_oc.append(xyz_oc)
            
            tpr_m= toSpherical(xyz_oc).add(Point3D(0.0,0.0,1.0e-4))
            xyz_m = toCartesian(tpr_m)
            xf_mantle.append(xyz_m)
            vf_mantle.append(sphericalHarmoAnalysis(xyz_m))
    print("finised computation")        
    #Vis.setVectorfield(xf_mantle,vf_mantle)
    Vis.setVectorfield2(xf_oc,vf_oc)
    return    

def DS_compared_Vecfield():
    loadGaussCoefSimu("../../Gauss_RE.dat","../../Gauss_ICB.dat")
    data = ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    data.builtKDTree()
    sph.setData(data)
    
    xf_oc = []
    for theta in range(10,2*3141,315):
        for phi in range(10,3141,158):
            tpr_oc = ds.Point3D(theta/1000.0, phi/1000.0, 1.537983852128 - 1.0e-3)
            xyz_oc = toCartesian(tpr_oc)
            xf_oc.append(xyz_oc)
    vf_oc_kd = []
    vf_oc_bf = []
    i=0
    for x_oc in xf_oc:
        i+=1
        #print("OC>>>Pos:",toSpherical(x_oc) )
        v_bf = data.getValue(x_oc,1.0)
        v_kd = data.getValueKDTree(x_oc,1.0)
        
        if((i*100.0/len(xf_oc))%10)==0: print("Evaluation reached: " +str(i*100.0/len(xf_oc)) +"%")
       # print("OC>>>Pos:",x_oc ,"Value:",v)
        vf_oc_bf.append(v_bf)
        vf_oc_kd.append(v_kd)
    #xf = xf_mantle + xf_oc
    #vf = vf_mantle + vf_oc    
    Vis.setVectorfield(xf_oc,vf_oc_bf)    #purple
    Vis.setVectorfield2(xf_oc,vf_oc_kd)   #light blue
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
 
def GridTestVis(num_Cells):
    data = ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    #Vis.initVertexBuffer(data._vertexList)
    cellList = []
    for i in range(1,len(data._cellList),len(data._cellList)/num_Cells):
        cellList.append(data._cellList[i])
    Vis.setCellList(cellList)


def VisGridConcave():
    data=ds.VTKData()
    data.loadFile('C:/out.1200.vtk')
    data.computeCellTopology()
    cellList = computeCellVolumes(data)
    visualCellList =[]
    for cell in cellList:
        cell.computeFaceNeighbours(data)
        if len(set(cell._verts))<8:
            print(set(cell._verts), cell._ID)
        #for nb in cell._faceNeighbours:
            #if nb != None:
                #visualCellList.append(data._cellList[nb-1])
    #Vis.setCellList(visualCellList)
    print("Concave test finished")
    return
        
def computeCellVolumes(data):
    vol=0.0
    vol_List =[]
    cell_List = []
    ##to avoid misunderstanding in Cell ID numberin and list numbering, Cell ID 0 doesn't exist
    #vol_List.append(0.0)
    print("compute Cell volumes ... ")
    for cell in data._cellList:
        v=[]
        for vert in cell._verts:
            v.append(vert._pos)
        def pV(a,b):
            return v[b].sub(v[a])
        vol =1.0/12.0* ds.dot(pV(2,4),(ds.cross(pV(7,5),pV(4,6) )).add(ds.cross(pV(5,0),pV(4,1)).add(ds.cross(pV(0,7),pV(4,3) ) ) ) )
        vol+=1.0/12.0* ds.dot(pV(7,2),(ds.cross(pV(5,2),pV(1,6) )).add(ds.cross(pV(7,6),pV(3,6))))
        vol+=1.0/12.0* ds.dot(pV(7,0),(ds.cross(pV(2,0),pV(1,3) )).add(ds.cross(pV(7,3),pV(4,3))))
        vol+=1.0/12.0* ds.dot(pV(7,5),(ds.cross(pV(0,5),pV(1,4) )).add(ds.cross(pV(7,4),pV(6,4))))
        if( vol < 0.0):
            verts=[]
            for ver in cell._verts: 
                verts.append(toSpherical(ver._pos))
            print(cell._ID,vol,verts)
            vol_List.append(vol)
            cell_List.append(cell) 
        #else:    
            #vol_List.append(vol)
    maxVol = max(vol_List)
    minVol = min(vol_List)
    print(" .. finished cell volumes")
    print("Max Volume: ", maxVol,"Min Volume: ", minVol)
    return cell_List        

def main():
    data= ds.VTKData()
    data = loadData('C:/out.1200.vtk')
    #Vis.built_cm_rainbow(data)
    data.builtKDTree()
    #VisGridConcave()
    #GridTestVis(5000)                              ##define number of Cells to be visible
    #DS_compared_Vecfield()
    #testBoundaryVecField()
    #testRK_Dipol_SL(2.8,0.2,2.8,"forward")
    """Test of Extrapolation Method """
    #for phi in range(10,2*3141,500):
       #testRK_Dipol_SL(2.8,phi/1000.0,2.8,"forward")
       # testRK_Dipol_SL(0.3,phi/1000.0,2.3,"forward")
     #  testRK_Dipol_SL(0.6,phi/1000.0,2.3,"forward")
    """Test of whole Streamline Vis""" 
    #for phi in range(200,2*3141,600):
     #   test_OC_only(1.4,phi/1000.0,0.7,"forward",data)
      #  test_OC_only(1.4,phi/1000.0,0.7,"backward",data)
    test_OC_only(1.4,0.6,0.7,"forward",data)
    #test_OC_only(1.0,0.6,1.0,"forward",data)
    ##----degenereted case
    #testRK_Whole_SL(1.4,0.6,1.1,"forward")
    ##
    #testRK_Whole_SL(0.8,0.6,1.1,"forward")
    """test of Integration Method"""
    #endpoint = testPerfectDipol(0.2,10/1000.0,2.9,"forward")
    #endpoint1 = testPerfectDipol(endpoint._x,endpoint._y,endpoint._z,"backward")
    #err=endpoint._x - endpoint1._x + endpoint._y - endpoint1._y +endpoint._z - endpoint1._z
    #print("Error for RK4: ", err)
    """end of Test"""
    NeHeGL.main()

if __name__ == "__main__":
    main()
