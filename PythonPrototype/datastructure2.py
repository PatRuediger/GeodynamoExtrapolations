# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 21:31:00 2014

@author: Patrick
"""

#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon May 12 10:46:55 2014

@author: Patrick

Data structure for Visualization
Using hexahedron Grid cells with 8 vertices
"""
import math
import numpy as np
from numpy.linalg import inv
from scipy.spatial import KDTree
import scipy.optimize as opt
from datetime import datetime
import timeit

def mainTest():
    print('Main startet')
    DS = VTKData()
    DS.loadFile('C:/out.1200.vtk')
    DS.builtKDTree()
    x = Point3D(0.5,0.4,0.75)
    tstart = datetime.now()
    bf = DS.getValue(x,1.0)
    tend = datetime.now()
    delta = tend -tstart
    #print("Search time for Brute Force",  dt)
    print(bf)
    print("Search time for Brute Force",  int(delta.total_seconds() * 1000))
    tstart = datetime.now()
    kdv = DS.getValueKDTree(x,1.0)
    tend = datetime.now()
    delta = tend -tstart
    print(kdv)
    print("Search time for KD Tree", int(delta.total_seconds() * 1000))



def testInterpolation():
     #           c000,c100,c110,c010,c001,c101,c111,c011
    v0 = Vertex(0,Point3D(0.0,0.0,0.0),Point3D(1,0,0))
    v1 = Vertex(1,Point3D(1.0,0.0,0.0),Point3D(1,1,0))
    v2 = Vertex(2,Point3D(1.0,1.0,0.0),Point3D(1,0,0))
    v3 = Vertex(3,Point3D(0.0,1.0,0.0),Point3D(1,0,0))
    v4 = Vertex(4,Point3D(0.0,0.0,1.0),Point3D(1,0,1))
    v5 = Vertex(5,Point3D(1.0,0.0,1.0),Point3D(1,0,0))
    v6 = Vertex(6,Point3D(1.0,2.0,1.0),Point3D(1,1,0))
    v7 = Vertex(7,Point3D(0.0,1.0,1.0),Point3D(1,0,0))

    cell1 = Cell(0,v0,v1,v2,v3,v4,v5,v6,v7)
    

    p= Point3D(0.6,0.5,0.5)
    output = cell1.isInside(p)
    print(output)
    #temp = cell1.trilinear(p)
    #print(temp)

"""Dot/Scalar Product for 3D Vectors"""
def dot(a,b):
        a0 = Point3D(0,0,0)
        b0 = Point3D(0,0,0)
        a0 = a
        b0 = b
        return a0._x * b0._x + a0._y * b0._y + a0._z * b0._z

def cross(a,b):
    c = Point3D(0.0,0.0,0.0)
    c._x = a._y*b._z - a._z*b._y
    c._y = a._z*b._x - a._x*b._z
    c._z = a._x*b._y - a._y*b._x
    return c

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

class Point3D:
#    _x = 0.0
#    _y = 0.0
#    _z = 0.0

    def _length(self):
        return math.sqrt(self._x**2.0+self._y**2.0+self._z**2.0)

    def toNumpyArray(self):    
        return np.array((self._x,self._y,self._z))
        
    def __init__(self,x,y,z):
        self._x = x
        self._y = y
        self._z = z
        
    def __getitem__(self,i):
        if (i==0): return self._x
        elif(i==1): return self._y
        elif(i==2): return self._z
        else: print("Point got only 3-Dimensions, use _length() ")
    
    def scalarProduct(self,x):
        return self._x*x._x+self._y*x._y +self._z*x._z
         
    def mult(self,scalar):
        return Point3D(self._x * scalar,self._y * scalar,self._z * scalar)

    def add(self,vec):
        return Point3D(self._x + vec._x, self._y + vec._y, self._z + vec._z)
        
    def sub(self,vec):
        return Point3D(self._x - vec._x, self._y - vec._y, self._z - vec._z)   
    
    def __repr__(self):
        return "("+str(self._x)+","+str(self._y)+","+str(self._z)+")"


class Vertex:
#---class Attributes should only be readed
    #Vertex ID
    _ID = 0
    #Postion
    _pos = Point3D(0,0,0)

    #Magnetic field
    _mag = Point3D(0,0,0)
    #Pointer to nex Vertex
    _nextVertex = None
   
    _partOfCell=[]


#---class Methods
    # Constructor method
    def __init__(self,ID,pos,mag):
        self._ID = ID
        self._pos = pos
        self._mag = mag
        self._partOfCell = []

    def addVertex(self,ID,pos,mag):
        newVertex = Vertex(ID,pos,mag)
        self._nextVertex = newVertex
        return self._nextVertex



#---- Cell is ment to be an hexahedron
class Cell:
#---class Attributes should only be readed
    _ID = 0
    # should be 8 everytime
    # ordering:c000,c100,c110,c010,c001,c101,c111,c011
    # leads to 0,1,2,3,4,5,6,7 in the array
    _verts =[]
    #Pointer to next Cell
    _nextCell = None
    
    #ordering: front,back,left,right,top,bottom
    #Face sharing neighbours
    _faceNeighbours =[None]*6
    _faceNormals = [None]*6
    #the order is not given
    #To access the cells in the list, you must use ID-1
    _neighbours=[]

#---class Methods
    # Constructor method
    def __init__(self,ID,v0,v1,v2,v3,v4,v5,v6,v7):
        self._verts = [v0,v1,v2,v3,v4,v5,v6,v7]
        self._ID = ID
       # self.computeCellNormals()

    def computeFaceNeighbours(self,data):
        front =[]
        back = []
        left = []
        right = []
        top = []
        bottom = []
        for nb in self._neighbours:
            face = set(data._cellList[nb-1]._verts).intersection(set([self._verts[0],self._verts[1],self._verts[2],self._verts[3]]))
            if len(face) ==4:
                bottom.append(nb)
            face = set(data._cellList[nb-1]._verts).intersection(set([self._verts[4],self._verts[5],self._verts[6],self._verts[7]]))
            if len(face) ==4:
                top.append(nb)
            face = set(data._cellList[nb-1]._verts).intersection(set([self._verts[0],self._verts[1],self._verts[4],self._verts[5]]))
            if len(face) ==4:
                left.append(nb)
            face = set(data._cellList[nb-1]._verts).intersection(set([self._verts[2],self._verts[3],self._verts[6],self._verts[7]]))
            if len(face) ==4:
                right.append(nb)
            face = set(data._cellList[nb-1]._verts).intersection(set([self._verts[1],self._verts[2],self._verts[5],self._verts[6]]))
            if len(face) ==4:
                front.append(nb)
            face = set(data._cellList[nb-1]._verts).intersection(set([self._verts[0],self._verts[3],self._verts[4],self._verts[7]]))
            if len(face) ==4:
                back.append(nb)
        front.remove(self._ID)        
        back.remove(self._ID)
        left.remove(self._ID)
        right.remove(self._ID)
        top.remove(self._ID)
        bottom.remove(self._ID)
            
        if len(front) >1 : print("multiple Front Face Neighbours for Cell", self._ID)        
        if len(back) >1 : print("multiple back Face Neighbours for Cell", self._ID)
        if len(left) >1 : print("multiple left Face Neighbours for Cell", self._ID)
        if len(right) >1 : print("multiple right Face Neighbours for Cell", self._ID)
        if len(top) >1 : print("multiple top Face Neighbours for Cell", self._ID)
        if len(bottom) >1 : print("multiple bottom Face Neighbours for Cell", self._ID)   
        #print(front,back,left,right,top,bottom)    
        if len(front) !=0:
            self._faceNeighbours.append(data._cellList[front[0]-1])
        if len(back) !=0:
            self._faceNeighbours.append(data._cellList[back[0]-1])
        if len(left) !=0:
            self._faceNeighbours.append(data._cellList[left[0]-1])
        if len(right) !=0:
            self._faceNeighbours.append(data._cellList[right[0]-1])
        if len(top) !=0:
            self._faceNeighbours.append(data._cellList[top[0]-1])
        if len(bottom) !=0:
            self._faceNeighbours.append(data._cellList[bottom[0]-1])
    
            
    def addCell(self,ID,v0,v1,v2,v3,v4,v5,v6,v7):
        newCell = Cell(ID,v0,v1,v2,v3,v4,v5,v6,v7)
        self._nextCell = newCell
        return self._nextCell


    ## out of date method!!!
    def computeCellNormals(self):
        """All normals are oriented to point inwards
        """
        self._faceNormals[0] = cross(self._verts[0]._pos.sub(self._verts[3]._pos) ,self._verts[2]._pos.sub(self._verts[3]._pos) ) 
        self._faceNormals[0] = self._faceNormals[0].mult(1.0/self._faceNormals[0]._length())
        
        self._faceNormals[1] = cross(self._verts[6]._pos.sub(self._verts[7]._pos) ,self._verts[4]._pos.sub(self._verts[7]._pos) ) 
        self._faceNormals[1] = self._faceNormals[1].mult(1.0/self._faceNormals[1]._length())
        
        self._faceNormals[2] = cross(self._verts[0]._pos.sub(self._verts[1]._pos) ,self._verts[5]._pos.sub(self._verts[1]._pos) ) 
        self._faceNormals[2] = self._faceNormals[2].mult(1.0/self._faceNormals[2]._length())    

        self._faceNormals[3] = cross(self._verts[2]._pos.sub(self._verts[3]._pos) ,self._verts[7]._pos.sub(self._verts[3]._pos) ) 
        self._faceNormals[3] = self._faceNormals[3].mult(1.0/self._faceNormals[3]._length())

        self._faceNormals[4] = cross(self._verts[1]._pos.sub(self._verts[2]._pos) ,self._verts[6]._pos.sub(self._verts[2]._pos) ) 
        self._faceNormals[4] = self._faceNormals[4].mult(1.0/self._faceNormals[4]._length())

        self._faceNormals[5] = cross(self._verts[3]._pos.sub(self._verts[0]._pos) ,self._verts[4]._pos.sub(self._verts[0]._pos) ) 
        self._faceNormals[5] = self._faceNormals[5].mult(1.0/self._faceNormals[5]._length())
        #print(self._faceNormals)
        return
       
    def inversejacobiTrilinear(self,x,a,b,c):
        """input trilinear interpolated point x at alpha beta gamma"""
    
        """J_a = self._verts[0]._mag.mult((1.0-b)*(1.0-c))
        J_a = J_a.add(self._verts[1]._mag.mult(((1.0-b)*(1.0-c)) ))
        J_a = J_a.add(self._verts[3]._mag.mult(b*(1.0-c)))
        J_a =J_a.add(self._verts[4]._mag.mult((1.0-b)*c))
        J_a =J_a.add(self._verts[5]._mag.mult((1.0-b)*c))
        J_a = J_a.add(self._verts[7]._mag.mult(b*c))
        J_a =J_a.add(self._verts[2]._mag.mult(b*(1.0-c)))
        J_a = J_a.add(self._verts[6]._mag.mult(b*c))      
        
        J_b = self._verts[0]._mag.mult((1.0-a)*(1.0-c))
        J_b = J_b.add(self._verts[1]._mag.mult((a*(1.0-c)) ))
        J_b = J_b.add(self._verts[3]._mag.mult((1.0-a)*(1.0-c)))
        J_b =J_b.add(self._verts[4]._mag.mult((1.0-a)*c))
        J_b =J_b.add(self._verts[5]._mag.mult(a*c))
        J_b = J_b.add(self._verts[7]._mag.mult((1.0-a)*c))
        J_b =J_b.add(self._verts[2]._mag.mult(a*(1.0-c)))
        J_b = J_b.add(self._verts[6]._mag.mult(a*c)) 
        
        J_c = self._verts[0]._mag.mult((1.0-a)*(1.0-b))
        J_c = J_c.add(self._verts[1]._mag.mult((a*(1.0-b)) ))
        J_c = J_c.add(self._verts[3]._mag.mult((1.0-a)*b))
        J_c = J_c.add(self._verts[4]._mag.mult((1.0-a)*(1.0-b)))
        J_c = J_c.add(self._verts[5]._mag.mult(a*(1.0-b)))
        J_c = J_c.add(self._verts[7]._mag.mult((1.0-a)*b))
        J_c = J_c.add(self._verts[2]._mag.mult(a*b))
        J_c = J_c.add(self._verts[6]._mag.mult(a*b))
                
        jacobi = np.array(((J_a[0],J_b[0],J_c[0]),(J_a[1],J_b[1],J_c[1]),(J_a[2],J_b[2],J_c[2])))"""
        
       # J0=Point3D(a +1.0,b,c).sub(x)  
      #  J1=Point3D(a,b +1.0,c).sub(x)  
       # J2=Point3D(a,b,c +1.0).sub(x)  
        
        J0=Point3D(a +1.0,b,c).sub(Point3D(a -1.0,b,c))
        J0= J0.mult(0.5)
        J1=Point3D(a,b +1.0,c).sub(Point3D(a,b -1.0,c)) 
        J1= J1.mult(0.5)
        J2=Point3D(a,b,c +1.0).sub(Point3D(a,b,c -1.0))
        J2= J2.mult(0.5)
        
        jacobi = np.array(((J0[0],J1[0],J2[0]),(J0[1],J1[1],J2[1]),(J0[2],J1[2],J2[2])))
        #print("Jacobi:",jacobi)
       # print("det ",np.linalg.det(jacobi))
        return inv(jacobi)
    
    
    ## @input is a Point3D 
    def isInside(self,p):
        alpha = 0.5
        beta = 0.5
        gamma = 0.5
        delta = 1.0e-9  ##error threshold
        max_iter = 100
        delta_xc = Point3D(1.0,1.0,1.0)
        i=0
        while(delta_xc._length()>delta):
            i+=1
            xp_star = self.trilinear_pos(alpha,beta,gamma)
           # print("xp_star: ",xp_star)
            delta_xp = xp_star.sub(p)
          #  print("delta_xp_star: ",delta_xp)
            delta_xc_array = self.inversejacobiTrilinear(xp_star,alpha,beta,gamma).dot(delta_xp.toNumpyArray())
          #  print ("delta_xc_array: ",delta_xc_array)
            alpha += delta_xc_array[0]
            beta += delta_xc_array[1]    
            gamma += delta_xc_array[2]
            if alpha >1.0 or beta >1.0 or gamma >1.0: 
                return False, None
            delta_xc = Point3D(delta_xc_array[0],delta_xc_array[1] ,delta_xc_array[2])
            if(delta_xc._length()<delta or i >= max_iter):
                return True,(alpha,beta,gamma)
    
    def gridSize(self):
        bounds = self.boundaries()
        gridlen = []
        gridlen.append(bounds[3]-bounds[0])
        gridlen.append(bounds[4]-bounds[1])
        gridlen.append(bounds[5]-bounds[2])
        return min(gridlen)

    def trilinear_pos(self,a,b,c):
         # ordering:c000,c100,c110,c010,c001,c101,c111,c011
        result = self._verts[0]._pos.mult((1.0-a)*(1.0-b)*(1.0-c))
        result = result.add(self._verts[1]._pos.mult((a*(1.0-b)*(1.0-c)) ))
        result = result.add(self._verts[3]._pos.mult((1.0-a)*b*(1.0-c)))
        result = result.add(self._verts[4]._pos.mult((1.0-a)*(1.0-b)*c))
        result = result.add(self._verts[5]._pos.mult(a*(1.0-b)*c))
        result = result.add(self._verts[7]._pos.mult((1.0-a)*b*c))
        result =result.add(self._verts[2]._pos.mult(a*b*(1.0-c)))
        result = result.add(self._verts[6]._pos.mult(a*b*c))        
        return result
            
    def trilinear_mag(self,a,b,c):
        result = self._verts[0]._mag.mult((1.0-a)*(1.0-b)*(1.0-c))
        result = result.add(self._verts[1]._mag.mult((a*(1.0-b)*(1.0-c)) ))
        result = result.add(self._verts[3]._mag.mult((1.0-a)*b*(1.0-c)))
        result = result.add(self._verts[4]._mag.mult((1.0-a)*(1.0-b)*c))
        result = result.add(self._verts[5]._mag.mult(a*(1.0-b)*c))
        result = result.add(self._verts[7]._mag.mult((1.0-a)*b*c))
        result = result.add(self._verts[2]._mag.mult(a*b*(1.0-c)))
        result = result.add(self._verts[6]._mag.mult(a*b*c))         
        return result
        
    ## @input is a Point3D 
    ## @output is a Point3D
    def trilinear(self,x):
        if(self.isInside(x)[0]):
            a = self.isInside(x)[1][0]
            b = self.isInside(x)[1][1]
            c = self.isInside(x)[1][2]
            result = self._verts[0]._mag.mult((1.0-a)*(1.0-b)*(1.0-c)).add(self._verts[1]._mag.mult((a*(1.0-b)*(1.0-c)) ))
            result.add(self._verts[3]._mag.mult((1.0-a)*b*(1.0-c)).add(self._verts[4]._mag.mult((1.0-a)*(1.0-b)*c)))
            result.add(self._verts[5]._mag.mult(a*(1.0-b)*c).add(self._verts[7]._mag.mult((1.0-a)*b*c)))
            result.add(self._verts[2]._mag.mult(a*b*(1.0-c)).add(self._verts[6]._mag.mult(a*b*c)))        
            return True,result,self
        else:
            return False,None,self



class VTKData:
    """class for VTK Data Files"""
    _currentCell=None
    _numVert =0
    _numCells = 0
    _dim = 0
    _vertexList =[]
    _cellList =[]
    _valueNames =[]
    _kdTree = None

#---class Methods
    # Constructor method
    def __init__(self):
          return
           
    
    def builtKDTree(self):
        print("Built KD- Tree ... ")
        vertexListTripples = []
        ## List Index is the same as in _vertexList
        for elem in self._vertexList:
            vertexListTripples.append((elem._pos._x,elem._pos._y,elem._pos._z))
        self._kdTree = KDTree(vertexListTripples)

    def getValueKDTree(self,x,dt):
        xtupple = (x._x,x._y,x._z)
        d,ni = self._kdTree.query(xtupple,50) ## return the indices of the nearest neigbhour, d and ni are arrays
        #print(max(d))
        for ver in ni:
            nnVertex = self._vertexList[ver]
            #print(nnVertex._partOfCell)
            for neighbourID in nnVertex._partOfCell:
                self._currentCell = self._cellList[neighbourID-1]
                isFound,intPoint, nextCell = self._currentCell.trilinear(x)
                if(isFound): 
                    return intPoint
        print("nearest neighbour didn't help",toSpherical(x))
        for cell in self._cellList:
            isFound,intPoint, nextCell = cell.trilinear(x)
            if(isFound):
                print("Cell Found", cell._ID, intPoint.sub(x)._length())
                return intPoint
        print("Cell not Found",x)        
 
    def trilinearInterpolation(self,verts,x):
        ##trilinear interpolation of x with lexikographically ordered input array verts
        x0=verts[0]._pos[0]
        x1=verts[4]._pos[0]
        y0=verts[0]._pos[1]
        y1=verts[3]._pos[1]
        z0=verts[0]._pos[2]
        z1=verts[1]._pos[2]

        if((x1-x0)==0.0): xd=0.0
        else: xd = (x[0]-x0)/(x1-x0)
        if((y1-y0)==0.0): yd =0.0
        else: yd = (x[1]-y0)/(y1-y0)
        if((z1-z0)==0.0): zd=0.0
        else: zd = (x[2]-z0)/(z1-z0)

        c00 =  Point3D(0,0,0)
        c10 =  Point3D(0,0,0)
        c01 =  Point3D(0,0,0)
        c11 =  Point3D(0,0,0)
        c0 =  Point3D(0,0,0)
        c1 =  Point3D(0,0,0)
        c =  Point3D(0,0,0)

        c00 = verts[0]._mag.mult(1.0-xd)
        c00=c00.add(verts[4]._mag.mult(xd))
        c10 = verts[2]._mag.mult(1.0-xd)
        c10=c10.add(verts[6]._mag.mult(xd))
        c01 = verts[1]._mag.mult(1.0-xd)
        c01=c01.add(verts[5]._mag.mult(xd))
        c11 = verts[3]._mag.mult(1.0-xd)
        c11=c11.add(verts[7]._mag.mult(xd))
        
        c0= c00.mult(1.0-yd)
        c0=c0.add(c10.mult(yd))
        c1=c01.mult(1.0-yd)
        c1=c1.add(c11.mult(yd))
       
        c = c0.mult(1.0-zd)
        c=c.add(c1.mult(zd))
        
        return c
        
    def getValueNNInt(self,x,dt):
        xtupple = (x._x,x._y,x._z)
        d,ni = self._kdTree.query(xtupple,8) ## return the indices of the nearest neigbhour, d and ni are arrays
        if ni == []:
            print("Cell not Found",x)
            return
        nn_list = []
        for i in ni:
            nn_list.append(self._vertexList[i])
        nn_list = sorted(nn_list)
        return self.trilinearInterpolation(nn_list,x)
                
        
    ## returns the interpolated value at position x
    def getValue(self,x,dt):
     #   print("Search for Cell with point: ", toSpherical(x))
        isFound,intPoint,self._currentCell = self._currentCell.trilinear(x)
        if(isFound): 
           # print(self._currentCell._ID,isFound,"Found in Current Cell")                    
            return intPoint        
       # print("Searching in neighbours")
        for ver in self._currentCell._verts:
            for neighbourID in ver._partOfCell:
                neighbour = self._cellList[neighbourID-1]
                isFound,intPoint,self._currentCell = neighbour.trilinear(x)
                if(isFound):
                    #print("Celllist of Brute Force",ver._partOfCell)
                    return intPoint
       #if not found in neighbouring cell - should only be the case, for the first time reaching the outer core
        for cell in self._cellList:
            #start where the last one was found, makes the next search more easier
            isFound,intPoint,self._currentCell = cell.trilinear(x)
            if(isFound): 
       #         print(cell._ID,isFound,"Found outside neighbour domain")                    
                return intPoint
        print("Point not found in dataset")
        
        ## read in data file
    def loadFile(self,path):
        """load a VTK File... Parser is not yet fully complete. Scalar Values in DataSet is not yet implemented"""
        with open(path,'r') as file:
            columnCount = 0
            ver_counter = 0
            cell_counter = 0
            isVertex = False
            isCell = False
            isValue1= False
            dim=0
            value_counter =0
            for line in file:
                columnCount +=1
                #if(columnCount % 100000)==0:print(columnCount)
                ### split by whitespaces
                entries = line.split()
                if not entries: 
                    isVertex = False
                    isCell = False
                    isValue1= False
                    continue
                elif entries[0] == "POINTS":
                    isVertex = True
                    isCell = False
                    isValue1= False                    
                    self._numVert=int(entries[1])
                    print('read ',self._numVert,' Verticies')
                elif entries[0] == "CELLS":
                    isVertex = False
                    isCell = True
                    isValue1= False
                    self._numCells=int(entries[1])
                    print('read ',self._numCells ,' Cells')
                elif entries[0] == "CELL_TYPES":
                    print("Skip Cell Types")
                    isVertex = False
                    isCell = False
                    isValue1= False
                    continue
                elif columnCount == 5929872:
                    print("Test.... end of file reached")
                elif entries[0] == "POINT_DATA":
                    print("Point Data")                    
                    isVertex = False
                    isCell = False
                    isValue1= False
                    continue
                elif entries[0] == "VECTORS":
                    ##skip the line
                    if entries[1] == "current_density":
                        print("stop file reading")
                        break
                    print("Reading Values")
                    isVertex = False
                    isCell = False
                    isValue1= True
                    dim = 3
                    value_counter =0
                    self._valueNames.append(entries[1])
                    
                elif isVertex and entries[0]!="POINTS":
                    ver_counter+=1
                    #print(entries[0])
                    pos = Point3D(float(entries[0]),float(entries[1]),float(entries[2]))
                    ## changed later
                    mag = Point3D(0.0,0.0,0.0)
                    current =  Vertex(ver_counter,pos,mag)
                    self._vertexList.append(current)
                    if ((ver_counter*100.0/self._numVert) % 10)==0:
                        print('Vertex loading reached: ' + str(ver_counter*100.0/self._numVert) + '%')
               #     print('read'+str(ver_counter)+'Vertex')
                elif isCell and entries[0]!="CELLS":                 
                    cell_counter+=1
                    currentCell = Cell(cell_counter,self._vertexList[int(entries[1])-1],self._vertexList[int(entries[2])-1],self._vertexList[int(entries[3])-1],self._vertexList[int(entries[4])-1],self._vertexList[int(entries[5])-1],self._vertexList[int(entries[6])-1],self._vertexList[int(entries[7])-1],self._vertexList[int(entries[8])-1])
                    self._cellList.append(currentCell)
                    ##add link from vertex to cell
                    for ki in range(0,8,1):
                        self._vertexList[int(entries[ki])-1]._partOfCell.append(cell_counter)
                        #print("lenght of part of Cell List for ID: ",int(entries[ki]), len(self._vertexList[int(entries[ki])-1]._partOfCell))
                    self._currentCell = currentCell
                    if ((cell_counter*100.0/self._numCells) % 10)==0:
                        print('Cell loading reached: ' + str(cell_counter*100.0/self._numCells) + '%')
               
                elif isValue1 and entries[0]!="VECTORS":
                  # print('read Values')
                   value_counter+=1
                   if(dim ==3):
                       self._vertexList[value_counter-1]._mag = Point3D(float(entries[0]),float(entries[1]),float(entries[2]))
                   if(dim==1):
                       ##add scalar if needed
                       skippy =2
                   if ((value_counter*100.0/self._numVert) % 10)==0:
                        print('Values loading reached: ' + str(value_counter*100.0/self._numVert) + '%')
     #   if(file.closed):self.computeCellTopology()

    def computeCellTopology(self):
        """ Computes a basic cell topology, based on cells which are sharing a vertex
            Takes very very long and consumes a lot of memory"""
        cellcount=1
        print("compute Cell topology .... ")
        for cell in self._cellList:
            #contains all neighbours, point, line and face neighbours
            neighbourList=[]
            for ver in cell._verts:
                neighbourList.extend(ver._partOfCell)
                #free memory and delete partofCell list
                #del(ver._partofCell[:])
            ##remove duplicates with converting list to a set and back
            #print(len(neighbourList))
            cell._neighbours=list(set(neighbourList))
            cellcount+=1
            #if(cellcount % 100)==0:print(cellcount)
            if((cellcount*100.0/self._numCells)%10) ==0:
                print('Cell topology computition reached: ' + str(cellcount*100.0/self._numCells) + '%')
                print("NeighbourList for Cell :", cell._ID,len(cell._neighbours))
        for cell in self._cellList:
            cell.computeFaceNeighbours(self)
            if len(set(cell._verts))<8:
                print(set(cell._verts), cell._ID)
        #free memory and delete partofCell list
        #li_length=[]
        #for vertex in self._vertexList:
         #   li_length.append(len(vertex._partOfCell))    
          #  del(vertex._partOfCell[:])
        #print(max(li_length))

    


class AvsUcdAscii:
#---class Attributes should only be readed
    _firstVertex = None
    _lastVertex = None
    _currentCell = None

    _lastCell = None
    _numVert =0
    _numCells = 0
    _dim = 0
    _vertexList =[]
    _cellList =[]
#---class Methods
    # Constructor method
    def __init__(self):
          return

            
    ## returns the interpolated value at position x
    def getValue(self,x,dt):
        print("Search for Cell with point: ", toSpherical(x))
        isFound,intPoint,self._currentCell = self._currentCell.trilinear(x)
        if(isFound): 
           # print(self._currentCell._ID,isFound,"Found in Current Cell")                    
            return intPoint        
       # print("Searching in neighbours")
        for neighbourID in self._currentCell._neighbours:
            neighbour = self._cellList[neighbourID-1]
        #for i in range(0,len(self._currentCell._neighbours)):
         #   neighbourID = self._currentCell._neighbours[i]
          #  neighbour = self._cellList[neighbourID-1]
            isFound,intPoint,self._currentCell = neighbour.trilinear(x)
            if(isFound): 
       #         print(neighbour._ID,isFound,"Found in neighbour List")                    
                return intPoint
      #  print("Searching in whole data set")
        #if not found in neighbouring cell - should only be the case, for the first time reaching the outer core
        for cell in self._cellList:
            #start where the last one was found, makes the next search more easier
            isFound,intPoint,self._currentCell = cell.trilinear(x)
            if(isFound): 
       #         print(cell._ID,isFound,"Found outside neighbour domain")                    
                return intPoint
        print("Point not found in dataset")
    
        

    ## read in data file
    def loadFile(self,path):
        with open(path,'r') as file:
            firstline = file.readline()
            firstentries = firstline.split()
            self._numVert = int(firstentries[0])
            self._numCells = int(firstentries[1])
            self._dim = int(firstentries[2])
            print('#Vertices:'+ str(self._numVert))
            print('#Cells:'+ str(self._numCells))

            current = Vertex(0,Point3D(0,0,0),Point3D(0,0,0))
            currentCell = Cell(0,current,current,current,current,current,current,current,current)
            ver_counter = 1
            cell_counter = 1
            value_counter = 1
            i = 1
            for line in file:
                ### split by whitespaces
                entries = line.split()
                if ver_counter is 1:
                    print('read first Vertex')
                    ID = int(entries[0])
                    pos = Point3D(float(entries[1]),float(entries[2]),float(entries[3]))
                    ## changed later
                    mag = Point3D(0.0,0.0,0.0)
                    current =  Vertex(ID,pos,mag)
                    self._firstVertex = current
                    self._vertexList.append(current)
                    ver_counter+=1
                elif ver_counter <= self._numVert:
                    if (ver_counter*100/self._numVert)>=(10*i):
                        print('Vertex loading reached: ' + str(10*i) + '%')
                        i+=1
               #     print('read'+str(ver_counter)+'Vertex')
                    ID = int(entries[0])
                    pos = Point3D(float(entries[1]),float(entries[2]),float(entries[3]))
                    ## changed later
                    mag = Point3D(0.0,0.0,0.0)
                    #print('before:'+str(current._ID))
                    self._vertexList.append(Vertex(ID,pos,mag))
                    #print('after: ' +str(current._ID))
                    ver_counter+=1
                elif cell_counter is 1:
                    i=1
                    print('read first Cell')
                    ID = int(entries[0])
                    currentCell = Cell(ID,self._vertexList[int(entries[3])-1],self._vertexList[int(entries[4])-1],self._vertexList[int(entries[5])-1],self._vertexList[int(entries[6])-1],self._vertexList[int(entries[7])-1],self._vertexList[int(entries[8])-1],self._vertexList[int(entries[9])-1],self._vertexList[int(entries[10])-1])
                    self._cellList.append(currentCell)
                    ##add link from vertex to cell
                    for ki in range(3,11,1):
                        self._vertexList[int(entries[ki])-1]._partOfCell.append(ID)
                        #print("lenght of part of Cell List for ID: ",int(entries[ki]), len(self._vertexList[int(entries[ki])-1]._partOfCell))
                    self._currentCell = currentCell

                    cell_counter +=1
                elif cell_counter <= self._numCells:
                    if (cell_counter*100/self._numCells)>=(10*i):
                        print('Cell loading reached: ' + str(10*i) + '%')
                        i+=1
                    ID = int(entries[0])
                    self._cellList.append(Cell(ID,self._vertexList[int(entries[3])-1],self._vertexList[int(entries[4])-1],self._vertexList[int(entries[5])-1],self._vertexList[int(entries[6])-1],self._vertexList[int(entries[7])-1],self._vertexList[int(entries[8])-1],self._vertexList[int(entries[9])-1],self._vertexList[int(entries[10])-1]))
                    ##add link from vertex to cell
                    for ki in range(3,11,1):
                        self._vertexList[int(entries[ki])-1]._partOfCell.append(ID)
                        ## Kill duplicates
                        tempList=list(set(self._vertexList[int(entries[ki])-1]._partOfCell))
                        self._vertexList[int(entries[ki])-1]._partOfCell=tempList
                        #print("lenght of part of Cell List for ID: ",int(entries[ki]), len(tempList))

                    cell_counter +=1
                elif value_counter is 1:
                    ## skip the line
                    value_counter +=1
                    i=1
                elif value_counter is 2:
                    ##skip the line
                    value_counter +=1
                    print('read first value')
                elif value_counter <=(self._numVert+2):
                    if (value_counter*100/self._numVert)>=(10*i):
                        print('Values loading reached: ' + str(10*i) + '%')
                        i+=1
                    v = self._vertexList[(int(entries[0]))-1]
                    v._mag = Point3D(float(entries[1]),float(entries[2]),float(entries[3]))
                    value_counter +=1
            self.computeCellTopology()

    def computeCellTopology(self):
        cellcount=1
        i=1
        for cell in self._cellList:
            if(cellcount*100/self._numCells)>=(10*i):
                print('Cell topology computition reached: ' + str(10*i) + '%')
                i+=1
            #contains all neighbours, point, line and face neighbours
            neighbourList=[]
            for ver in cell._verts:
                neighbourList.extend(ver._partOfCell)
                #free memory and delete partofCell list
                #del(ver._partofCell[:])
            ##remove duplicates with converting list to a set and back
            cell._neighbours=list(set(neighbourList))
            cellcount+=1
           # print("NeighbourList for Cell :", cell._ID,len(cell._neighbours))
        #free memory and delete partofCell list
        li_length=[]
        for vertex in self._vertexList:
            li_length.append(len(vertex._partOfCell))    
            del(vertex._partOfCell[:])
        print(max(li_length))
    ## return Vertex from his ID
    def getVertex(self,ID):
        vertex = self._firstVertex
        while vertex._ID != int(ID):
            vertex = vertex._nextVertex
        return vertex


