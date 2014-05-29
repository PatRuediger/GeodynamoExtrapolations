#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Mon May 12 10:46:55 2014

@author: Patrick

Data structure for Visualization
Using hexahedron Grid cells with 8 vertices
"""
import math


def mainTest():
    print('Main startet')
    DS = AvsUcdAscii()
    DS.loadFile('C:\out.1550.inp')

def testInterpolation():
    v0 = Vertex(0,Point3D(0,0,0),Point3D(0,0,0))
    v1 = Vertex(1,Point3D(0,1,0),Point3D(0,1,0))
    v2 = Vertex(2,Point3D(1,1,0),Point3D(1,1,0))
    v3 = Vertex(3,Point3D(1,0,0),Point3D(1,0,0))
    v4 = Vertex(4,Point3D(0,0,1),Point3D(0,0,1))
    v5 = Vertex(5,Point3D(0,1,1),Point3D(0,0,1))
    v6 = Vertex(6,Point3D(1,1,1),Point3D(1,1,1))
    v7 = Vertex(7,Point3D(1,1,0),Point3D(1,1,0))

    cell = Cell(0,v0,v1,v2,v3,v4,v5,v6,v7)
    p= [0.5,0.5,0.5]
    output = cell.isInside(p)
    print(output)
    temp = cell.trilinear(p)
    output = [temp._x,temp._y, temp._z]
    print(output)

"""Dot/Scalar Product for 3D Vectors"""
def dot(a,b):
        a0 = Point3D(0,0,0)
        b0 = Point3D(0,0,0)
        a0 = a
        b0 = b
        return a0._x*b0._x + a0._y*b0._y + a0._z*b0._z


class Point3D:
#    _x = 0.0
#    _y = 0.0
#    _z = 0.0

    def _length(self):
        return math.sqrt(self._x**2.0+self._y**2.0+self._z**2.0)

    def __init__(self,x,y,z):
        self._x = x
        self._y = y
        self._z = z

    def mult(self,scalar):
        return Point3D(self._x * scalar,self._y * scalar,self._z * scalar)

    def add(self,vec):
        return Point3D(self._x + vec._x, self._y + vec._y, self._z + vec._z)

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

#---class Methods
    # Constructor method
    def __init__(self,ID,pos,mag):
        self._ID = ID
        self._pos = pos
        self._mag = mag

    def addVertex(self,ID,pos,mag):
        newVertex = Vertex(ID,pos,mag)
        self._nextVertex = newVertex
        return self._nextVertex



#---- Cell is ment to be an hexahedron
class Cell:
#---class Attributes should only be readed
    _ID = 0
    # should be 8 everytime
    # ordering: c000,c001,c101,c100,c010,c011,c111,c110
    _verts =[]
    #Pointer to next Cell
    _nextCell = None

#---class Methods
    # Constructor method
    def __init__(self,ID,v0,v1,v2,v3,v4,v5,v6,v7):
        self._verts = [v0,v1,v2,v3,v4,v5,v6,v7]
        self._ID = ID

    def addCell(self,ID,v0,v1,v2,v3,v4,v5,v6,v7):
        newCell = Cell(ID,v0,v1,v2,v3,v4,v5,v6,v7)
        self._nextCell = newCell
        return self._nextCell


    ## @input is an array  with 3 entries
    def isInside(self,p):
        bounds = self.boundaries()
        if bounds[0]<=p[0] and bounds[1]<=p[1] and bounds[2]<=p[2] and bounds[3]>=p[0] and bounds[4]>=p[1] and bounds[5]>=p[2]:
            return True
        else:
            return False


    ## @input is an array  with 3 entries
    ## @output is a Point3D
    def trilinear(self,x):
        if(self.isInside(x)):
            xd = (x[0]-self._verts[0]._pos._x)/(self._verts[6]._pos._x-self._verts[0]._pos._x)
            yd = (x[1]-self._verts[0]._pos._y)/(self._verts[6]._pos._y-self._verts[0]._pos._y)
            zd = (x[2]-self._verts[0]._pos._z)/(self._verts[6]._pos._z-self._verts[0]._pos._z)

            c00 =  Point3D(0,0,0)
            c10 =  Point3D(0,0,0)
            c01 =  Point3D(0,0,0)
            c11 =  Point3D(0,0,0)
            c0 =  Point3D(0,0,0)
            c1 =  Point3D(0,0,0)
            c =  Point3D(0,0,0)

            c00 = self._verts[0]._mag.mult(1.0-xd)
            c00.add(self._verts[3]._mag.mult(xd))
            c10 = self._verts[4]._mag.mult(1.0-xd)
            c10.add(self._verts[7]._mag.mult(xd))
            c01 = self._verts[1]._mag.mult(1.0-xd)
            c01.add(self._verts[2]._mag.mult(xd))
            c11 = self._verts[5]._mag.mult(1.0-xd)
            c11.add(self._verts[6]._mag.mult(xd))

            c0= c00.mult(1.0-yd)
            c0.add(c10.mult(yd))
            c1=c01.mult(1.0-yd)
            c1.add(c11.mult(yd))

            c = c0.mult(1.0-zd)
            c.add(c1.mult(zd))
            return c
        elif(self._nextCell is not None):
            self._nextCell.trilinear(x)
        else:
            return None

    def boundaries(self):
        min_x = 10000000000000000.0
        min_y = 10000000000000000.0
        min_z = 10000000000000000.0
        max_x = 0.0
        max_y = 0.0
        max_z = 0.0
        for v in self._verts:
            if min_x >= v._pos._x:
                min_x = v._pos._x
            if min_y >= v._pos._y:
                min_y = v._pos._y
            if min_z >= v._pos._z:
                min_z = v._pos._z
            if max_x <= v._pos._x:
                max_x = v._pos._x
            if max_y <= v._pos._y:
                max_y = v._pos._y
            if max_z <= v._pos._z:
                max_z = v._pos._z
        bounds = [min_x,min_y,min_z,max_x,max_y,max_z]
        #print(bounds)
        return bounds


class AvsUcdAscii:
#---class Attributes should only be readed
    _firstVertex = None
    _lastVertex = None
    _firstCell = None
    _lastCell = None
    _numVert =0
    _numCells = 0
    _dim = 0
    _vertexList =[]
    _cellList =[]
#---class Methods
    # Constructor method
#    def __init__(self):
#          return
    ## returns the interpolated value at position x
    def getValue(self,x):
        if self._firstCell is not None:
            return self._firstCell.trilinear(x)

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
                    nextV = current.addVertex(ID,pos,mag)
                    self._vertexList.append(nextV)
                    current = nextV
                    #print('after: ' +str(current._ID))
                    ver_counter+=1
                elif cell_counter is 1:
                    i=1
                    print('read first Cell')
                    ID = int(entries[0])
                  #  currentCell = Cell(ID,self.getVertex(entries[3]),self.getVertex(entries[4]),self.getVertex(entries[5]),self.getVertex(entries[6]),self.getVertex(entries[7]),self.getVertex(entries[8]),self.getVertex(entries[9]),self.getVertex(entries[10]))
                    currentCell = Cell(ID,self._vertexList[int(entries[3])-1],self._vertexList[int(entries[4])-1],self._vertexList[int(entries[5])-1],self._vertexList[int(entries[6])-1],self._vertexList[int(entries[7])-1],self._vertexList[int(entries[8])-1],self._vertexList[int(entries[9])-1],self._vertexList[int(entries[10])-1])
                    self._firstCell = currentCell
                    cell_counter +=1
                elif cell_counter <= self._numCells:
                    if (cell_counter*100/self._numCells)>=(10*i):
                        print('Cell loading reached: ' + str(10*i) + '%')
                        i+=1
                    ID = int(entries[0])
                    #nextC= currentCell.addCell(ID,self.getVertex(entries[3]),self.getVertex(entries[4]),self.getVertex(entries[5]),self.getVertex(entries[6]),self.getVertex(entries[7]),self.getVertex(entries[8]),self.getVertex(entries[9]),self.getVertex(entries[10]))
                    nextC= currentCell.addCell(ID,self._vertexList[int(entries[3])-1],self._vertexList[int(entries[4])-1],self._vertexList[int(entries[5])-1],self._vertexList[int(entries[6])-1],self._vertexList[int(entries[7])-1],self._vertexList[int(entries[8])-1],self._vertexList[int(entries[9])-1],self._vertexList[int(entries[10])-1])
                    currentCell = nextC
                    cell_counter +=1
                elif value_counter is 1:
                    ## skip the line
                    value_counter +=1
                    i=1
                elif value_counter is 2:
                    ##skip the line
                    value_counter +=1
                    print('read first value')
                elif value_counter <=(self._numCells+2):
                    if (value_counter*100/self._numVert)>=(10*i):
                        print('Values loading reached: ' + str(10*i) + '%')
                        i+=1
                    v = self._vertexList[(int(entries[0]))-1]
                    v._mag = Point3D(float(entries[1]),float(entries[2]),float(entries[3]))
                    value_counter +=1

    ## return Vertex from his ID
    def getVertex(self,ID):
        vertex = self._firstVertex
        while vertex._ID != int(ID):
            vertex = vertex._nextVertex
        return vertex


