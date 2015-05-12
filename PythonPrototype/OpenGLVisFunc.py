
from datastructure import Point3D
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys
import copy
from math import cos, sin
from sklearn.externals.joblib import Parallel, delayed
import multiprocessing
import buffers
import numpy

from ArcBall import * 				# ArcBallT 

PI2 = 2.0*3.1415926535			# 2 * PI (not squared!) 		// PI Squared

# *********************** Globals *********************** 
# Python 2.2 defines these directly
try:
	True
except NameError:
	True = 1==1
	False = 1==0

g_Transform = Matrix4fT ()
g_LastRot = Matrix3fT ()
g_ThisRot = Matrix3fT ()

g_ArcBall = ArcBallT (1024,768)
g_isDragging = False
g_quadratic = None

lmax = []

g_cm_rb = (0.0,0.0)
g_streamLinePos = []
g_streamLineValue = []

g_vecFieldPos = []
g_vecFieldValue = []
g_vecFieldPos2 = []
g_vecFieldValue2 = []

g_streamLineList = []
g_cellList = []
g_VertexBuffer = None
g_CellBuffer = None


# A general OpenGL initialization function.  Sets all of the initial parameters. 
def Initialize (Width, Height):				# We call this right after our OpenGL window is created.
	global g_quadratic

	glClearColor(0.0, 0.0, 0.0, 1.0)					# This Will Clear The Background Color To Black
	glClearDepth(1.0)									# Enables Clearing Of The Depth Buffer
	glDepthFunc(GL_LEQUAL)								# The Type Of Depth Test To Do
	glEnable(GL_DEPTH_TEST)								# Enables Depth Testing
	glShadeModel (GL_FLAT);								# Select Flat Shading (Nice Definition Of Objects)
	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST) 	# Really Nice Perspective Calculations

	g_quadratic = gluNewQuadric();
	gluQuadricNormals(g_quadratic, GLU_SMOOTH);
	gluQuadricDrawStyle(g_quadratic, GLU_FILL); 


	glEnable (GL_LIGHT0)
	glEnable (GL_LIGHTING)

	glEnable (GL_COLOR_MATERIAL)

	return True

def initVertexBuffer(vertexList):
    global g_VertexBuffer
    vertexListTripples = []
    ## List Index is the same as in _vertexList
    for elem in vertexList:
        vertexListTripples.append([elem._pos._x,elem._pos._y,elem._pos._z])
    numpy_verts = numpy.array(vertexListTripples,dtype=numpy.float32)
    g_VertexBuffer = vertexBuffer(numpy_verts,GL_STATIC_DRAW)
      
def Upon_Drag (cursor_x, cursor_y):
	""" Mouse cursor is moving
		Glut calls this function (when mouse button is down)
		and pases the mouse cursor postion in window coords as the mouse moves.
	"""
	global g_isDragging, g_LastRot, g_Transform, g_ThisRot

	if (g_isDragging):
		mouse_pt = Point2fT (cursor_x, cursor_y)
		ThisQuat = g_ArcBall.drag (mouse_pt)						# // Update End Vector And Get Rotation As Quaternion
		g_ThisRot = Matrix3fSetRotationFromQuat4f (ThisQuat)		# // Convert Quaternion Into Matrix3fT
		# Use correct Linear Algebra matrix multiplication C = A * B
		g_ThisRot = Matrix3fMulMatrix3f (g_LastRot, g_ThisRot)		# // Accumulate Last Rotation Into This One
		g_Transform = Matrix4fSetRotationFromMatrix3f (g_Transform, g_ThisRot)	# // Set Our Final Transform's Rotation From This One
	return

def Upon_Click (button, button_state, cursor_x, cursor_y):
	""" Mouse button clicked.
		Glut calls this function when a mouse button is
		clicked or released.
	"""
	global g_isDragging, g_LastRot, g_Transform, g_ThisRot

	g_isDragging = False
	if (button == GLUT_RIGHT_BUTTON and button_state == GLUT_UP):
		# Right button click
		g_LastRot = Matrix3fSetIdentity ();							# // Reset Rotation
		g_ThisRot = Matrix3fSetIdentity ();							# // Reset Rotation
		g_Transform = Matrix4fSetRotationFromMatrix3f (g_Transform, g_ThisRot);	# // Reset Rotation
	elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_UP):
		# Left button released
		g_LastRot = copy.copy (g_ThisRot);							# // Set Last Static Rotation To Last Dynamic One
	elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_DOWN):
		# Left button clicked down
		g_LastRot = copy.copy (g_ThisRot);							# // Set Last Static Rotation To Last Dynamic One
		g_isDragging = True											# // Prepare For Dragging
		mouse_pt = Point2fT (cursor_x, cursor_y)
		g_ArcBall.click (mouse_pt);								# // Update Start Vector And Prepare For Dragging

	return



def Torus(MinorRadius, MajorRadius):		
    # // Draw A Torus With Normals
    glBegin( GL_TRIANGLE_STRIP );									# // Start A Triangle Strip
    for i in xrange (20): 											# // Stacks
        for j in xrange (-1, 20): 										# // Slices
			# NOTE, python's definition of modulus for negative numbers returns
			# results different than C's
			#       (a / d)*d  +  a % d = a
            if (j < 0):
                wrapFrac = (-j%20)/20.0
                wrapFrac *= -1.0
            else:
                wrapFrac = (j%20)/20.0;
            phi = PI2*wrapFrac;
            sinphi = sin(phi);
            cosphi = cos(phi);
            r = MajorRadius + MinorRadius*cosphi;
            glNormal3f (sin(PI2*(i%20+wrapFrac)/20.0)*cosphi, sinphi, cos(PI2*(i%20+wrapFrac)/20.0)*cosphi);
            glVertex3f (sin(PI2*(i%20+wrapFrac)/20.0)*r, MinorRadius*sinphi, cos(PI2*(i%20+wrapFrac)/20.0)*r);
            
            glNormal3f (sin(PI2*(i+1%20+wrapFrac)/20.0)*cosphi, sinphi, cos(PI2*(i+1%20+wrapFrac)/20.0)*cosphi);
            glVertex3f (sin(PI2*(i+1%20+wrapFrac)/20.0)*r, MinorRadius*sinphi, cos(PI2*(i+1%20+wrapFrac)/20.0)*r);
    glEnd();														# // Done Torus
    return
 
def Sphere(center,radius):
    return

def built_cm_rainbow(data):
    le =[]
    global g_cm_rb
    verts = data._vertexList
    for v in verts:
        le.append(v._mag._length())
    max_v = max(le)
    min_v = min(le)
    g_cm_rb=(max_v,min_v)
    return

def cm_rainbow(x):
    i = x/(g_cm_rb[0]-g_cm_rb[1])
    r = sin(0.3*i + 0)*0.49 + 0.51
    g = sin(0.3*i + PI2/3.0)*0.49 + 0.51
    b = sin(0.3*i + 2*PI2/3.0)*0.49 + 0.51
    #print("Value: ",x , "Color: ", r,g,b)
    return r,g,b
    
    
def arrowGlyph(pos,value,scale):
    #Normalize the vector
    maxl=max(lmax)
    if value == None:
        return
    value_n = Point3D(value._x/value._length(),value._y/value._length(),value._z/value._length())
    glPointSize( 5.0 );
    glBegin(GL_POINTS);
    glVertex3f(pos._x,pos._y,pos._z);
    glEnd();

    vlen = 0.1 +( scale * value._length() / maxl)   
   # vlen =scale    
    
    glLineWidth(1.0);
    glBegin(GL_LINES);
    glVertex3f(pos._x,pos._y,pos._z);
    glVertex3f(pos._x+value_n._x*vlen,pos._y+value_n._y*vlen,pos._z+value_n._z*vlen);
    glEnd();
    return

def streamLine(pos,value):
    #// Pos,value are Point3D
    glLineWidth(0.5);
    for i in range(len(pos)-1):
        glBegin( GL_LINES );
        glVertex3f (pos[i]._x,pos[i]._y,pos[i]._z);
        glVertex3f (pos[i+1]._x,pos[i+1]._y,pos[i+1]._z);
        glEnd();
    return

def drawVectorfield(xf,vf,scale):
    """input are two arrays of Point3D and a scaling factor"""
    for i in range(len(xf)-1):
        arrowGlyph(xf[i],vf[i],scale)
    return

def setStreamLine(x,v):
    global g_streamLinePos
    g_streamLinePos = x
    global g_streamLineValue
    g_streamLineValue = v
    return

"""adding a whole streamline as a tupple of positions and values"""
def addStreamLine(xf,vf):
    sl = (xf,vf)
    global g_streamLineList
    g_streamLineList.append(sl)
    return

def drawStreamLines():
    glLineWidth(0.8)
    global g_streamLineList
    for sl in g_streamLineList:
        for i in range(len(sl[0])-1):
           # val = sl[1][i]._length()
           # r,g,b = cm_rainbow(val)
           # glColor4f(r,g,b,1.0);
            glColor4f(1.0,0.1,0.1,1.0)
            glBegin( GL_LINES );
            glVertex3f (sl[0][i]._x,sl[0][i]._y,sl[0][i]._z);
            glVertex3f (sl[0][i+1]._x,sl[0][i+1]._y,sl[0][i+1]._z);
            glEnd();
    return

def setCellList(cellList):
    global g_cellList
    g_cellList = cellList

def drawSingleWire(cell):
    #print("Drawing cell with id: ", cell._ID)
    glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    glColor4f(1.0,1.0,1.0,0.9)
    glBegin(GL_QUADS);
    glVertex3f(cell._verts[0]._pos._x,cell._verts[0]._pos._y,cell._verts[0]._pos._z)
    glVertex3f(cell._verts[1]._pos._x,cell._verts[1]._pos._y,cell._verts[1]._pos._z)
    glVertex3f(cell._verts[2]._pos._x,cell._verts[2]._pos._y,cell._verts[2]._pos._z)
    glVertex3f(cell._verts[3]._pos._x,cell._verts[3]._pos._y,cell._verts[3]._pos._z)
    glEnd();    
    glBegin(GL_QUADS);
    glVertex3f(cell._verts[5]._pos._x,cell._verts[5]._pos._y,cell._verts[5]._pos._z)
    glVertex3f(cell._verts[1]._pos._x,cell._verts[1]._pos._y,cell._verts[1]._pos._z)
    glVertex3f(cell._verts[0]._pos._x,cell._verts[0]._pos._y,cell._verts[0]._pos._z)
    glVertex3f(cell._verts[4]._pos._x,cell._verts[4]._pos._y,cell._verts[4]._pos._z)
    glEnd();
    glBegin(GL_QUADS);
    glVertex3f(cell._verts[6]._pos._x,cell._verts[6]._pos._y,cell._verts[6]._pos._z)
    glVertex3f(cell._verts[7]._pos._x,cell._verts[7]._pos._y,cell._verts[7]._pos._z)
    glVertex3f(cell._verts[4]._pos._x,cell._verts[4]._pos._y,cell._verts[4]._pos._z)
    glVertex3f(cell._verts[5]._pos._x,cell._verts[5]._pos._y,cell._verts[5]._pos._z)
    glEnd();
    glBegin(GL_QUADS);
    glVertex3f(cell._verts[2]._pos._x,cell._verts[2]._pos._y,cell._verts[2]._pos._z)
    glVertex3f(cell._verts[6]._pos._x,cell._verts[6]._pos._y,cell._verts[6]._pos._z)
    glVertex3f(cell._verts[7]._pos._x,cell._verts[7]._pos._y,cell._verts[7]._pos._z)
    glVertex3f(cell._verts[3]._pos._x,cell._verts[3]._pos._y,cell._verts[3]._pos._z)
    glEnd();
    glBegin(GL_QUADS);
    glVertex3f(cell._verts[1]._pos._x,cell._verts[1]._pos._y,cell._verts[1]._pos._z)
    glVertex3f(cell._verts[5]._pos._x,cell._verts[5]._pos._y,cell._verts[5]._pos._z)
    glVertex3f(cell._verts[6]._pos._x,cell._verts[6]._pos._y,cell._verts[6]._pos._z)
    glVertex3f(cell._verts[2]._pos._x,cell._verts[2]._pos._y,cell._verts[2]._pos._z)
    glEnd();
    glBegin(GL_QUADS);
    glVertex3f(cell._verts[0]._pos._x,cell._verts[0]._pos._y,cell._verts[0]._pos._z)
    glVertex3f(cell._verts[3]._pos._x,cell._verts[3]._pos._y,cell._verts[3]._pos._z)
    glVertex3f(cell._verts[7]._pos._x,cell._verts[7]._pos._y,cell._verts[7]._pos._z)
    glVertex3f(cell._verts[4]._pos._x,cell._verts[4]._pos._y,cell._verts[4]._pos._z)
    glEnd();


def drawWireFrame():
    global g_cellList
    glLineWidth(0.1)
   # print("Drawing Wireframe")
    for i in range (len(g_cellList)):
        drawSingleWire(g_cellList[i])
    #num_cores = multiprocessing.cpu_count()
    #Parallel(n_jobs=num_cores)(delayed(drawSingleWire)(cell) for cell in g_cellList)
"""    for cell in g_cellList:
        print("Drawing cell with id: ", cell._ID)
        glColor4f(1.0,0.1,0.1,1.0)
        glBegin(GL_LINES);
        glVertex3f(cell._verts[0]._pos._x,cell._verts[0]._pos._y,cell._verts[0]._pos._z)
        glVertex3f(cell._verts[1]._pos._x,cell._verts[1]._pos._y,cell._verts[1]._pos._z)
        glEnd();    
        glBegin(GL_LINES);
        glVertex3f(cell._verts[0]._pos._x,cell._verts[0]._pos._y,cell._verts[0]._pos._z)
        glVertex3f(cell._verts[3]._pos._x,cell._verts[3]._pos._y,cell._verts[3]._pos._z)
        glEnd();
        glBegin(GL_LINES);
        glVertex3f(cell._verts[0]._pos._x,cell._verts[0]._pos._y,cell._verts[0]._pos._z)
        glVertex3f(cell._verts[4]._pos._x,cell._verts[4]._pos._y,cell._verts[4]._pos._z)
        glEnd();
        glBegin(GL_LINES);
        glVertex3f(cell._verts[1]._pos._x,cell._verts[1]._pos._y,cell._verts[1]._pos._z)
        glVertex3f(cell._verts[5]._pos._x,cell._verts[5]._pos._y,cell._verts[5]._pos._z)
        glEnd();
        glBegin(GL_LINES);
        glVertex3f(cell._verts[1]._pos._x,cell._verts[1]._pos._y,cell._verts[1]._pos._z)
        glVertex3f(cell._verts[2]._pos._x,cell._verts[2]._pos._y,cell._verts[2]._pos._z)
        glEnd();
        glBegin(GL_LINES);
        glVertex3f(cell._verts[3]._pos._x,cell._verts[3]._pos._y,cell._verts[3]._pos._z)
        glVertex3f(cell._verts[2]._pos._x,cell._verts[2]._pos._y,cell._verts[2]._pos._z)
        glEnd();
        glBegin(GL_LINES);
        glVertex3f(cell._verts[3]._pos._x,cell._verts[3]._pos._y,cell._verts[3]._pos._z)
        glVertex3f(cell._verts[7]._pos._x,cell._verts[7]._pos._y,cell._verts[7]._pos._z)
        glEnd();
        glBegin(GL_LINES);
        glVertex3f(cell._verts[2]._pos._x,cell._verts[2]._pos._y,cell._verts[2]._pos._z)
        glVertex3f(cell._verts[6]._pos._x,cell._verts[6]._pos._y,cell._verts[6]._pos._z)
        glEnd();"""
         
def setVectorfield(xf,vf):
    global g_vecFieldPos
    g_vecFieldPos = xf
    global g_vecFieldValue
    g_vecFieldValue = vf
    global lmax
    le =[]
    for v in vf:
        if v != None:
            le.append(v._length())
    lmax.append(max(le))
    return
    
def setVectorfield2(xf,vf):
    global g_vecFieldPos2
    g_vecFieldPos2 = xf
    global g_vecFieldValue2
    g_vecFieldValue2 = vf
    global lmax2
    le =[]
    for v in vf:
        if v != None:
            le.append(v._length())
    lmax.append(max(le))
    return

def Draw ():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				# // Clear Screen And Depth Buffer
    glLoadIdentity();												# // Reset The Current Modelview Matrix
    glTranslatef(0.0,0.0,-15.0);									
    glPushMatrix();													# // NEW: Prepare Dynamic Transform
    glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );
<<<<<<< HEAD
    glColor4f(0.8,0.1,1.0,1.0); #purple
    drawVectorfield(g_vecFieldPos,g_vecFieldValue,0.3)
    glColor4f(0.1,0.8,1.0,1.0); #light blue
    drawVectorfield(g_vecFieldPos2,g_vecFieldValue2,0.3)
    glColor4f(1.0,0.1,0.1,1.0)
   # drawStreamLines()
#    drawWireFrame()
=======
    glColor4f(0.8,0.1,1.0,1.0);
    drawVectorfield(g_vecFieldPos,g_vecFieldValue,0.5)
    glColor4f(0.1,0.8,1.0,1.0);
    drawVectorfield(g_vecFieldPos2,g_vecFieldValue2,0.5)
    #glColor4f(1.0,0.1,0.1,1.0)
    #drawStreamLines()
    #glEnable(GL_CULL_FACE)
    #drawWireFrame()
>>>>>>> origin/test01
    #streamLine(g_streamLinePos,g_streamLineValue)
    glPopMatrix();													# // NEW: Unapply Dynamic Transform
    
    glLoadIdentity();												# // Reset The Current Modelview Matrix
    glTranslatef(0.0,0.0,-15.0);										
    glPushMatrix();													# // NEW: Prepare Dynamic Transform
    glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform
    
    #Define Sphere here!!
    glEnable(GL_LINE_SMOOTH)
    glLineWidth(0.01)
    glColor4f(1.0,0.75,0.75,0.3);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
    glLineWidth(0.01)
    gluSphere(g_quadratic,2.92,20,20);
    
    #Define Sphere here!!
    glEnable(GL_LINE_SMOOTH)
    glLineWidth(0.01)
    glColor4f(1.0,0.5,0.75,0.3);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
    glLineWidth(0.01)
    gluSphere(g_quadratic,1.53,20,20);
    
    #Define Sphere here!!
    glEnable(GL_LINE_SMOOTH)
    glColor4f(1.0,0.1,0.1,0.8);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
    gluSphere(g_quadratic,0.53,20,20);
    
    glPopMatrix();													# // NEW: Unapply Dynamic Transform
    glFlush ();														# // Flush The GL Rendering Pipeline
    glutSwapBuffers()
    return

    