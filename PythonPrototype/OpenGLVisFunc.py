
from datastructure import Point3D
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys
import copy
from math import cos, sin

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

g_streamLinePos = []
g_streamLineValue = []

g_vecFieldPos = []
g_vecFieldValue = []

g_streamLineList = []


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

def arrowGlyph(pos,value,scale):
    #Normalize the vector
    value_n = Point3D(value._x/value._length(),value._y/value._length(),value._z/value._length())
    glPointSize( 5.0 );
    glBegin(GL_POINTS);
    glVertex3f(pos._x,pos._y,pos._z);
    glEnd();
    
    glLineWidth(1.0);
    glBegin(GL_LINES);
    glVertex3f(pos._x,pos._y,pos._z);
    glVertex3f(pos._x+value_n._x*scale,pos._y+value_n._y*scale,pos._z+value_n._z*scale);
    glEnd();
    return

def streamLine(pos,value):
    #// Pos,value are Point3D
    glLineWidth(2.0);
    for i in range(len(pos)-1):
        glBegin( GL_LINES );
        glVertex3f (pos[i]._x,pos[i]._y,pos[i]._z);
        glVertex3f (pos[i+1]._x,pos[i+1]._y,pos[i+1]._z);
        glEnd();
    return

def drawVectorfield(xf,vf,scale=5):
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
    glLineWidth(1.0)
    global g_streamLineList
    for sl in g_streamLineList:
        for i in range(len(sl[0])-1):
            glBegin( GL_LINES );
            glVertex3f (sl[0][i]._x,sl[0][i]._y,sl[0][i]._z);
            glVertex3f (sl[0][i+1]._x,sl[0][i+1]._y,sl[0][i+1]._z);
            glEnd();
    return
            
def setVectorfield(xf,vf):
    global g_vecFieldPos
    g_vecFieldPos = xf
    global g_vecFieldValue
    g_vecFieldValue = vf
    return

def Draw ():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				# // Clear Screen And Depth Buffer
    glLoadIdentity();												# // Reset The Current Modelview Matrix
    glTranslatef(0.0,0.0,-15.0);									
    glPushMatrix();													# // NEW: Prepare Dynamic Transform
    glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform
    glColor3f(0.8,0.75,1.0);
    drawVectorfield(g_vecFieldPos,g_vecFieldValue,0.5)
    drawStreamLines()
    #streamLine(g_streamLinePos,g_streamLineValue)
    glPopMatrix();													# // NEW: Unapply Dynamic Transform
    
    glLoadIdentity();												# // Reset The Current Modelview Matrix
    glTranslatef(0.0,0.0,-15.0);										
    glPushMatrix();													# // NEW: Prepare Dynamic Transform
    glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform
    
    #Define Sphere here!!
    glColor3f(1.0,0.75,0.75);
    glPolygonMode(GL_FRONT_AND_BACK,GL_LINE)
    glLineWidth(1.0)
    gluSphere(g_quadratic,2.92,20,20);
    
    #Define Sphere here!!
    glColor3f(1.0,0.1,0.1);
    glPolygonMode(GL_FRONT_AND_BACK,GL_FILL)
    gluSphere(g_quadratic,0.53,20,20);
    
    glPopMatrix();													# // NEW: Unapply Dynamic Transform
    glFlush ();														# // Flush The GL Rendering Pipeline
    glutSwapBuffers()
    return

