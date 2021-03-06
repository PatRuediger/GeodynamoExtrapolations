# /***********************************************
# *                                              *
# *    Jeff Molofee's Revised OpenGL Basecode    *
# *  Huge Thanks To Maxwell Sayles & Peter Puck  *
# *            Adapted to PyOpenGL GLUT          *
# *                                              *
# ***********************************************/


from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GLUT.freeglut import*
import sys
import math

from OpenGLVisFunc import *		# Draw (), Initialize () and all the real OpenGL work.
from ArcBall import *		# // *NEW* ArcBall header


# *********************** Globals *********************** 
# Python 2.2 defines these directly
try:
	True
except NameError:
	True = 1==1
	False = 1==0


# Some api in the chain is translating the keystrokes to this octal string
# so instead of saying: ESCAPE = 27, we use the following.
ESCAPE = '\033'

# Number of the glut window.
window = 0
g_zoomfactor = 45.0
g_Width = 1024
g_Height = 768





# A general OpenGL initialization function.  Sets all of the initial parameters. 
def InitGL(Width, Height):				# We call this right after our OpenGL window is created.

	glShadeModel(GL_SMOOTH)				# Enables Smooth Color Shading
	glClearColor(0.0, 0.0, 0.0, 0.5)	# This Will Clear The Background Color To Black
	glClearDepth(1.0)					# Enables Clearing Of The Depth Buffer
	glEnable(GL_DEPTH_TEST)				# Enables Depth Testing
	glDepthFunc(GL_LEQUAL)				# The Type Of Depth Test To Do
	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST) # Really Nice Perspective Calculations



	return True									# // Initialization Went OK


# Reshape The Window When It's Moved Or Resized
def ReSizeGLScene(Width, Height):
	if Height == 0:						# Prevent A Divide By Zero If The Window Is Too Small 
		Height = 1

	glViewport(0, 0, Width, Height)		# Reset The Current Viewport And Perspective Transformation
	glMatrixMode(GL_PROJECTION)			# // Select The Projection Matrix
	glLoadIdentity()					# // Reset The Projection Matrix
	# // field of view, aspect ratio, near and far
	# This will squash and stretch our objects as the window is resized.
	# Note that the near clip plane is 1 (hither) and the far plane is 1000 (yon)
    	global g_Height
   	global g_Width
   	global g_zoomfactor
 	g_Width = Width
   	g_Height = Height
	gluPerspective(g_zoomfactor, float(g_Width)/float(g_Height), 1, 100.0)
         
	glMatrixMode (GL_MODELVIEW);		# // Select The Modelview Matrix
	glLoadIdentity ();					# // Reset The Modelview Matrix
	g_ArcBall.setBounds (Width, Height)	# //*NEW* Update mouse bounds for arcbal
	return

def mouseWheel(button,  direction ,  x,  y):
    global g_zoomfactor
    if (direction > 0): #zoom in
        g_zoomfactor= g_zoomfactor*0.8
    elif (direction < 0):   #zoom out
        g_zoomfactor=g_zoomfactor*1.2
    ReSizeGLScene(g_Width,g_Height)
    return

# The function called whenever a key is pressed. Note the use of Python tuples to pass in: (key, x, y)  
def keyPressed(*args):
	global g_quadratic
	# If escape is pressed, kill everything.
	key = args [0]
	if key == ESCAPE:
		gluDeleteQuadric (g_quadratic)
		sys.exit ()



def main():
	global window
	# pass arguments to init
	glutInit(sys.argv)
      
	# Select type of Display mode:   
	#  Double buffer 
	#  RGBA color
	# Alpha components supported 
	# Depth buffer
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH)
	
	# get a 640 x 480 window 
	glutInitWindowSize(1024, 768)
	
	# the window starts at the upper left corner of the screen 
	glutInitWindowPosition(0, 0)
	
	# Okay, like the C version we retain the window id to use when closing, but for those of you new
	# to Python, remember this assignment would make the variable local and not global
	# if it weren't for the global declaration at the start of main.
	window = glutCreateWindow("Rotate sphere")

   	# Register the drawing function with glut, BUT in Python land, at least using PyOpenGL, we need to
	# set the function pointer and invoke a function to actually register the callback, otherwise it
	# would be very much like the C version of the code.	
	glutDisplayFunc(Draw)
	
	# Uncomment this line to get full screen.
	#glutFullScreen()

	# When we are doing nothing, redraw the scene.
	glutIdleFunc(Draw)
	
	# Register the function called when our window is resized.
	glutReshapeFunc(ReSizeGLScene)
	
	# Register the function called when the keyboard is pressed.  
	glutKeyboardFunc(keyPressed)


	# GLUT When mouse buttons are clicked in window
	glutMouseFunc (Upon_Click)
      # GLUT When wheel used
	glutMouseWheelFunc (mouseWheel)

	# GLUT When the mouse mvoes
	glutMotionFunc (Upon_Drag)
 
	# We've told Glut the type of window we want, and we've told glut about
	# various functions that we want invoked (idle, resizing, keyboard events).
	# Glut has done the hard work of building up thw windows DC context and 
	# tying in a rendering context, so we are ready to start making immediate mode
	# GL calls.
	# Call to perform inital GL setup (the clear colors, enabling modes
	Initialize (1024, 768)

	# Start Event Processing Engine	
	glutMainLoop()

# Print message to console, and kick off the main to get it rolling.
if __name__ == "__main__":
	print "Hit ESC key to quit."
	main()

