from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys

#  from pyquaternion import Quaternion    ## would be useful for 3D simulation
import numpy as np
import math as math

window = 0  # number of the glut window
theta = 0.0
simTime = 0
dT = 0.01
simRun = True
RAD_TO_DEG = 180.0 / 3.1416
multi = 0
boardY = 0


#####################################################
#### Link class, i.e., for a rigid body
#####################################################

class Link:
    color = [0, 0, 0]  ## draw color
    size = [1, 1, 1]  ## dimensions
    mass = 1.0  ## mass in kg
    izz = 1.0  ## moment of inertia about z-axis
    theta = 0  ## 2D orientation  (will need to change for 3D)
    omega = 0  ## 2D angular velocity
    posn = np.array([0.0, 0.0, 0.0])  ## 3D position (keep z=0 for 2D)
    vel = np.array([0.0, 0.0, 0.0])  ## initial velocity

    def draw(self):  ### steps to draw a link
        glPushMatrix()  ## save copy of coord frame
        glTranslatef(self.posn[0], self.posn[1], self.posn[2])  ## move
        glRotatef(self.theta * RAD_TO_DEG, 0, 0, 1)  ## rotate
        glScale(self.size[0], self.size[1], self.size[2])  ## set size
        glColor3f(self.color[0], self.color[1], self.color[2])  ## set colour
        DrawCube()  ## draw a scaled cube
        glPopMatrix()  ## restore old coord frame

#####################################################
# #### Bar class, i.e., for dynamic illustration
#####################################################
class Bar:
    color = [0, 0, 0]
    size = [1, 1, 1]
    posn = np.array([5.0, 0.0, 0.0])

    def draw(self):
        glPushMatrix()  ## save copy of coord frame
        glTranslatef(self.posn[0], self.posn[1], self.posn[2])  ## move
        glScale(self.size[0], self.size[1], self.size[2])  ## set size
        glColor3f(self.color[0], self.color[1], self.color[2])  ## set colour
        DrawCube()  ## draw a scaled cube
        glPopMatrix()  ## restore old coord frame

#####################################################
#  #### Bar class, i.e., for dynamic illustration
#####################################################
class Ground:
    color = [0.1, 0.1, 0.7]
    size = [7, 0.02, 7]
    posn = np.array([0.0, -3.5, 0.0])

    def draw(self):
        glPushMatrix()  ## save copy of coord frame
        glTranslatef(self.posn[0], self.posn[1], self.posn[2])  ## move
        glScale(self.size[0], self.size[1], self.size[2])  ## set size
        glColor3f(self.color[0], self.color[1], self.color[2])  ## set colour
        DrawCube()  ## draw a scaled cube
        glPopMatrix()  ## restore old coord frame


#####################################################
#### main():   launches app
#####################################################

def main():
    global window
    global link1, link2, link3, link4, alone
    global kBar, pBar
    global ground
    global multi

    multi = 0

    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH)  # display mode
    glutInitWindowSize(640, 480)  # window size
    glutInitWindowPosition(0, 0)  # window coords for mouse start at top-left
    window = glutCreateWindow("Lotus' Pendulum")
    glutDisplayFunc(DrawWorld)  # register the function to draw the world
    # glutFullScreen()               # full screen
    glutIdleFunc(SimWorld)  # when doing nothing, redraw the scene
    glutReshapeFunc(ReSizeGLScene)  # register the function to call when window is resized
    glutKeyboardFunc(keyPressed)  # register the function to call when keyboard is pressed
    InitGL(640, 480)  # initialize window

    link1 = Link();
    link2 = Link();
    link3 = Link();
    link4 = Link();

    alone = Link();

    kBar = Bar()
    pBar = Bar()

    ground = Ground();

    resetSim()

    glutMainLoop()  # start event processing loop


#####################################################
#### keyPressed():  called whenever a key is pressed
#####################################################

def resetSim():
    global link1, link2, link3, link4
    global kBar, pBar
    global simTime, simRun
    global ground
    global boardY

    printf("Simulation reset\n")
    simRun = True
    simTime = 0

    alone.size = [0.04, 1.0, 0.12]
    alone.color = [0.5, 0.9, 0.9]
    alone.vel = np.array([0.0, 0.0, 0.0])
    alone.theta = math.pi * 3/4
    alone.posn = np.array([-0.5 * math.sin(alone.theta), 0.5 * math.cos(alone.theta), 0.0])
    alone.omega = 0.0  ## radians per second

    link1.size = [0.04, 1.0, 0.12]
    link1.color = [0.85, 0.9, 1.0]
    link1.vel = np.array([0.0, 0.0, 0.0])
    link1.theta = math.pi * 3/4
    link1.posn = [0,0,0] + np.array([-0.5 * math.sin(link1.theta), 0.5 * math.cos(link1.theta), 0.0])
    link1.omega = 0.0  ## radians per second

    link2.size = [0.04, 1.0, 0.12]
    link2.color = [0.65, 0.79, 0.99]
    link2.theta = math.pi * 3/4
    link2.vel = np.array([0.0, 0.0, 0.0])
    link2.posn = link1.posn + np.array([-math.sin(link1.theta), math.cos(link1.theta), 0.0])
    link2.omega = 0.0  ## radians per second

    link3.size = [0.04, 1.0, 0.12]
    link3.color = [0.45, 0.66, 0.98]
    link3.theta = math.pi * 3/4
    link3.posn = link2.posn + np.array([-math.sin(link2.theta), math.cos(link2.theta), 0.0])
    link3.vel = np.array([0.0, 0.0, 0.0])
    link3.omega = 0.0  ## radians per second

    link4.size = [0.04, 1.0, 0.12]
    link4.color = [0.3, 0.55, 0.95]
    link4.theta = math.pi * 3 / 4
    link4.posn = link3.posn + np.array([-math.sin(link3.theta), math.cos(link3.theta), 0.0])
    link4.vel = np.array([0.0, 0.0, 0.0])
    link4.omega = 0.0  ## radians per second

    kBar.color = [0.6, 0.7, 0.8]
    kBar.size = [0.2, 0.4, 0.2]
    kBar.posn = np.array([1.0, -3.5, 6.5])

    pBar.color = [0.6, 0.9, 0.65]
    pBar.size = [0.2, 0.4, 0.2]
    pBar.posn = np.array([kBar.posn[0], kBar.size[1], kBar.posn[2]])

    ground.color = [0.3, 0.3, 0.7]
    boardY = ground.posn[1]

#####################################################
#### keyPressed():  called whenever a key is pressed
#####################################################

def keyPressed(key, x, y):
    global simRun
    global multi
    global boardY

    ch = key.decode("utf-8")
    if ch == ' ':  #### toggle the simulation
        if (simRun == True):
            simRun = False
        else:
            simRun = True
    elif ch == chr(27):  #### ESC key
        sys.exit()
    elif ch == 'q':  #### quit
        sys.exit()
    elif ch == 'r':  #### reset simulation
        resetSim()
    elif ch == 'n':  #### reset simulation
        multi = 0
        resetSim()
    elif ch == 'm':  #### reset simulation
        multi = 1
        resetSim()
    elif ch == 'w':  #### reset simulation
        boardY += 1
    elif ch == 's':  #### reset simulation
        boardY -= 1



#####################################################
#### SimWorld():  simulates a time step
#####################################################

def SimWorld():
    global simTime, dT, simRun
    global link1, link2, link3, link4, alone
    global boardY

    I = 1.0/12.0

    deltaTheta = 2.4
    if (simRun == False):  ## is simulation stopped?
        return

###############################
#### For 4-link
###############################

    r1 = [0.5 * math.sin(link1.theta),-0.5 * math.cos(link1.theta),0]
    r1x = r1[0]
    r1y = r1[1]

    r2 = [0.5 * math.sin(link2.theta),-0.5 * math.cos(link2.theta),0]
    r2x = r2[0]
    r2y = r2[1]

    r3 = [0.5 * math.sin(link3.theta),-0.5 * math.cos(link3.theta),0]
    r3x = r3[0]
    r3y = r3[1]

    r4 = [0.5 * math.sin(link4.theta),-0.5 * math.cos(link4.theta),0]
    r4x = r4[0]
    r4y = r4[1]

    constraint1 =1*(link1.posn + r1 - [0,0,0])+ 0.01*(np.cross([0,0,link1.omega], r1) - link1.vel)
    constraint2 =1*(link2.posn + r2 - (link1.posn - r1))+ 0.01*(np.cross([0,0,link2.omega], r2) - link2.vel)
    constraint3 =1*(link3.posn + r3 - (link2.posn - r2))+ 0.01*(np.cross([0,0,link3.omega], r3) - link3.vel)
    constraint4 =1*(link4.posn + r4 - (link3.posn - r3))+ 0.01*(np.cross([0,0,link4.omega], r4) - link4.vel)

    t1x = -r1x * link1.omega * link1.omega + constraint1[0]
    t1y = -r1y * link1.omega * link1.omega + constraint1[1]

    t2x = r1x * link1.omega * link1.omega + r2x * link2.omega * link2.omega - constraint2[0]
    t2y = r1y * link1.omega * link1.omega + r2y * link2.omega * link2.omega - constraint2[1]

    t3x = r2x * link2.omega * link2.omega + r3x * link3.omega * link3.omega - constraint3[0]
    t3y = r2y * link2.omega * link2.omega + r3y * link3.omega * link3.omega - constraint3[1]

    t4x = r3x * link3.omega * link3.omega + r4x * link4.omega * link4.omega - constraint4[0]
    t4y = r3y * link3.omega * link3.omega + r4y * link4.omega * link4.omega - constraint4[1]

    damp1 = link1.omega * 0.5
    damp2 = link2.omega * 0.5
    damp3 = link3.omega * 0.5
    damp4 = link4.omega * 0.5

    ground.posn[1] = boardY
    bottom = link4.posn - r4 - 0.2

    penalty = 80.0*(ground.posn[1] - bottom[1]) - 20.0* link4.vel[1]
    if ground.posn[1] > bottom[1]:
        ground.color = [0.9,0.9,0.9]
    else:
        ground.color = [0.85, 0.85, 0.85]
        penalty = 0
    # print(penalty)

    #                1      2      3      4        A      B       C        D
    a = np.array([[1,0,0, 0,0,0, 0,0,0, 0,0,0,   -1,0,   1,0,    0,0,    0,0], #1
                  [0,1,0, 0,0,0, 0,0,0, 0,0,0,   0,-1,   0,1,    0,0,    0,0],
                  [0,0,I, 0,0,0, 0,0,0, 0,0,0, r1y,-r1x,r1y,-r1x,0,0,    0,0],
                  [0,0,0, 1,0,0, 0,0,0, 0,0,0,    0,0,  -1,0,    1,0,    0,0], #2
                  [0,0,0, 0,1,0, 0,0,0, 0,0,0,    0,0,  0,-1,    0,1,    0,0],
                  [0,0,0, 0,0,I, 0,0,0, 0,0,0,    0,0,r2y,-r2x,r2y,-r2x, 0,0],
                  [0,0,0, 0,0,0, 1,0,0, 0,0,0,    0,0,   0,0,   -1,0,    1,0], #3
                  [0,0,0, 0,0,0, 0,1,0, 0,0,0,    0,0,   0,0,   0,-1,    0,1],
                  [0,0,0, 0,0,0, 0,0,I, 0,0,0,    0,0,   0,0,r3y,-r3x,r3y,-r3x],
                  [0,0,0, 0,0,0, 0,0,0, 1,0,0,    0,0,   0,0,    0,0,   -1,0], #4
                  [0,0,0, 0,0,0, 0,0,0, 0,1,0,    0,0,   0,0,    0,0,   0,-1],
                  [0,0,0, 0,0,0, 0,0,0, 0,0,I,    0,0,   0,0,    0,0,r4y,-r4x],

                  [-1,0,r1y,  0,0,0,  0,0,0,  0,0,0,    0,0,   0,0,    0,0,    0,0], #A
                  [0,-1,-r1x, 0,0,0,  0,0,0,  0,0,0,    0,0,   0,0,    0,0,    0,0],
                  [-1,0,-r1y, 1,0,-r2y,0,0,0, 0,0,0,    0,0,   0,0,    0,0,    0,0], #B
                  [0,-1,r1x,  0,1,r2x, 0,0,0, 0,0,0,    0,0,   0,0,    0,0,    0,0],
                  [0,0,0,   -1,0,-r2y,1,0,-r3y,0,0,0,   0,0,   0,0,    0,0,    0,0], #C
                  [0,0,0,   0,-1,r2x,0,1,r3x, 0,0,0,    0,0,   0,0,    0,0,    0,0],
                  [0,0,0,   0,0,0,  -1,0,-r3y,1,0,-r4y, 0,0,   0,0,    0,0,    0,0], #D
                  [0,0,0,   0,0,0,  0,-1,r3x, 0,1,r4x,  0,0,   0,0,    0,0,    0,0]])

    b = np.array([0,-10,-damp1, 0,-10,-damp2, 0,-10,-damp3, 0,-10+penalty, -damp4 - r4x * penalty, t1x,t1y, t2x,t2y, t3x,t3y, t4x,t4y])

    x = np.linalg.solve(a, b)

        #### solve for the equations of motion (simple in this case!)
    acc1 = np.array([x[0], x[1], 0])  ### linear acceleration = [0, -G, 0]
    acc2 = np.array([x[3], x[4], 0])  ### linear acceleration = [0, -G, 0]
    acc3 = np.array([x[6], x[7], 0])  ### linear acceleration = [0, -G, 0]
    acc4 = np.array([x[9], x[10], 0])  ### linear acceleration = [0, -G, 0]

    omega_dot1 = x[2]  ### assume no angular acceleration
    omega_dot2 = x[5]  ### assume no angular acceleration
    omega_dot3 = x[8]  ### assume no angular acceleration
    omega_dot4 = x[11]  ### assume no angular acceleration

    #### explicit Euler integration to update the state
    link1.posn += link1.vel * dT
    link1.vel += acc1 * dT
    link1.theta += link1.omega * dT
    link1.omega += omega_dot1 * dT

    link2.posn += link2.vel * dT
    link2.vel += acc2 * dT
    link2.theta += link2.omega * dT
    link2.omega += omega_dot2 * dT

    link3.posn += link3.vel * dT
    link3.vel += acc3 * dT
    link3.theta += link3.omega * dT
    link3.omega += omega_dot3 * dT

    link4.posn += link4.vel * dT
    link4.vel += acc4 * dT
    link4.theta += link4.omega * dT
    link4.omega += omega_dot4 * dT

#################################
#######   For Alone:
#################################

    ra = [0.5 * math.sin(alone.theta), -0.5 * math.cos(alone.theta), 0]
    rax = ra[0]
    ray = ra[1]

    dampA = alone.omega * 0.03

    constraintA = 1 * (alone.posn + ra - [0, 0, 0]) + 0.05 * (np.cross([0, 0, alone.omega], ra) - alone.vel)
    tax = -rax * alone.omega * alone.omega + constraintA[0]
    tay = -ray * alone.omega * alone.omega + constraintA[1]

    aa = np.array([[1, 0, 0, -1, 0],
                   [0, 1, 0, 0, -1],
                   [0, 0, I, -ray, rax],
                   [-1, 0, ray, 0, 0],
                   [0, -1, -rax, 0, 0]])
    ab = np.array([0, -10, dampA, tax, tay])
    ax = np.linalg.solve(aa, ab)

    accA = np.array([ax[0], ax[1], 0])  ### linear acceleration = [0, -G, 0]
    omega_dotA = ax[2]  ### assume no angular acceleration

    #### explicit Euler integration to update the state
    alone.posn += alone.vel * dT
    alone.vel += accA * dT
    alone.theta += alone.omega * dT
    alone.omega += omega_dotA * dT

    # PE = mgh
    pe = 10*(alone.posn[1] +0.55)
    # KE = 0.5(I + md^2)w^2
    ke = 0.5*(0.25/3.0 + 0.25)*alone.omega*alone.omega

    pBar.size[1] = 0.1* pe
    pBar.posn[1] = pBar.size[1] * 0.5 -2.8

    kBar.size[1] = 0.1* ke
    kBar.posn[1] = pBar.size[1] + kBar.size[1] * 0.5 -2.8

    simTime += dT

    #### draw the updated state
    DrawWorld()
    printf("simTime=%.2f\n", simTime)


#####################################################
#### DrawWorld():  draw the world
#####################################################

def DrawWorld():
    global link1, link2, link3,link4

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  # Clear The Screen And The Depth Buffer
    glLoadIdentity();
    gluLookAt(2.5, -2.0, 9.0,   0, -2.5, 0,    0, 1, 0)

    DrawOrigin()

    if multi == 0:
        alone.draw()
        kBar.draw()
        pBar.draw()
    else:
        link1.draw()
        link2.draw()
        link3.draw()
        link4.draw()
        ground.draw()

    glutSwapBuffers()  # swap the buffers to display what was just drawn


#####################################################
#### initGL():  does standard OpenGL initialization work
#####################################################

def InitGL(Width, Height):  # We call this right after our OpenGL window is created.
    glClearColor(1.0, 1.0, 1.0, 0.0)  # This Will Clear The Background Color To Black
    glClearDepth(1.0)  # Enables Clearing Of The Depth Buffer
    glDepthFunc(GL_LESS)  # The Type Of Depth Test To Do
    glEnable(GL_DEPTH_TEST)  # Enables Depth Testing
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glShadeModel(GL_SMOOTH)  # Enables Smooth Color Shading
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()  # Reset The Projection Matrix
    gluPerspective(45.0, float(Width) / float(Height), 0.1, 100.0)
    glMatrixMode(GL_MODELVIEW)


#####################################################
#### ReSizeGLScene():    called when window is resized
#####################################################

def ReSizeGLScene(Width, Height):
    if Height == 0:  # Prevent A Divide By Zero If The Window Is Too Small
        Height = 1
    glViewport(0, 0, Width, Height)  # Reset The Current Viewport And Perspective Transformation
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45.0, float(Width) / float(Height), 0.1,
                   100.0)  ## 45 deg horizontal field of view, aspect ratio, near, far
    glMatrixMode(GL_MODELVIEW)


#####################################################
#### DrawOrigin():  draws RGB lines for XYZ origin of coordinate system
#####################################################

def DrawOrigin():
    glLineWidth(3.0);

    glColor3f(0.65, 0.65, 0.65)  ## light red x-axis
    glBegin(GL_LINES)
    glVertex3f(0, 0, 0)
    glVertex3f(1, 0, 0)
    glEnd()

    glColor3f(0.5, 0.5, 0.5)  ## light green y-axis
    glBegin(GL_LINES)
    glVertex3f(0, 0, 0)
    glVertex3f(0, 1, 0)
    glEnd()

    glColor3f(0.8, 0.8, 0.8)  ## light blue z-axis
    glBegin(GL_LINES)
    glVertex3f(0, 0, 0)
    glVertex3f(0, 0, 1)
    glEnd()


#####################################################
#### DrawCube():  draws a cube that spans from (-1,-1,-1) to (1,1,1)
#####################################################

def DrawCube():
    glScalef(0.5, 0.5, 0.5);  # dimensions below are for a 2x2x2 cube, so scale it down by a half first
    glBegin(GL_QUADS);  # Start Drawing The Cube

    glVertex3f(1.0, 1.0, -1.0);  # Top Right Of The Quad (Top)
    glVertex3f(-1.0, 1.0, -1.0);  # Top Left Of The Quad (Top)
    glVertex3f(-1.0, 1.0, 1.0);  # Bottom Left Of The Quad (Top)
    glVertex3f(1.0, 1.0, 1.0);  # Bottom Right Of The Quad (Top)

    glVertex3f(1.0, -1.0, 1.0);  # Top Right Of The Quad (Bottom)
    glVertex3f(-1.0, -1.0, 1.0);  # Top Left Of The Quad (Bottom)
    glVertex3f(-1.0, -1.0, -1.0);  # Bottom Left Of The Quad (Bottom)
    glVertex3f(1.0, -1.0, -1.0);  # Bottom Right Of The Quad (Bottom)

    glVertex3f(1.0, 1.0, 1.0);  # Top Right Of The Quad (Front)
    glVertex3f(-1.0, 1.0, 1.0);  # Top Left Of The Quad (Front)
    glVertex3f(-1.0, -1.0, 1.0);  # Bottom Left Of The Quad (Front)
    glVertex3f(1.0, -1.0, 1.0);  # Bottom Right Of The Quad (Front)

    glVertex3f(1.0, -1.0, -1.0);  # Bottom Left Of The Quad (Back)
    glVertex3f(-1.0, -1.0, -1.0);  # Bottom Right Of The Quad (Back)
    glVertex3f(-1.0, 1.0, -1.0);  # Top Right Of The Quad (Back)
    glVertex3f(1.0, 1.0, -1.0);  # Top Left Of The Quad (Back)

    glVertex3f(-1.0, 1.0, 1.0);  # Top Right Of The Quad (Left)
    glVertex3f(-1.0, 1.0, -1.0);  # Top Left Of The Quad (Left)
    glVertex3f(-1.0, -1.0, -1.0);  # Bottom Left Of The Quad (Left)
    glVertex3f(-1.0, -1.0, 1.0);  # Bottom Right Of The Quad (Left)

    glVertex3f(1.0, 1.0, -1.0);  # Top Right Of The Quad (Right)
    glVertex3f(1.0, 1.0, 1.0);  # Top Left Of The Quad (Right)
    glVertex3f(1.0, -1.0, 1.0);  # Bottom Left Of The Quad (Right)
    glVertex3f(1.0, -1.0, -1.0);  # Bottom Right Of The Quad (Right)
    glEnd();  # Done Drawing The Quad

    ### Draw the wireframe edges
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(1.0);

    glBegin(GL_LINE_LOOP);
    glVertex3f(1.0, 1.0, -1.0);  # Top Right Of The Quad (Top)
    glVertex3f(-1.0, 1.0, -1.0);  # Top Left Of The Quad (Top)
    glVertex3f(-1.0, 1.0, 1.0);  # Bottom Left Of The Quad (Top)
    glVertex3f(1.0, 1.0, 1.0);  # Bottom Right Of The Quad (Top)
    glEnd();  # Done Drawing The Quad

    glBegin(GL_LINE_LOOP);
    glVertex3f(1.0, -1.0, 1.0);  # Top Right Of The Quad (Bottom)
    glVertex3f(-1.0, -1.0, 1.0);  # Top Left Of The Quad (Bottom)
    glVertex3f(-1.0, -1.0, -1.0);  # Bottom Left Of The Quad (Bottom)
    glVertex3f(1.0, -1.0, -1.0);  # Bottom Right Of The Quad (Bottom)
    glEnd();  # Done Drawing The Quad

    glBegin(GL_LINE_LOOP);
    glVertex3f(1.0, 1.0, 1.0);  # Top Right Of The Quad (Front)
    glVertex3f(-1.0, 1.0, 1.0);  # Top Left Of The Quad (Front)
    glVertex3f(-1.0, -1.0, 1.0);  # Bottom Left Of The Quad (Front)
    glVertex3f(1.0, -1.0, 1.0);  # Bottom Right Of The Quad (Front)
    glEnd();  # Done Drawing The Quad

    glBegin(GL_LINE_LOOP);
    glVertex3f(1.0, -1.0, -1.0);  # Bottom Left Of The Quad (Back)
    glVertex3f(-1.0, -1.0, -1.0);  # Bottom Right Of The Quad (Back)
    glVertex3f(-1.0, 1.0, -1.0);  # Top Right Of The Quad (Back)
    glVertex3f(1.0, 1.0, -1.0);  # Top Left Of The Quad (Back)
    glEnd();  # Done Drawing The Quad

    glBegin(GL_LINE_LOOP);
    glVertex3f(-1.0, 1.0, 1.0);  # Top Right Of The Quad (Left)
    glVertex3f(-1.0, 1.0, -1.0);  # Top Left Of The Quad (Left)
    glVertex3f(-1.0, -1.0, -1.0);  # Bottom Left Of The Quad (Left)
    glVertex3f(-1.0, -1.0, 1.0);  # Bottom Right Of The Quad (Left)
    glEnd();  # Done Drawing The Quad

    glBegin(GL_LINE_LOOP);
    glVertex3f(1.0, 1.0, -1.0);  # Top Right Of The Quad (Right)
    glVertex3f(1.0, 1.0, 1.0);  # Top Left Of The Quad (Right)
    glVertex3f(1.0, -1.0, 1.0);  # Bottom Left Of The Quad (Right)
    glVertex3f(1.0, -1.0, -1.0);  # Bottom Right Of The Quad (Right)
    glEnd();  # Done Drawing The Quad


####################################################
# printf()
####################################################

def printf(format, *args):
    sys.stdout.write(format % args)


################################################################################
# start the app

print("Hit ESC key to quit.")
main()
