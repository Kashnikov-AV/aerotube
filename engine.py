import numpy as np
#import time
#import matplotlib.pyplot
#import matplotlib.animation
import pyqtgraph as pg
from dataclasses import dataclass

@dataclass
class Initials:
    # Define constants:
    height = 80  # lattice dimensions
    width = 200
    viscosity = 0.02  # fluid viscosity
    omega = 1 / (3 * viscosity + 0.5)  # "relaxation" parameter
    u0 = 0.1  # initial and in-flow speed
    four9ths = 4.0 / 9.0  # abbreviations for lattice-Boltzmann weight factors
    one9th = 1.0 / 9.0
    one36th = 1.0 / 36.0
    performanceData = False  # set to True if performance data is desired

@dataclass
class Data:
    # init data
    initials = Initials()
    # Initialize all the arrays to steady rightward flow:
    n0 = initials.four9ths * (np.ones((initials.height, initials.width)) - 1.5 * initials.u0 ** 2)  # particle densities along 9 directions
    nN = initials.one9th * (np.ones((initials.height, initials.width)) - 1.5 * initials.u0 ** 2)
    nS = initials.one9th * (np.ones((initials.height, initials.width)) - 1.5 * initials.u0 ** 2)
    nE = initials.one9th * (np.ones((initials.height, initials.width)) + 3 * initials.u0 + 4.5 * initials.u0 ** 2 - 1.5 * initials.u0 ** 2)
    nW = initials.one9th * (np.ones((initials.height, initials.width)) - 3 * initials.u0 + 4.5 * initials.u0 ** 2 - 1.5 * initials.u0 ** 2)
    nNE = initials.one36th * (np.ones((initials.height, initials.width)) + 3 * initials.u0 + 4.5 * initials.u0 ** 2 - 1.5 * initials.u0 ** 2)
    nSE = initials.one36th * (np.ones((initials.height, initials.width)) + 3 * initials.u0 + 4.5 * initials.u0 ** 2 - 1.5 * initials.u0 ** 2)
    nNW = initials.one36th * (np.ones((initials.height, initials.width)) - 3 * initials.u0 + 4.5 * initials.u0 ** 2 - 1.5 * initials.u0 ** 2)
    nSW = initials.one36th * (np.ones((initials.height, initials.width)) - 3 * initials.u0 + 4.5 * initials.u0 ** 2 - 1.5 * initials.u0 ** 2)
    # Initialize barriers:
    barrier = np.zeros((initials.width, initials.height), bool)  # True wherever there's a barrier
    
    @property
    def rho(self):
        return self.n0 + self.nN + self.nS + self.nE + self.nW + self.nNE + self.nSE + self.nNW + self.nSW  # macroscopic density

    @property
    def ux(self):
        return (self.nE + self.nNE + self.nSE - self.nW - self.nNW - self.nSW) / self.rho  # macroscopic x velocity

    @property
    def uy(self):
        return(self.nN + self.nNE + self.nNW - self.nS - self.nSE - self.nSW) / self.rho  # macroscopic y velocity

    def update_U0(self, data):
        self.initials.U0 = data / 1000

    def update_viscosity(self, data):
        self.initials.viscosity = data / 1000
    
    def set_barrier(self):
        self.barrier[90:111, 29] = True
        self.barrier[90:111, 49]  = True
        self.barrier[90, 29:49] = True
        self.barrier[110, 29:49] = True
        # self.barrier[int((self.initials.height / 2) - 8) : int((self.initials.height / 2) + 8, self.initials.height / 2)] = True  # simple linear barrier
        self.barrierN  = np.roll(self.barrier,   1, axis=0)  # sites just north of barriers
        self.barrierS  = np.roll(self.barrier,  -1, axis=0)  # sites just south of barriers
        self.barrierE  = np.roll(self.barrier,   1, axis=1)  # etc.
        self.barrierW  = np.roll(self.barrier,  -1, axis=1)
        self.barrierNE = np.roll(self.barrierN,  1, axis=1)
        self.barrierNW = np.roll(self.barrierN, -1, axis=1)
        self.barrierSE = np.roll(self.barrierS,  1, axis=1)
        self.barrierSW = np.roll(self.barrierS, -1, axis=1)
        
    def draw_barrier(self):
        bImageArray = np.zeros((self.initials.width, self.initials.height, 4), np.uint8)  # an RGBA image
        bImageArray[self.barrier, 3] = 255  # set alpha=255 only at barrier sites
        barrierImage = pg.ImageItem(self.barrier)
        return barrierImage
    
    # Move all particles by one step along their directions of motion (pbc):
    def stream(self):
        self.nN  = np.roll(self.nN, 1, axis=0)  # axis 0 is north-south; + direction is north
        self.nNE = np.roll(self.nNE, 1, axis=0)
        self.nNW = np.roll(self.nNW, 1, axis=0)
        self.nS  = np.roll(self.nS, -1, axis=0)
        self.nSE = np.roll(self.nSE, -1, axis=0)
        self.nSW = np.roll(self.nSW, -1, axis=0)
        self.nE  = np.roll(self.nE, 1, axis=1)  # axis 1 is east-west; + direction is east
        self.nNE = np.roll(self.nNE, 1, axis=1)
        self.nSE = np.roll(self.nSE, 1, axis=1)
        self.nW  = np.roll(self.nW, -1, axis=1)
        self.nNW = np.roll(self.nNW, -1, axis=1)
        self.nSW = np.roll(self.nSW, -1, axis=1)
        # Use tricky boolean arrays to handle barrier collisions (bounce-back):
        self.nN[self.barrierN] = self.nS[self.barrier]
        self.nS[self.barrierS] = self.nN[self.barrier]
        self.nE[self.barrierE] = self.nW[self.barrier]
        self.nW[self.barrierW] = self.nE[self.barrier]
        self.nNE[self.barrierNE] = self.nSW[self.barrier]
        self.nNW[self.barrierNW] = self.nSE[self.barrier]
        self.nSE[self.barrierSE] = self.nNW[self.barrier]
        self.nSW[self.barrierSW] = self.nNE[self.barrier]




# def circle_mask(center_x, center_y, radius, array_shape):
#     height, width = array_shape
#     y, x = np.ogrid[:height, :width]
#
#     # Вычисляем расстояние до центра
#     dist = np.sqrt((x - center_x)**2 + (y - center_y)**2)
#
#     # Окружность: точки на расстоянии ≈ radius
#     # Используем толщину в 1 пиксель
#     thickness = 1
#     return (dist >= radius - thickness/2) & (dist <= radius + thickness/2)
#
# barrier = circle_mask(50, 50, 30, (height, width))

    # Collide particles within each cell to redistribute velocities 
    # (could be optimized a little more):
    def collide(self):
        ux2 = self.ux * self.ux  # pre-compute terms used repeatedly...
        uy2 = self.uy * self.uy
        u2 = ux2 + uy2
        omu215 = 1 - 1.5 * u2  # "one minus u2 times 1.5"
        uxuy = self.ux * self.uy
        self.n0 = (1 - self.initials.omega) * n0 + omega * four9ths * rho * omu215
        self.nN = (1 - self.initials.omega) * nN + omega * one9th * rho * (omu215 + 3 * uy + 4.5 * uy2)
        self.nS = (1 - self.initials.omega) * nS + omega * one9th * rho * (omu215 - 3 * uy + 4.5 * uy2)
        self.nE = (1 - self.initials.omega) * nE + omega * one9th * rho * (omu215 + 3 * ux + 4.5 * ux2)
        self.nW = (1 - self.initials.omega) * nW + omega * one9th * rho * (omu215 - 3 * ux + 4.5 * ux2)
        self.nNE = (1 - omega) * nNE + omega * one36th * rho * (omu215 + 3 * (ux + uy) + 4.5 * (u2 + 2 * uxuy))
        self.nNW = (1 - omega) * nNW + omega * one36th * rho * (omu215 + 3 * (-ux + uy) + 4.5 * (u2 - 2 * uxuy))
        self.nSE = (1 - omega) * nSE + omega * one36th * rho * (omu215 + 3 * (ux - uy) + 4.5 * (u2 - 2 * uxuy))
        self.nSW = (1 - omega) * nSW + omega * one36th * rho * (omu215 + 3 * (-ux - uy) + 4.5 * (u2 + 2 * uxuy))
        # Force steady rightward flow at ends (no need to set 0, N, and S components):
        self.nE[:, 0] = one9th * (1 + 3 * u0 + 4.5 * u0 ** 2 - 1.5 * u0 ** 2)
        self.nW[:, 0] = one9th * (1 - 3 * u0 + 4.5 * u0 ** 2 - 1.5 * u0 ** 2)
        self.nNE[:, 0] = one36th * (1 + 3 * u0 + 4.5 * u0 ** 2 - 1.5 * u0 ** 2)
        self.nSE[:, 0] = one36th * (1 + 3 * u0 + 4.5 * u0 ** 2 - 1.5 * u0 ** 2)
        self.nNW[:, 0] = one36th * (1 - 3 * u0 + 4.5 * u0 ** 2 - 1.5 * u0 ** 2)
        self.nSW[:, 0] = one36th * (1 - 3 * u0 + 4.5 * u0 ** 2 - 1.5 * u0 ** 2)

# # Compute curl of the macroscopic velocity field:
# def curl(ux, uy):
#     return np.roll(uy, -1, axis=1) - np.roll(uy, 1, axis=1) - np.roll(ux, -1, axis=0) + np.roll(ux, 1,
#                                                                                                             axis=0)
#
#
# # Here comes the graphics and animation...
# theFig = matplotlib.pyplot.figure(figsize=(10, 4))
# # FIXME
# fluidImage = matplotlib.pyplot.imshow(curl(ux, uy), origin='lower', norm=matplotlib.pyplot.Normalize(-.1, .1),
#                                       cmap=matplotlib.pyplot.get_cmap(cmaps[8]), interpolation='none')
# # See http://www.loria.fr/~rougier/teaching/matplotlib/#colormaps for other cmap options
# bImageArray = np.zeros((height, width, 4), np.uint8)  # an RGBA image
# bImageArray[barrier, 3] = 255  # set alpha=255 only at barrier sites
# barrierImage = matplotlib.pyplot.imshow(bImageArray, origin='lower', interpolation='none')
#
# # Function called for each successive animation frame:
# startTime = time.perf_counter()
#
#
# # frameList = open('frameList.txt','w')		# file containing list of images (to make movie)
# def nextFrame(arg):  # (arg is the frame number, which we don't need)
#     global startTime
#     if performanceData and (arg % 100 == 0) and (arg > 0):
#         endTime = time.perf_counter()
#         print
#         "%1.1f" % (100 / (endTime - startTime)), 'frames per second'
#         startTime = endTime
#     # frameName = "frame%04d.png" % arg
#     # matplotlib.pyplot.savefig(frameName)
#     # frameList.write(frameName + '\n')
#     for step in range(20):  # adjust number of steps for smooth animation
#         stream()
#         collide()
#     fluidImage.set_array(curl(ux, uy))
#     return (fluidImage, barrierImage)  # return the figure elements to redraw
#
#
# animate = matplotlib.animation.FuncAnimation(theFig, nextFrame, interval=1, blit=True)
# matplotlib.pyplot.show()
