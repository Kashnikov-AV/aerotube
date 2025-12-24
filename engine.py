import numpy as np
import pyqtgraph as pg
from dataclasses import dataclass

@dataclass
class Initials:
    height = 80  # решетка
    width = 200
    viscosity = 0.02  # вязкость
    omega = 1 / (3 * viscosity + 0.5)  # "relaxation" parameter
    u0 = 0.002  # initial and in-flow speed
    four9ths = 4.0 / 9.0  # abbreviations for lattice-Boltzmann weight factors
    one9th = 1.0 / 9.0
    one36th = 1.0 / 36.0

@dataclass
class Data:
    # init data
    ini = Initials()
    # Initialize all the arrays to steady rightward flow:
    n0 = ini.four9ths * (np.ones((ini.width, ini.height)) - 1.5 * ini.u0 ** 2)
    nN = ini.one9th *   (np.ones((ini.width, ini.height)) - 1.5 * ini.u0 ** 2)
    nS = ini.one9th *   (np.ones((ini.width, ini.height)) - 1.5 * ini.u0 ** 2)
    nE = ini.one9th *   (np.ones((ini.width, ini.height)) + 3 * ini.u0 + 4.5 * ini.u0 ** 2 - 1.5 * ini.u0 ** 2)
    nW = ini.one9th *   (np.ones((ini.width, ini.height)) - 3 * ini.u0 + 4.5 * ini.u0 ** 2 - 1.5 * ini.u0 ** 2)
    nNE = ini.one36th * (np.ones((ini.width, ini.height)) + 3 * ini.u0 + 4.5 * ini.u0 ** 2 - 1.5 * ini.u0 ** 2)
    nSE = ini.one36th * (np.ones((ini.width, ini.height)) + 3 * ini.u0 + 4.5 * ini.u0 ** 2 - 1.5 * ini.u0 ** 2)
    nNW = ini.one36th * (np.ones((ini.width, ini.height)) - 3 * ini.u0 + 4.5 * ini.u0 ** 2 - 1.5 * ini.u0 ** 2)
    nSW = ini.one36th * (np.ones((ini.width, ini.height)) - 3 * ini.u0 + 4.5 * ini.u0 ** 2 - 1.5 * ini.u0 ** 2)
    # инициализация барьера
    barrier = np.zeros((ini.width, ini.height), bool)  # True баррьер
    barrierN  = np.zeros((ini.width, ini.height), bool)
    barrierS  = np.zeros((ini.width, ini.height), bool)
    barrierE  = np.zeros((ini.width, ini.height), bool)
    barrierW  = np.zeros((ini.width, ini.height), bool)
    barrierNE = np.zeros((ini.width, ini.height), bool)
    barrierNW = np.zeros((ini.width, ini.height), bool)
    barrierSE = np.zeros((ini.width, ini.height), bool)
    barrierSW = np.zeros((ini.width, ini.height), bool)

    @property
    def rho(self):
        return self.n0 + self.nN + self.nS + self.nE + self.nW + self.nNE + self.nSE + self.nNW + self.nSW  # macroscopic density

    @property
    def ux(self):
        return (self.nE + self.nNE + self.nSE - self.nW - self.nNW - self.nSW) / self.rho  # macroscopic x velocity

    @property
    def uy(self):
        return (self.nN + self.nNE + self.nNW - self.nS - self.nSE - self.nSW) / self.rho  # macroscopic y velocity

    def update_U0(self, data):
        self.ini.U0 = data / 1000

    def update_viscosity(self, data):
        self.ini.viscosity = data / 1000

    def reset_barrier(self):
        self.barrier = np.zeros((self.ini.width, self.ini.height), bool)

    def other_barriers(self):
        self.barrierN  = np.roll(self.barrier,  1, axis=1)  # север
        self.barrierS  = np.roll(self.barrier, -1, axis=1)  # юг
        self.barrierE  = np.roll(self.barrier,  1, axis=0)  # etc.
        self.barrierW  = np.roll(self.barrier, -1, axis=0)
        self.barrierNE = np.roll(self.barrierN, 1, axis=0)
        self.barrierNW = np.roll(self.barrierN,-1, axis=0)
        self.barrierSE = np.roll(self.barrierS, 1, axis=0)
        self.barrierSW = np.roll(self.barrierS,-1, axis=0)

    def set_simple_barrier(self):
        ''' простой барьер '''
        self.reset_barrier()
        self.barrier[60, 29:49] = True
        self.other_barriers()

    def set_square_barrier(self):
        ''' квадратный барьер '''
        self.reset_barrier()
        self.barrier[40:61, 29] = True
        self.barrier[40:61, 49] = True
        self.barrier[40, 29:49] = True
        self.barrier[60, 29:49] = True
        self.other_barriers()

    def set_circle_barrier(self):
        ''' цилиндрический барьер '''
        self.reset_barrier()
        x, y = np.ogrid[:self.ini.width, :self.ini.height]
        dist_from_center = np.sqrt((x - 50)**2 + (y - 39)**2)
        # Булева маска: True внутри круга
        mask = (dist_from_center > 9) & (dist_from_center < 11)
        self.barrier = mask
        self.other_barriers()
        
    def draw_barrier(self):
        bImageArray = np.zeros((self.ini.width, self.ini.height, 4), np.uint8)  # an RGBA image
        bImageArray[self.barrier] = [0, 0, 0, 255]
        bImageArray[~self.barrier] = [1, 1, 1, 0]
        barrierImage = pg.ImageItem(bImageArray)
        return barrierImage
    
    # перемещение частиц вдоль их направления движения
    def stream(self):
        # Сохраняем СТАРЫЕ значения для bounce-back
        nN_old = self.nN.copy()
        nS_old = self.nS.copy()
        nE_old = self.nE.copy()
        nW_old = self.nW.copy()
        nNE_old = self.nNE.copy()
        nNW_old = self.nNW.copy()
        nSE_old = self.nSE.copy()
        nSW_old = self.nSW.copy()

        # Распространение (правильно для вашей конвенции)
        self.nN  = np.roll(self.nN,  1, axis=1)  # север: +Y
        self.nNE = np.roll(self.nNE, 1, axis=1)
        self.nNW = np.roll(self.nNW, 1, axis=1)

        self.nS  = np.roll(self.nS,  -1, axis=1)  # юг: -Y
        self.nSE = np.roll(self.nSE, -1, axis=1)
        self.nSW = np.roll(self.nSW, -1, axis=1)

        self.nE  = np.roll(self.nE,  1, axis=0)  # восток: +X
        self.nNE = np.roll(self.nNE, 1, axis=0)
        self.nSE = np.roll(self.nSE, 1, axis=0)

        self.nW  = np.roll(self.nW,  -1, axis=0)  # запад: -X
        self.nNW = np.roll(self.nNW, -1, axis=0)
        self.nSW = np.roll(self.nSW, -1, axis=0)

        # bounce-back (используем старые значения!)
        self.nN[self.barrier] = nS_old[self.barrier]    # север ← юг
        self.nS[self.barrier] = nN_old[self.barrier]    # юг ← север
        self.nE[self.barrier] = nW_old[self.barrier]    # восток ← запад
        self.nW[self.barrier] = nE_old[self.barrier]    # запад ← восток
        self.nNE[self.barrier] = nSW_old[self.barrier]  # северо-восток ← юго-запад
        self.nNW[self.barrier] = nSE_old[self.barrier]  # северо-запад ← юго-восток
        self.nSE[self.barrier] = nNW_old[self.barrier]  # юго-восток ← северо-запад
        self.nSW[self.barrier] = nNE_old[self.barrier]  # юго-запад ← северо-восток

    def collide(self):
        ux2 = self.ux * self.ux  # pre-compute terms used repeatedly...
        uy2 = self.uy * self.uy
        u2 = ux2 + uy2
        omu215 = 1 - 1.5 * u2  # "one minus u2 times 1.5"
        uxuy = self.ux * self.uy
        self.n0 = (1 - self.ini.omega) * self.n0 + self.ini.omega * self.ini.four9ths * self.rho * omu215
        self.nN = (1 - self.ini.omega) * self.nN + self.ini.omega * self.ini.one9th * self.rho * (omu215 + 3 * self.uy + 4.5 * uy2)
        self.nS = (1 - self.ini.omega) * self.nS + self.ini.omega * self.ini.one9th * self.rho * (omu215 - 3 * self.uy + 4.5 * uy2)
        self.nE = (1 - self.ini.omega) * self.nE + self.ini.omega * self.ini.one9th * self.rho * (omu215 + 3 * self.ux + 4.5 * ux2)
        self.nW = (1 - self.ini.omega) * self.nW + self.ini.omega * self.ini.one9th * self.rho * (omu215 - 3 * self.ux + 4.5 * ux2)
        self.nNE= (1 - self.ini.omega) * self.nNE + self.ini.omega * self.ini.one36th * self.rho * (omu215 + 3 * (self.ux + self.uy) + 4.5 * (u2 + 2 * uxuy))
        self.nNW = (1 - self.ini.omega) * self.nNW + self.ini.omega * self.ini.one36th * self.rho * (omu215 + 3 * (-self.ux + self.uy) + 4.5 * (u2 - 2 * uxuy))
        self.nSE = (1 - self.ini.omega) * self.nSE + self.ini.omega * self.ini.one36th * self.rho * (omu215 + 3 * (self.ux - self.uy) + 4.5 * (u2 - 2 * uxuy))
        self.nSW = (1 - self.ini.omega) * self.nSW + self.ini.omega * self.ini.one36th * self.rho * (omu215 + 3 * (-self.ux - self.uy) + 4.5 * (u2 + 2 * uxuy))
        # Force steady rightward flow at ends (no need to set 0, N, and S components):
        self.nE[ 0, :] = self.ini.one9th  * (1 + 3 * self.ini.u0 + 4.5 * self.ini.u0 ** 2 - 1.5 * self.ini.u0 ** 2)
        self.nW[ 0, :] = self.ini.one9th  * (1 - 3 * self.ini.u0 + 4.5 * self.ini.u0 ** 2 - 1.5 * self.ini.u0 ** 2)
        self.nNE[0, :] = self.ini.one36th * (1 + 3 * self.ini.u0 + 4.5 * self.ini.u0 ** 2 - 1.5 * self.ini.u0 ** 2)
        self.nSE[0, :] = self.ini.one36th * (1 + 3 * self.ini.u0 + 4.5 * self.ini.u0 ** 2 - 1.5 * self.ini.u0 ** 2)
        self.nNW[0, :] = self.ini.one36th * (1 - 3 * self.ini.u0 + 4.5 * self.ini.u0 ** 2 - 1.5 * self.ini.u0 ** 2)
        self.nSW[0, :] = self.ini.one36th * (1 - 3 * self.ini.u0 + 4.5 * self.ini.u0 ** 2 - 1.5 * self.ini.u0 ** 2)
        # ----------------------------------------------------
        self.nE[-1, :] = self.nE[-2, :]
        self.nW[-1, :] = self.nW[-2, :]
        self.nNE[-1, :] = self.nNE[-2, :]
        self.nSE[-1, :] = self.nSE[-2, :]
        self.nNW[-1, :] = self.nNW[-2, :]
        self.nSW[-1, :] = self.nSW[-2, :]
        # ------------------------------------------------------

    # Считаем вихри
    def curl(self):
        dv_dx = np.roll(self.uy, -1, axis=0) - np.roll(self.uy, 1, axis=0)
        # ∂u/∂y = производная ux по Y (axis=1)
        du_dy = np.roll(self.ux, -1, axis=1) - np.roll(self.ux, 1, axis=1)
        # Вихрь: ω = ∂v/∂x - ∂u/∂y
        vorticity = dv_dx - du_dy
        return vorticity
