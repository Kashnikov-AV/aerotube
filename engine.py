import numpy as np
import pyqtgraph as pg
from dataclasses import dataclass

@dataclass
class Initials:
    height = 80  # решетка
    width = 200
    viscosity = 0.005  # вязкость
    omega = 1 / (3 * viscosity + 0.5)  # "relaxation" parameter
    u0 = 0.03  # initial and in-flow speed
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
        self.barrier[80, 29:49] = True
        self.other_barriers()

    def set_square_barrier(self):
        ''' квадратный барьер '''
        self.reset_barrier()
        self.barrier[60:81, 29] = True
        self.barrier[60:81, 49] = True
        self.barrier[60, 29:49] = True
        self.barrier[80, 29:49] = True
        self.other_barriers()

    def set_angle_barrier(self):
        ''' угол барьер '''
        self.reset_barrier()
        self.barrier[60:81, 29] = True
        self.barrier[60, 29:39] = True
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
        # Сохраняем старые значения
        nN_old = self.nN.copy()
        nS_old = self.nS.copy()
        nE_old = self.nE.copy()
        nW_old = self.nW.copy()
        nNE_old = self.nNE.copy()
        nNW_old = self.nNW.copy()
        nSE_old = self.nSE.copy()
        nSW_old = self.nSW.copy()

        # Восток-Запад (X)
        self.nE  = np.roll(self.nE,   1, axis=0)
        self.nNE = np.roll(self.nNE,  1, axis=0)
        self.nSE = np.roll(self.nSE,  1, axis=0)

        self.nW  = np.roll(self.nW,  -1, axis=0)
        self.nNW = np.roll(self.nNW, -1, axis=0)
        self.nSW = np.roll(self.nSW, -1, axis=0)

        # Север-Юг (Y)
        self.nN  = np.roll(self.nN,   1, axis=1)
        self.nNE = np.roll(self.nNE,  1, axis=1)
        self.nNW = np.roll(self.nNW,  1, axis=1)

        self.nS  = np.roll(self.nS,  -1, axis=1)
        self.nSE = np.roll(self.nSE, -1, axis=1)
        self.nSW = np.roll(self.nSW, -1, axis=1)

        # BOUNCE-BACK (используем маски!):
        self.nN[self.barrierN] = nS_old[self.barrier]    # северные ← южные
        self.nS[self.barrierS] = nN_old[self.barrier]    # южные ← северные
        self.nE[self.barrierE] = nW_old[self.barrier]    # восточные ← западные
        self.nW[self.barrierW] = nE_old[self.barrier]    # западные ← восточные
        self.nNE[self.barrierNE] = nSW_old[self.barrier] # северо-восточные ← юго-западные
        self.nNW[self.barrierNW] = nSE_old[self.barrier] # северо-западные ← юго-восточные
        self.nSE[self.barrierSE] = nNW_old[self.barrier] # юго-восточные ← северо-западные
        self.nSW[self.barrierSW] = nNE_old[self.barrier] # юго-западные ← северо-восточные

    def collide(self):
        # Вычисляем макроскопические параметры
        current_rho = self.rho
        current_ux = self.ux
        current_uy = self.uy

        # Столкновение
        ux2 = current_ux * current_ux
        uy2 = current_uy * current_uy
        u2 = ux2 + uy2
        omu215 = 1 - 1.5 * u2
        uxuy = current_ux * current_uy

        vis = self.ini.viscosity # вязкость
        omega = 1 / (3 * vis + 0.5)  # "relaxation" parameter

        # Релаксация к равновесию
        self.n0 = (1 - omega) * self.n0 + omega * self.ini.four9ths * current_rho * omu215
        self.nN = (1 - omega) * self.nN + omega * self.ini.one9th * current_rho * (omu215 + 3 * current_uy + 4.5 * uy2)
        self.nS = (1 - omega) * self.nS + omega * self.ini.one9th * current_rho * (omu215 - 3 * current_uy + 4.5 * uy2)
        self.nE = (1 - omega) * self.nE + omega * self.ini.one9th * current_rho * (omu215 + 3 * current_ux + 4.5 * ux2)
        self.nW = (1 - omega) * self.nW + omega * self.ini.one9th * current_rho * (omu215 - 3 * current_ux + 4.5 * ux2)
        self.nNE = (1 - omega) * self.nNE + omega * self.ini.one36th * current_rho * (omu215 + 3 * (current_ux + current_uy) + 4.5 * (u2 + 2 * uxuy))
        self.nNW = (1 - omega) * self.nNW + omega * self.ini.one36th * current_rho * (omu215 + 3 * (-current_ux + current_uy) + 4.5 * (u2 - 2 * uxuy))
        self.nSE = (1 - omega) * self.nSE + omega * self.ini.one36th * current_rho * (omu215 + 3 * (current_ux - current_uy) + 4.5 * (u2 - 2 * uxuy))
        self.nSW = (1 - omega) * self.nSW + omega * self.ini.one36th * current_rho * (omu215 + 3 * (-current_ux - current_uy) + 4.5 * (u2 + 2 * uxuy))
        # -----------------------
        rho_in = 1.0
        u_in = self.ini.u0
        u2_in = u_in * u_in
        omu215_in = 1 - 1.5 * u2_in

        # Левая граница (x=0) - ВХОД
        self.nE[0, :] = self.ini.one9th * rho_in * (omu215_in + 3*u_in + 4.5*u2_in)
        self.nW[0, :] = self.ini.one9th * rho_in * (omu215_in - 3*u_in + 4.5*u2_in)
        self.nNE[0, :] = self.ini.one36th * rho_in * (omu215_in + 3*u_in + 4.5*u2_in)
        self.nSE[0, :] = self.ini.one36th * rho_in * (omu215_in + 3*u_in + 4.5*u2_in)
        self.nNW[0, :] = self.ini.one36th * rho_in * (omu215_in - 3*u_in + 4.5*u2_in)
        self.nSW[0, :] = self.ini.one36th * rho_in * (omu215_in - 3*u_in + 4.5*u2_in)
        self.nN[0, :] = self.ini.one9th * rho_in * omu215_in
        self.nS[0, :] = self.ini.one9th * rho_in * omu215_in
        self.n0[0, :] = self.ini.four9ths * rho_in * omu215_in
        # ------------------------------------------------------

    # Считаем вихри
    def curl(self):
        dv_dx = np.roll(self.uy, -1, axis=0) - np.roll(self.uy, 1, axis=0)
        # ∂u/∂y = производная ux по Y (axis=1)
        du_dy = np.roll(self.ux, -1, axis=1) - np.roll(self.ux, 1, axis=1)
        # Вихрь: ω = ∂v/∂x - ∂u/∂y
        vorticity = dv_dx - du_dy
        return vorticity

    def calculate_force(self):
        barrierFx = 0.0
        barrierFy = 0.0
        barrier_count = 0
        barrier_x_sum = 0
        barrier_y_sum = 0

        for x in range(1, self.ini.width - 1):
            for y in range(1, self.ini.height - 1):
                if self.barrier[x, y]:
                    # Сила по X (drag)
                    barrierFx += (
                            self.nE[x, y] +
                            self.nNE[x, y] +
                            self.nSE[x, y] -
                            self.nW[x, y] -
                            self.nNW[x, y] -
                            self.nSW[x, y]
                    )

                    # Сила по Y (lift)
                    barrierFy += (
                            self.nN[x, y] +
                            self.nNE[x, y] +
                            self.nNW[x, y] -
                            self.nS[x, y] -
                            self.nSE[x, y] -
                            self.nSW[x, y]
                    )

                    # Для центра силы
                    barrier_count += 1
                    barrier_x_sum += x
                    barrier_y_sum += y

        # Средняя позиция силы
        if barrier_count > 0:
            force_center_x = barrier_x_sum / barrier_count
            force_center_y = barrier_y_sum / barrier_count
        else:
            force_center_x = force_center_y = 0

        return barrierFx, barrierFy, force_center_x, force_center_y