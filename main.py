import sys

import numpy as np
from PyQt6 import QtWidgets as widgets, uic
from PyQt6.QtCore import QTimer
from PyQt6.QtGui import QPainter
import pyqtgraph as pg

from utils import cmaps
from engine import Data

# по умолчанию фон графиков pyqtgraph чёрный, делаем его белым
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

Design, _ = uic.loadUiType('design.ui')  # Подгружаем список виджетов из design.ui


class ExampleApp(widgets.QMainWindow, Design):
    def __init__(self):
        super().__init__()
        self.setupUi(self)  # Этот метод из класса Design, он инициализирует виджеты
        
        # класс с данными и массивами
        self.data = Data()
        # цвета кнопок
        self.bt_stop.setStyleSheet('background-color: rgb(255, 174, 174);')
        self.bt_stop.setEnabled(False)
        self.bt_start.setStyleSheet('background-color: rgb(191, 255, 190);')
        self.bt_reset.setStyleSheet('background-color: rgb(255, 255, 161);')

        self.label_vel.setText(f'скорость потока = {self.data.ini.u0}')
        self.label_vis.setText(f'вязкость = {self.data.ini.viscosity}')

        # barrier type
        self.select_prop.addItems(['без барьера', 'простой барьер', 'квадрат', 'цилиндр', 'профиль крыла']) # barrier type
        self.select_prop.currentIndexChanged.connect(self.set_barrier_type)
        self.select_graph.addItems(['плотность', 'x скорость', 'у скорость', 'скорость', 'завихрение']) # plot type
        self.select_graph.currentIndexChanged.connect(self.set_plot_type)
        self.select_cmap.addItems(cmaps) # цветовые схемы
        self.select_cmap.currentIndexChanged.connect(self.set_cmap) # обработка выбора меню cmap

        # self.data.set_simple_barrier()
        self.img_item = pg.ImageItem(self.set_plot_type())
        self.canvas.addItem(self.img_item)

        # self.data.set_barriers()
        # self.canvas.addItem(self.data.draw_barriers())
        
        # timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.step)
        self.timer.stop()

        self.bt_start.clicked.connect(self.set_start)
        self.bt_stop.clicked.connect(self.set_stop)
        self.bt_reset.clicked.connect(self.set_reset)
        self.bt_step.clicked.connect(self.step)
        self.slider_velocity.valueChanged.connect(self.set_velocity)
        self.slider_viscosity.valueChanged.connect(self.set_viscosity)

        # рисование стрелочки с силой 
        # drawForceArrow(barrierxSum/barrierCount, barrierySum/barrierCount, barrierFx, barrierFy);

    def set_velocity(self):
        d = float(self.slider_velocity.value())
        self.label_vel.setText(f'скорость потока = {d / 1000}')
        self.data.update_U0(d)
        # FIXME change data

    def set_viscosity(self):
        d = float(self.slider_viscosity.value())
        self.label_vis.setText(f'вязкость = {d / 1000}')
        self.data.update_viscosity(d)
        # FIXME change data

    def set_cmap(self):
        cmap = self.select_cmap.currentText()
        self.img_item.setColorMap(pg.colormap.get(cmap, source='matplotlib'))
        
    def set_barrier_type(self):
        barr = self.select_prop.currentText()
        if barr == 'простой барьер':
            self.data.set_simple_barrier()
            self.draw_plot()
        elif barr == 'квадрат':
            self.data.set_square_barrier()
            self.draw_plot()
        elif barr == 'цилиндр':
            self.data.set_circle_barrier()
            self.draw_plot()
        elif barr == 'профиль крыла':
            pass
        elif barr == 'без барьера':
            self.data.reset_barrier()
            self.draw_plot()
    
    def set_plot_type(self):
        plot = self.select_graph.currentText()
        if plot == 'Плотность':
            return self.data.rho
        elif plot == 'x скорость':
            return self.data.ux
        elif plot == 'у скорость':
            return self.data.uy
        elif plot == 'скорость':
            return np.sqrt(self.data.ux**2 + self.data.uy**2)
        elif plot == 'завихрение':
            return self.data.curl()



    def step(self):
        # Эта функция вызывается АВТОМАТИЧЕСКИ каждые 50 мс
        print('Таймер сработал!')
        self.data.collide()
        self.data.stream()
        self.draw_plot()

        # Проверка средних скоростей
        mean_ux = np.mean(self.data.ux)
        mean_uy = np.mean(self.data.uy)
        print(f"Шаг: средняя ux={mean_ux:.6f}, uy={mean_uy:.6f}")

        # Проверка на левой и правой границах
        left_flow = np.mean(self.data.nE[0, :] - self.data.nW[0, :])
        right_flow = np.mean(self.data.nE[-1, :] - self.data.nW[-1, :])
        print(f"  Поток слева: {left_flow:.6f}, справа: {right_flow:.6f}")

        print(f"u0 = {self.data.ini.u0}")

    def set_start(self):
        self.bt_step.setEnabled(False)
        self.bt_start.setEnabled(False)
        self.bt_reset.setEnabled(False)
        self.bt_stop.setEnabled(True)
        # Шаг 3: Запускаем таймер (например, каждые 50 мс)
        # Теперь каждые 50 миллисекунд будет вызываться self.step()
        self.timer.start(50)
        print('start')
        self.step()

    def set_reset(self):
        self.data.n0  = self.data.ini.four9ths *(np.ones((self.data.ini.width, self.data.ini.height)) - 1.5 * self.data.ini.u0 ** 2)
        self.data.nN  = self.data.ini.one9th *  (np.ones((self.data.ini.width, self.data.ini.height)) - 1.5 * self.data.ini.u0 ** 2)
        self.data.nS  = self.data.ini.one9th *  (np.ones((self.data.ini.width, self.data.ini.height)) - 1.5 * self.data.ini.u0 ** 2)
        self.data.nE  = self.data.ini.one9th *  (np.ones((self.data.ini.width, self.data.ini.height)) + 3  *  self.data.ini.u0 + 4.5 * self.data.ini.u0 ** 2 - 1.5 * self.data.ini.u0 ** 2)
        self.data.nW  = self.data.ini.one9th *  (np.ones((self.data.ini.width, self.data.ini.height)) - 3  *  self.data.ini.u0 + 4.5 * self.data.ini.u0 ** 2 - 1.5 * self.data.ini.u0 ** 2)
        self.data.nNE = self.data.ini.one36th * (np.ones((self.data.ini.width, self.data.ini.height)) + 3  *  self.data.ini.u0 + 4.5 * self.data.ini.u0 ** 2 - 1.5 * self.data.ini.u0 ** 2)
        self.data.nSE = self.data.ini.one36th * (np.ones((self.data.ini.width, self.data.ini.height)) + 3  *  self.data.ini.u0 + 4.5 * self.data.ini.u0 ** 2 - 1.5 * self.data.ini.u0 ** 2)
        self.data.nNW = self.data.ini.one36th * (np.ones((self.data.ini.width, self.data.ini.height)) - 3  *  self.data.ini.u0 + 4.5 * self.data.ini.u0 ** 2 - 1.5 * self.data.ini.u0 ** 2)
        self.data.nSW = self.data.ini.one36th * (np.ones((self.data.ini.width, self.data.ini.height)) - 3  *  self.data.ini.u0 + 4.5 * self.data.ini.u0 ** 2 - 1.5 * self.data.ini.u0 ** 2)
        self.draw_plot()


    def set_stop(self):
        self.bt_step.setEnabled(True)
        self.bt_start.setEnabled(True)
        self.bt_reset.setEnabled(True)
        self.bt_stop.setEnabled(False)
        self.timer.stop()
        print('stop')

    def draw_plot(self):
        self.img_item = pg.ImageItem(self.set_plot_type())
        self.canvas.clear()
        self.canvas.addItem(self.img_item)
        br = self.data.draw_barrier()
        br.setCompositionMode(mode=QPainter.CompositionMode.CompositionMode_Multiply)
        self.canvas.addItem(br)
        self.set_cmap()

if __name__ == '__main__':  # Если мы запускаем файл напрямую, а не импортируем
    app = widgets.QApplication(sys.argv)  # Новый экземпляр QApplication
    app.setStyle(widgets.QStyleFactory.create('Fusion'))  # Более современная тема оформления
    app.setPalette(widgets.QApplication.style().standardPalette())  # Берём цвета из темы оформления
    window = ExampleApp()  # Создаём объект класса ExampleApp
    window.show()  # Показываем окно
    app.exec()  # и запускаем приложение
