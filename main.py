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

        self.label_vel.setText(f'скорость потока = {self.data.initials.u0}')
        self.label_vis.setText(f'вязкость = {self.data.initials.viscosity}')

        # barrier type
        self.select_prop.addItems(['Квадрат', 'Цилиндр', 'Профиль крыла']) # barrier type
        self.select_prop.currentIndexChanged.connect(self.set_barrier_type)
        self.select_graph.addItems(['Плотность', 'x скорость', 'у скорость', 'скорость', 'завихрение']) # plot type
        self.select_graph.currentIndexChanged.connect(self.set_plot_type)
        self.select_cmap.addItems(cmaps) # цветовые схемы
        self.select_cmap.currentIndexChanged.connect(self.set_cmap) # обработка выбора меню cmap
        
        # рисование градиента
        gradient_array = np.linspace(0, 1, 200).reshape(-1, 1) # 1 строка, 200 столбцов
        gradient_array = np.repeat(gradient_array, 80, axis=1)  # повторяем 80 строк
        
        self.data.set_barrier()
        br = self.data.draw_barrier()
        br.setCompositionMode(mode=QPainter.CompositionMode.CompositionMode_Multiply)
        
        self.img_item = pg.ImageItem(gradient_array)
        self.img_item.setColorMap(pg.colormap.get('plasma', source='matplotlib'))
        self.canvas.clear()
        self.canvas.addItem(self.img_item)
        self.canvas.addItem(br)
        # self.canvas.setCompositionMode(QPainter.CompositionMode.CompositionMode_Darken)
        self.canvas.showGrid(x=False, y=False)
        self.canvas.hideAxis('left')
        self.canvas.hideAxis('bottom')
        
        
        # self.data.set_barriers()
        # self.canvas.addItem(self.data.draw_barriers())
        
        # timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.step)
        self.timer.stop()


        self.bt_start.clicked.connect(self.set_start)
        self.bt_stop.clicked.connect(self.set_stop)
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
        barrier = self.select_prop.currentText()
        if barrier == 'Квадрат':
            pass
        elif barrier == 'Цилиндр':
            pass
        elif barrier == 'Профиль крыла':
            pass
    
    def set_plot_type(self):
        plot = self.select_graph.currentText()
        print(plot)
        if plot == 'Плотность':
            return self.data.rho.T
        elif plot == 'x скорость':
            return self.data.ux.T
        elif plot == 'у скорость':
            return self.data.uy.T
        elif plot == 'скорость':
            return self.data.ux.T + self.data.uy.T
        elif plot == 'завихрение':
            return self.data.curl()
        
        
        # db_spin_velocity
    def step(self):
        # Эта функция вызывается АВТОМАТИЧЕСКИ каждые 50 мс
        # print('Таймер сработал!')
        self.data.stream();
        self.data.collide();
        graphic = self.set_plot_type()
        self.img_item = pg.ImageItem(graphic)
        self.canvas.clear()
        self.canvas.addItem(self.img_item)
        self.set_cmap()
        # self.canvas.addItem(br)

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
        

    def set_stop(self):
        self.bt_step.setEnabled(True)
        self.bt_start.setEnabled(True)
        self.bt_reset.setEnabled(True)
        self.bt_stop.setEnabled(False)
        self.timer.stop()
        print('stop')
            





if __name__ == '__main__':  # Если мы запускаем файл напрямую, а не импортируем
    app = widgets.QApplication(sys.argv)  # Новый экземпляр QApplication
    app.setStyle(widgets.QStyleFactory.create('Fusion'))  # Более современная тема оформления
    app.setPalette(widgets.QApplication.style().standardPalette())  # Берём цвета из темы оформления
    window = ExampleApp()  # Создаём объект класса ExampleApp
    window.show()  # Показываем окно
    app.exec()  # и запускаем приложение
