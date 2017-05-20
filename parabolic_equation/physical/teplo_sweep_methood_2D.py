import numpy as np
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D

# region Параметры стали
lmbd = 46  # Вт/(м*К)
rho = 7800  # кг/м^3
c = 460  # Дж/(кг*К)
a = lmbd / (rho * c)
# endregion
# region Размер и количество разбиений
lx = 0.4  # м
ly = 0.4  # м
NX = 40
NY = 40
# endregion

# Интервал
t_max = 500  # с

# region Задание шага
hx = lx / NX
hy = ly / NY
# tau = (hx ** 2)/(4 * a)
tau = 5
# endregion

# Массив температур
T = np.zeros((NX, NY))

# region Начальные условия
T0 = 20  # К
T[:, :] = T0
# endregion
# region Род граничных условия
kindXL = 2
kindXR = 2
kindYU = 2
kindYD = 2
# endregion
# region Граничные условия по оси X
if kindXL == 1:
    TL = 100  # К
    T[0, :] = TL
elif kindXL == 2:
    # При QL > 0 происходит нагревание
    QL = 3000  # Вт
    QL /= rho * c

if kindXR == 1:
    TR = 100  # К
    T[NX - 1, :] = TR
elif kindXR == 2:
    # При QR > 0 происходит охлаждение
    QR = -6000  # Вт
    QR /= rho * c
# endregion
# region Граничные условия по оси Y
if kindYU == 1:
    TU = 100  # К
    T[:, NY - 1] = TU
elif kindYU == 2:
    # При QU > 0 происходит нагревание
    QU = -20000  # Вт
    QU /= rho * c

if kindYD == 1:
    TD = 100  # К
    T[:, 0] = TD
elif kindYD == 2:
    # При QD > 0 происходит охлаждение
    QD = 3000  # Вт
    QD /= rho * c
# endregion

# Тип вывода (contour, 3d)
type = 'contour'

# region Формирование сетки для вывода
x = np.arange(0, lx, hx)
y = np.arange(0, ly, hy)
x, y = np.meshgrid(y, x)
# endregion
# region Подготовка вывода на экран
try:
    fig = plt.figure()
    if type == '3d':
        ax = fig.add_subplot(111, projection='3d')
        li = ax.plot_surface(x, y, T)
    elif type == 'contour':
        ax = fig.add_subplot(111)
        li = ax.contour(x, y, T)
# Костыль
except ValueError:
    pass

plt.show(block=False)
# endregion


def compute():
    t = 0
    global T, li

    while t < t_max:
        # region Проход вдоль оси X
        alpha = np.zeros((NX, NY), dtype=np.double)
        beta = np.zeros((NX, NY), dtype=np.double)

        # Коэффициенты для граничного условия первого рода
        if kindXL == 1:
            alpha[0, :] = 0
            beta[0, :] = TL
            downLimit = 0
        # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        elif kindXL == 2:
            alpha[0, :] = 1
            beta[0, :] = hx * QL / a
            downLimit = -1

        # Расчет коэффициентов метода прогонки
        ai = ci = a / (hx ** 2)
        bi = 2 * a / (hx ** 2) + 1 / tau

        # Нахождение всех коэффициентов alpha, beta
        for i in range(0, NY, 1):
            for j in range(1, NX - 1, 1):
                fi = -T[j, i] / tau
                alpha[j, i] = ai / (bi - ci * alpha[j - 1, i])
                beta[j, i] = (ci * beta[j - 1, i] - fi) / (bi - ci * alpha[j - 1, i])

        # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        if kindXR == 2:
            for i in range(0, NY, 1):
                T[NX - 1, i] = (a * beta[NX - 2, i] - hx * QR) / (a * (1 - alpha[NX - 2, i]))

        # Обратный ход метода прогонки
        for i in range(0, NY, 1):
            for j in range(NX - 2, downLimit, -1):
                T[j, i] = alpha[j, i] * T[j + 1, i] + beta[j, i]
        # endregion
        # region Проход вдоль оси Y
        alpha = np.zeros((NX, NY), dtype=np.double)
        beta = np.zeros((NX, NY), dtype=np.double)

        # Коэффициенты для граничного условия первого рода
        if kindYD == 1:
            alpha[:, 0] = 0
            beta[:, 0] = TD
            downLimit = 0
        # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        elif kindYD == 2:
            alpha[:, 0] = 1
            beta[:, 0] = hy * QD / a
            downLimit = -1

        # Расчет коэффициентов метода прогонки
        ai = ci = a / (hy ** 2)
        bi = 2 * a / (hy ** 2) + 1 / tau

        # Нахождение всех коэффициентов alpha, beta
        for i in range(1, NY - 1, 1):
            for j in range(0, NX, 1):
                fi = -T[j, i] / tau
                alpha[j, i] = ai / (bi - ci * alpha[j, i - 1])
                beta[j, i] = (ci * beta[j, i - 1] - fi) / (bi - ci * alpha[j, i - 1])

        # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        if kindYU == 2:
            for j in range(0, NX, 1):
                T[j, NY - 1] = (a * beta[j, NY - 2] - hy * QU) / (a * (1 - alpha[j, NY - 2]))

        # Обратный ход метода прогонки
        for i in range(NY - 2, downLimit, -1):
            for j in range(0, NX, 1):
                T[j, i] = alpha[j, i] * T[j, i + 1] + beta[j, i]
        # endregion
        # region Вывод на экран
        ax.clear()
        ax.set_title('Время: {0:.2f} с.'.format(t))
        plt.xlabel('y')
        plt.ylabel('x')

        if type == '3d':
            li = ax.plot_surface(x, y, T)
        elif type == 'contour':
            li = ax.contour(x, y, T)
            ax.clabel(li, fmt='%.1f')
        fig.canvas.draw()
        plt.pause(0.001)
        # endregion
        t += tau


start_time = time.time()
compute()
print('Время выполнения: {}'.format(time.time() - start_time))
plt.pause(-1)