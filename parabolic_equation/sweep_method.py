import numpy as np
import matplotlib.pyplot as plt
import time

# Коэффициент a ** 2
def a(x, t):
    return 1/100000

# Источник
def f(x, t):
    return np.sin(t)

# Начальные условия
def init_cond(x):
    return 0

# Граничное условия слева
def bound_cond_l(t):
    return 1

# Граничное условие справа
def bound_cond_r(t):
    return 2


# Константы
l = 0.1  # длина
n = 100  # количество разбиений
t_max = 200  # конечное время
h = l / n  # шаг по координате
tau = 0.01  # шаг по времени
kind = 1  # род граничных условий

# Массив значений
u = np.zeros(n, dtype=np.double)

# Разбиение по x
x = np.arange(0, l, h)

# Начальные условия
u[:] = init_cond(x)

fig = plt.figure()
ax = fig.add_subplot(111)
li, = ax.plot(x, u)
plt.grid()
plt.show(block=False)
plt.xlim(0, l)
plt.ylim(0, 3)

def compute():
    global u, x
    t = 0
    alpha = np.zeros(n, dtype=np.double)
    beta = np.zeros(n, dtype=np.double)

    while t < t_max:
        # Коэффициенты для граничного условия первого рода
        if kind == 1:
            alpha[0] = 0
            u[0] = bound_cond_l(t)
            beta[0] = u[0]
        # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        elif kind == 2:
            alpha[0] = 1
            beta[0] = h * bound_cond_l(t) / a(0, t)

        for j in range(1, n - 1, 1):
            ai = ci = a(x[j], t) / (h ** 2)
            bi = 2 * a(x[j], t) / (h ** 2) + 1 / tau
            fi = -u[j] / tau - f(x[j], t)
            alpha[j] = ai / (bi - ci * alpha[j - 1])
            beta[j] = (ci * beta[j - 1] - fi) / (bi - ci * alpha[j - 1])

        # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        if kind == 1:
            downLimit = 0
            u[n - 1] = bound_cond_r(t)
        elif kind == 2:
            u[n - 1] = (a(l, t) * beta[n - 2] - h * bound_cond_r(t)) / (a(l, t) * (1 - alpha[n - 2]))
            downLimit = -1

        for j in range(n - 2, downLimit, -1):
            u[j] = alpha[j] * u[j + 1] + beta[j]
        ax.set_title('Время: {0:.3f} с.'.format(t))
        li.set_ydata(u)
        fig.canvas.draw()
        # plt.pause(0.001)

        t += tau


start_time = time.time()
compute()
print('Время выполнения: {}'.format(time.time() - start_time))
plt.pause(0)
