import numpy as np
import matplotlib.pyplot as plt
import time

# Коэффициент a ** 2
def a(x, t):
    return 1

# Источник
def f(x, t):
    return 0  # np.sin(t)

# Начальные условия
def init_cond_u(x):
    return 0  # np.sin(x/l*2*np.pi)

# Начальные условия
def init_cond_du(x):
    return 0

# Граничное условия слева
def bound_cond_l(t):
    return 0

# Граничное условие справа
def bound_cond_r(t):
    return 0


# Константы
l = 0.1  # длина
n = 1000  # количество разбиений
t_max = 200  # конечное время
h = l / n  # шаг по координате
tau = 0.01  # шаг по времени
kind = 1  # род граничных условий
sigma = 0.25  # коэффициент схемы

# Массив значений
u_old = np.zeros(n + 1, dtype=np.double)
u_now = np.zeros(n + 1, dtype=np.double)
u_new = np.zeros(n + 1, dtype=np.double)

# Разбиение по x
x = np.arange(0, l + h, h)

# Начальные условия
u_old[:] = init_cond_u(x)
u_now[:] = u_old + init_cond_du(x)*tau
u_now[500] = 1
u_old[500] = 1

# Граничные условия
u_old[0] = bound_cond_l(0)
u_old[n] = bound_cond_r(0)
u_now[0] = bound_cond_l(tau)
u_now[n] = bound_cond_r(tau)

fig = plt.figure()
ax = fig.add_subplot(111)
li, = ax.plot(x, u_now)
plt.grid()
plt.show(block=False)
plt.xlim(0, l)
plt.ylim(-3, 3)

def compute():
    global u_old, u_now, u_new, x, sigma
    t = 0
    alpha = np.zeros(n, dtype=np.double)
    beta = np.zeros(n, dtype=np.double)

    while t < t_max:
        # Коэффициенты для граничного условия первого рода
        if kind == 1:
            alpha[0] = 0
            u_new[0] = bound_cond_l(t)
            beta[0] = u_new[0]
        # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        elif kind == 2:
            alpha[0] = 1
            beta[0] = h * bound_cond_l(t) / a(0, t)

        for j in range(1, n - 1, 1):
            gamma = a(x[j], t) / (h ** 2)
            ai = ci = gamma*sigma
            bi = 2*gamma*sigma + 1/(tau ** 2)
            fi = 1/(tau ** 2)*(u_old[j] - 2*u_now[j]) + gamma*((2*sigma - 1)*(u_now[j + 1] - 2*u_now[j] + u_now[j - 1]) - \
                 sigma*(u_old[j + 1] - 2*u_old[j] + u_old[j - 1])) - f(x[j], t)
            alpha[j] = ai / (bi - ci * alpha[j - 1])
            beta[j] = (ci * beta[j - 1] - fi) / (bi - ci * alpha[j - 1])

        # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        if kind == 1:
            downLimit = 0
            u_new[n - 1] = bound_cond_r(t)
        elif kind == 2:
            u_new[n - 1] = (a(l, t) * beta[n - 2] - h * bound_cond_r(t)) / (a(l, t) * (1 - alpha[n - 2]))
            downLimit = -1

        for j in range(n - 2, downLimit, -1):
            u_new[j] = alpha[j] * u_new[j + 1] + beta[j]
        ax.set_title('Время: {0:.3f} с.'.format(t))
        li.set_ydata(u_now)
        fig.canvas.draw()
        # plt.pause(0.001)

        u_old = u_now.copy()
        u_now = u_new.copy()

        t += tau


start_time = time.time()
compute()
print('Время выполнения: {}'.format(time.time() - start_time))
plt.pause(0)
