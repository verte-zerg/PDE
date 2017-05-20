import numpy as np
import matplotlib.pyplot as plt
import time

# Коэффициент a ** 2
def a(x, t):
    return 1

# Источник
def f(x, t):
    return np.sin(t)

# Начальные условия для функции
def init_cond_1(x):
    return np.sin(x/l*2*np.pi)

# Начальные условия для производной
def init_cond_2(x):
    return 0

# Граничное условия слева
def bound_cond_l(t):
    return 0

# Граничное условие справа
def bound_cond_r(t):
    return 0


# Константы
l = 0.1  # длина
n = 100  # количество разбиений
t_max = 200  # конечное время
h = l / n  # шаг по координате
tau = h/2  # шаг по времени
kind = 1  # род граничных условий

# Массивы значения на временных слоях
u_old = np.zeros(n + 1, dtype=np.double)
u_now = np.zeros(n + 1, dtype=np.double)
u_new = np.zeros(n + 1, dtype=np.double)

# Разбиение по x
x = np.arange(0, l + h, h)

# Начальные условия
u_old[:] = init_cond_1(x)
u_now[:] = u_old + init_cond_2(x)*tau

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
plt.ylim(-2, 2)

def compute():
    global u_old, u_now, u_new, x
    t = 2*tau

    while t < t_max:
        u_new[0] = bound_cond_l(t)
        u_new[n] = bound_cond_r(t)

        for i in range(1, n):
            alpha = a(x[i], t)*(tau ** 2)/(h ** 2)
            u_new[i] = alpha*(u_now[i + 1] + u_now[i - 1]) + 2*(1 - alpha)*u_now[i] - u_old[i]

        u_old = u_now.copy()
        u_now = u_new.copy()

        ax.set_title('Время: {0:.3f} с.'.format(t))
        li.set_ydata(u_now)
        fig.canvas.draw()
        # plt.pause(0.001)

        t += tau


start_time = time.time()
compute()
print('Время выполнения: {}'.format(time.time() - start_time))
plt.pause(0)
