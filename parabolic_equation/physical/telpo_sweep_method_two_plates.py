import numpy as np
import matplotlib.pyplot as plt
import time

l = 0.1  # м

# Параметры стали
lmbd1 = 46  # Вт/(м*С) - коэффициент теплопроводности
rho1 = 7800  # кг/м^3 - плотность
c1 = 460  # Дж/(кг*С) - удельная теплоемкость
l1 = 0.05  # м - толщина первого слоя

# Параметры золота
lmbd2 = 320  # Вт/(м*К) - коэффициент теплопроводности
rho2 = 19320  # кг/м^3 - плотность
c2 = 129  # Дж/(кг*К) - удельная теплоемкость

# Параметры граничных и начальных условий
T0 = 20  # С
TL = 300  # С
TR = 100  # C

# Интервал
t_max = 300  # с

# Количество разбиений
N = 100

h = l / N
a = lmbd1 / (rho1 * c1)
# tau = (h ** 2)/(4 * a)
tau = 0.5

# Номер ячейки раздела двух сред
N12 = round(l1 / h)

T = np.zeros(N, dtype=np.double)

T[:] = T0
T[0] = TL
T[N - 1] = TR

x = np.arange(0, l, h)

fig = plt.figure()
ax = fig.add_subplot(111)
li, = ax.plot(x, T)
plt.grid()
plt.show(block=False)


def compute():
    global T
    t = 0
    alpha = np.zeros(N, dtype=np.double)
    beta = np.zeros(N, dtype=np.double)

    alpha[0] = 0
    beta[0] = TL

    ai1 = ci1 = lmbd1 / (h ** 2)
    bi1 = 2 * lmbd1 / (h ** 2) + rho1 * c1 / tau

    ai2 = ci2 = lmbd2 / (h ** 2)
    bi2 = 2 * lmbd2 / (h ** 2) + rho2 * c2 / tau

    while t < t_max:
        for j in range(1, N12, 1):
            fi1 = -rho1 * c1 * T[j] / tau
            alpha[j] = ai1 / (bi1 - ci1 * alpha[j - 1])
            beta[j] = (ci1 * beta[j - 1] - fi1) / (bi1 - ci1 * alpha[j - 1])

        alpha[N12] = lmbd2 / (lmbd2 + lmbd1 * (1 - alpha[N12 - 1]))
        beta[N12] = lmbd1 * beta[N12 - 1] / (lmbd2 + lmbd1 * (1 - alpha[N12 - 1]))

        for j in range(N12 + 1, N - 1, 1):
            fi2 = -rho2 * c2 * T[j] / tau
            alpha[j] = ai2 / (bi2 - ci2 * alpha[j - 1])
            beta[j] = (ci2 * beta[j - 1] - fi2) / (bi2 - ci2 * alpha[j - 1])

        for j in range(N - 2, 0, -1):
            T[j] = alpha[j] * T[j + 1] + beta[j]

        ax.set_title('Время: {0:.2f} с.'.format(t))
        li.set_ydata(T)
        fig.canvas.draw()
        # plt.pause(0.001)

        t += tau


start_time = time.time()
compute()
print('Время выполнения: {}'.format(time.time() - start_time))
plt.pause(2)
