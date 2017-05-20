import numpy as np
import matplotlib.pyplot as plt
import time

l = 0.1  # м
lmbd = 46  # Вт/(м*С)
rho = 7800  # кг/м^3
c = 460  # Дж/(кг*С)
T0 = 20  # С
TL = 300  # С
TR = 100  # C
tmax = 10  # с

N = 100

h = l / N
a = lmbd / (rho * c)
tau = (h ** 2) / (4 * a)

T = np.zeros(N, dtype=np.double)

T[:] = T0
T[0] = TL
T[N - 1] = TR


def getNewT(j):
    return T[j] + a * tau * (T[j + 1] - 2 * T[j] + T[j - 1]) / (h ** 2)

x = np.arange(0, l, h)

fig = plt.figure()
ax = fig.add_subplot(111)
li, = ax.plot(x, T)
plt.grid()
plt.show(block=False)


def compute():
    global T
    t = 0
    newT = T
    while t < tmax:
        for j in range(1, N - 1, 1):
            newT[j] = getNewT(j)
        T = newT.copy()

        ax.set_title('Время: {0:.2f} с.'.format(t))
        li.set_ydata(T)
        fig.canvas.draw()
        # plt.pause(0.001)

        t += tau


start_time = time.time()
compute()
print(time.time() - start_time)
plt.pause(0)
