import numpy as np
import matplotlib.pyplot as plt
import time

l = 0.1 #м

#Параметры стали
lmbd = 46 #Вт/(м*К)
rho = 7800 #кг/м^3
c = 460 #Дж/(кг*К)
a = lmbd/(rho * c)

#Количество разбиений
N = 100

#Интервал
tmax = 20 #с

#Задание шага
h = l/N
#tau = (h ** 2)/(4 * a)
tau = 0.5

#Массив температур
T = np.zeros(N, dtype=np.double)

#Начальные условия
T0 = 20 #К

#Род граничных условия
kind = 2

#Граничные условия первого рода
if kind == 1:
    TL = 300 #К
    TR = 100 #К
    T[:] = T0
    T[0] = TL
    T[N - 1] = TR
#Граничные условия второго рода
elif kind == 2:
    QL = 30 #Вт
    QR = 10 #Вт
    QL /= rho*c
    QR /= rho*c

#Разбиение по x
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

    #Коэффициенты для граничного условия первого рода
    if kind == 1:
        alpha[0] = 0
        beta[0] = TL
    # Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
    elif kind == 2:
        alpha[0] = 1
        beta[0] = h*QL/a

    ai = ci = a/(h ** 2)
    bi = 2*a/(h ** 2) + 1/tau

    while t < tmax:
        for j in range(1, N - 1, 1):
            fi = -T[j]/tau - a/(2*(h ** 2))*(T[j - 1] - 2*T[j] + T[j + 1])
            alpha[j] = ai/(bi - ci*alpha[j - 1])
            beta[j] = (ci*beta[j - 1] - fi)/(bi - ci*alpha[j - 1])

        #Коэффициенты для граничного условия второго рода (первый порядок апроксимации)
        if kind == 1:
            downLimit = 0
        elif kind == 2:
            T[N - 1] = (a*beta[N - 2] - h*QR)/(a*(1 - alpha[N - 2]))
            downLimit = -1

        for j in range(N - 2, downLimit, -1):
            T[j] = alpha[j]*T[j + 1] + beta[j]
        ax.set_title('Время: {0:.2f} с.'.format(t))
        li.set_ydata(T)
        fig.canvas.draw()
        #plt.pause(0.001)

        t += tau

start_time = time.time()
compute()
print('Время выполнения: {}'.format(time.time() - start_time))
plt.pause(2)