import numpy as np
import pandas as pd
import math
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d
import copy

r, n, d, m = 10**(-6), 377, 20, 4
E0, U0 = 8.85 * 10**(-12), 4 * math.pi * 10**(-7)
dt = (E0 * U0/2)**(1/2)
t, N, f, L, C, M = dt*100, 100, 3*10**7, 100, 1/(E0 * U0)**(1/2), 100
dx = dy = L / M
sxmax = -(m + 1) * math.log(r) / (2 * n * d * dx)
sx0 = -math.log(r) * math.log(3)/(2 * n * dx * (3**d - 1))
sx = []
sxm = []
sym = []
sy = []
data = []
H1 = []
Exnew = []
Eynew = []
S = []
Sm = []
phy = []
phx = []
pey = []
pex = []
aex = []
aey = []
ahx = []
ahy = []
bex = []
bey = []
bhx = []
bhy = []
for i in range(M+1):
    S.append([])
    Sm.append([])
    sx.append(0.0)
    sxm.append(0.0)
    sym.append(0.0)
    sy.append(0.0)
    H1.append([])
    Exnew.append([])
    Eynew.append([])
    phy.append([])
    phx.append([])
    pex.append([])
    pey.append([])
    aex.append([])
    aey.append([])
    ahx.append([])
    ahy.append([])
    bex.append([])
    bey.append([])
    bhx.append([])
    bhy.append([])
    for j in range(M+1):
        S[i].append(0.0)
        Sm[i].append(0.0)
        H1[i].append(0.0)
        Exnew[i].append(0.0)
        Eynew[i].append(0.0)
        phy[i].append(0.0)
        phx[i].append(0.0)
        pex[i].append(0.0)
        pey[i].append(0.0)

for i in range(0, 20):  # 电导率进行初始化
    sx[i + 1] = ((20 - i - 0.5) / d)**m * sxmax
    sxm[i] = ((20 - i) / d)**m * sxmax
    sy[i + 1] = ((20 - i - 0.5) / d)**m * sxmax
    sym[i] = ((20 - i) / d)**m * sxmax
    sx[M - i] = ((20 - i - 0.5)/d)**m * sxmax
    sxm[M - i] = ((20 - i)/d)**m * sxmax
    sy[M - i] = ((20 - i - 0.5)/d)**m * sxmax
    sym[M - i] = ((20 - i)/d)**m * sxmax

for i in range(M+1):
    aex[i] = math.e**(-sx[i] * dt/E0) - 1
    bex[i] = aex[i] + 1
    aey[i] = math.e**(-sy[i] * dt/E0) - 1
    bey[i] = aey[i] + 1
    ahx[i] = math.e**(-sxm[i] * dt/E0) - 1
    bhx[i] = ahx[i] + 1
    ahy[i] = math.e**(-sym[i] * dt/E0) - 1
    bhy[i] = ahy[i] + 1

for i in range(1, N):
    H1[50][50] = math.sin(2 * math.pi * f * i * dt) # 刷新点源，电源的位置可以自己设置（>20&<80)
    for j in range(1, M):
        for k in range(1, M):
            pey[j][k] = bey[j] * pey[j][k] + aey[j] * (H1[j][k] - H1[j-1][k])/dy
            pex[j][k] = bex[k] * pex[j][k] + aex[k] * (H1[j][k-1] - H1[j][k])/dx
            Exnew[j][k] = Exnew[j][k] + dt/(dy * E0) * (H1[j][k] - H1[j-1][k]) + dt/E0*pey[j][k]
            Eynew[j][k] = Eynew[j][k] + dt/(dx * E0) * (H1[j][k-1] - H1[j][k]) + dt/E0*pex[j][k]

    for j in range(0, M):
        for k in range(0, M):
            phy[j][k] = bhy[j] * phy[j][k] + ahy[j] * (Exnew[j+1][k] - Exnew[j][k])/dy
            phx[j][k] = bhx[k] * phx[j][k] + ahx[k] * (Eynew[j][k+1] - Eynew[j][k])/dx
            H1[j][k] = H1[j][k] + dt/(U0 * dy) * (Exnew[j+1][k] - Exnew[j][k]) - \
                       dt/(U0 * dx) * (Eynew[j][k+1] - Eynew[j][k]) + dt/U0 * (phy[j][k] - phx[j][k])
    data.append(copy.deepcopy(H1))
# 此处可以设置数据的输出方式，如果只想输出某个时刻的图像，去掉循环，将i（小于99）改为具体的值即可
# for i in range(0, 98, 5):
plt.imshow(data[95])
plt.show()



