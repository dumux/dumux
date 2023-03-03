import matplotlib.pyplot as plt
from matplotlib import cm
from numpy import loadtxt, meshgrid, ndarray

x = loadtxt("x.txt")
y = loadtxt("y.txt")
z = loadtxt("z.txt")
zexact = loadtxt("zexact.txt")
X, Y = meshgrid(x, y)
Z = ndarray(shape=(len(y), len(x)))
ZExact = ndarray(shape=(len(y), len(x)))
ZDiff = ndarray(shape=(len(y), len(x)))
for j in range(len(y)):
    for i in range(len(x)):
        Z[j][i] = z[j*len(x) + i]
        ZExact[j][i] = zexact[j*len(x) + i]
        ZDiff[j][i] = ZExact[j][i] - Z[j][i]
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
surf = ax.plot_surface(X, Y, ZExact, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# surf = ax.plot_surface(X, Y, ZDiff, cmap=cm.coolwarm, linewidth=0, antialiased=False)
plt.show()
