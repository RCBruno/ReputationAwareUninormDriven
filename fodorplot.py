import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def Ufodor(x, y, neutro):
    if 0.0 <= x <= neutro and 0.0 <= y <= neutro:
        return float(2*x*y)
    elif neutro <= x <= 1.0 and neutro <= y <= 1.0:
        return float(2.0*(x+y-x*y-neutro))
    else:
        return float((x+y)/2.0)

neutro = 0.5 # You can adjust the neutro value here
x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)
X, Y = np.meshgrid(x, y)
Z = np.array([[Ufodor(x, y, neutro) for x, y in zip(x_row, y_row)] for x_row, y_row in zip(X, Y)])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Ufodor')

#plt.title('Fodor Uninorm with e={}'.format(neutro))
plt.savefig("fodor.png", dpi=300, bbox_inches='tight', pad_inches=.1)
plt.show()
