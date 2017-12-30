import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d



# Plot data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X = np.zeros( (11,11) )
Y = np.zeros( (11,11) )
Z = np.zeros( (11,11) )

for i in range(11):
    for j in range(11):
        X[i,j] = i - 10
        Y[i,j] = j - 10
        Z[i,j] = math.sin(i + j)
    
print(str(X))
print(str(Y))
print(str(Z))
ax.plot_wireframe(X, Y, Z)

plt.show()
