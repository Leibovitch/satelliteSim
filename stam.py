from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt 
import numpy as np


ax = plt.axes(projection='3d')
A = np.linspace(0,15, 1000)
ax.plot3D(A, A, A)
# Axes3D.plot(, np.arange(0,10), np.arange(0,10))
# thet = np.arange(0,2*np.pi,1)
# s = np.sin(thet).reshape(len(thet),1)
# c = np.cos(thet).reshape(len(thet),1)
# z = np.zeros([len(thet), 1])

# A = np.array([s, c ,z])
# first_row = np.vstack((s,c,z))
# second_row = np.vstack((s,c,z))
# third_row = np.vstack((s,c,z))
# MAT = np.array([first_row, second_row, third_row])

print('shape')




