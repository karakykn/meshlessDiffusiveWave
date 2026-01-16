import numpy as np
import matplotlib.pyplot as plt

filepath = f'../segment0/geo/'

for i in range(31):
    xss = np.loadtxt(f'{filepath}xs{i}')
    print(xss[0,1])
    plt.plot([xss[0,1], xss[0,2]], [0, 0], 'k')
    for j in range(1, len(xss)):
        plt.plot([xss[j,1], xss[j-1,1]], [xss[j,0], xss[j-1,0]], 'k')
        plt.plot([xss[j, 2], xss[j-1, 2]], [xss[j, 0], xss[j-1, 0]], 'k')
    plt.show()