import numpy as np
import matplotlib.pyplot as plt


def method_one():
    rows_num = 100
    cols_num = 100
    iteration = 1000
    k = 2  # k*L(T) = delta(T)/delta(t)
    h = 0.1  # Euler method step length

    T = np.zeros((rows_num, cols_num))
    colourMap = plt.cm.jet

    # boundary condition
    plt.figure()
    T[0, :] = 20  # top=20
    T[rows_num-1, :] = 40  # bottom=40
    T[:, 0] = 40  # left=40
    T[:, cols_num-1] = 40  # right=40

    for i in range(iteration):
        LT = T.copy()
        for x in range(1, rows_num-1):
            for y in range(1, cols_num-1):
                LT[x, y] = T[x+1][y] + T[x-1][y] + \
                    T[x][y+1] + T[x][y-1] - 4*T[x][y]
        T[1:rows_num-1, 1: cols_num-1] += h*k*LT[1:rows_num-1, 1: cols_num-1]
        plt.imshow(T, cmap=colourMap, vmin=0, vmax=40)
        plt.pause(0.0001)

    plt.imshow(T, cmap=colourMap, vmin=0, vmax=40)
    plt.colorbar()
    plt.show()


# method_one()


def method_two():
    rows_num = 100
    cols_num = 100
    iteration = 1000
    # k = 2  # k*L(T) = delta(T)/delta(t)
    # h = 0.1  # Euler method step length

    T = np.zeros((rows_num, cols_num))
    colourMap = plt.cm.jet

    # boundary condition
    plt.figure()
    T[0, :] = 20  # top=20
    T[rows_num-1, :] = 40  # bottom=40
    T[:, 0] = 40  # left=40
    T[:, cols_num-1] = 40  # right=40

    for i in range(iteration):
        for x in range(1, rows_num-1):
            for y in range(1, cols_num-1):
                T[x, y] = (T[x+1][y] + T[x-1][y] + T[x][y+1] + T[x][y-1])*0.25
        plt.imshow(T, cmap=colourMap, vmin=0, vmax=40)
        plt.pause(0.0001)

    plt.imshow(T, cmap=colourMap, vmin=0, vmax=40)
    plt.colorbar()
    plt.show()


method_two()

# reference:
# https://www.youtube.com/watch?v=QWtNLqF2Omc
# https://zhuanlan.zhihu.com/p/30312878
