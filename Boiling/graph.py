import matplotlib.pyplot as plt
import numpy as np
import sys

def ReadData(filename : str):
    x = set()
    y = set()
    v = []

    with open (filename, 'r') as file:
        for line in file:
            l = line.replace(',', '.')
            xCoord, yCoord, value = l.split(' ')
            x.add(float(xCoord))
            y.add(float(yCoord))
            v.append(float(value))

    return np.asarray(sorted(x)), np.asarray(sorted(y)), np.asarray(v)

def main():
    filename = sys.argv[1]

    x, y, v = ReadData(filename)
    X, Y = np.meshgrid(x, y)
    Z = v.reshape(len(y), len(x))

    plt.figure(1)
    cb = plt.contourf(X, Y, Z, levels = 10, cmap='jet')
    plt.colorbar(cb)
    plt.show()

if __name__ == "__main__":
    main()