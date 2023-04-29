import numpy as np
import re
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# Note: set working directory to root of repo

def performance_plot():
    ax = plt.axes()
    ax.set_facecolor('#f0f0f0')
    plt.xlabel(xlabel="Pixel count")
    plt.ylabel(ylabel="Performance [F/C]", loc="top", rotation=0, labelpad=-126)
    plt.grid(axis = 'y', color = 'white')

    # don't forget to modify these if you run on your pc
    title_cpu = "AMD Ryzen 7 5800H @3.2GHz"
    title_compiler = "Compiler: GCC 12.2.0"
    plt.title(title_cpu + ', ' + title_compiler, loc='left', pad=15)

    x = []
    y = []
    file = open("./gallery/results.txt", "r")
    for line in file:
        x.append(float(line.split('n=')[1].split(',')[0]))
        y.append(float(line.split('performance=')[1]))

    x.sort()
    y.sort()

    plt.xscale('log', base=2)
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    # plt.minorticks_off()

    plt.plot(x, y, '-o', label='random blocks')
    plt.legend()
    plt.savefig("./gallery/performance_plot.png")
    plt.show()


if __name__ == "__main__":
    performance_plot()