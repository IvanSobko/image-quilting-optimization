import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# Note: set working directory to root of repo

# DO NOT modify here
cpu_ivan = "AMD Ryzen 7 5800H @3.2GHz"
cpu_baptiste = "Intel(R) Core(TM) i7-10510U @1.8GHz"
compiler_ivan = "Compiler: GCC 12.2.0"
compiler_baptiste = "Compiler: GCC 10.2.1"

# Modify here
cpu = cpu_baptiste
compiler = compiler_baptiste

def performance_plot(blockdiv):
    ax = plt.axes()
    ax.set_facecolor('#f0f0f0')
    plt.xlabel(xlabel="Pixel count")
    plt.ylabel(ylabel="Performance [F/C]", loc="top", rotation=0, labelpad=-126)
    plt.grid(axis = 'y', color = 'white')
    plt.title(cpu + ', ' + compiler, loc='left', pad=15)

    for filename in ["timing_results_-BASE_O3-ffast-math-march=native.txt", "timing_results_-O3-ffast-math-march=native.txt", "timing_results_-O3-fno-tree-vectorize.txt", "timing_results_-O1.txt"]:
        data = []
        file = open("./results/blockdiv_" + str(blockdiv) + "/" + filename, "r")
        for line in file:
            x_value = int(line.split('n=')[1].split(',')[0])
            y_value = float(line.split('performance=')[1].split(',')[0])
            data.append((x_value, y_value))

        data.sort(key=lambda x:x[0])

        plt.xscale('log', base=2)
        plt.xticks([2 ** np.floor((np.log2(i[0]))) for i in data])
        label = filename.split('_-')[1].split('.txt')[0]
        plt.plot(*zip(*data), '-o', label=label)
    plt.legend()
    plt.savefig("./results/blockdiv_" + str(blockdiv) + "/" + "performance_plot.png")

def runtime_plot(blockdiv):
    ax = plt.axes()
    ax.set_facecolor('#f0f0f0')
    plt.xlabel(xlabel="Pixel count")
    plt.ylabel(ylabel="Runtime [Cycles]", loc="top", rotation=0, labelpad=-126)
    plt.grid(axis = 'y', color = 'white')
    plt.title(cpu + ', ' + compiler, loc='left', pad=15)

    for filename in ["timing_results_-BASE_O3-ffast-math-march=native.txt", "timing_results_-O3-ffast-math-march=native.txt", "timing_results_-O3-fno-tree-vectorize.txt", "timing_results_-O1.txt"]:
        data = []
        file = open("./results/blockdiv_" + str(blockdiv) + "/" + filename, "r")
        for line in file:
            x_value = int(line.split('n=')[1].split(',')[0])
            y_value = float(line.split('cycles=')[1])
            data.append((x_value, y_value))

        data.sort(key=lambda x:x[0])

        plt.xscale('log', base=2)
        plt.xticks([2 ** np.floor((np.log2(i[0]))) for i in data])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        label = filename.split('_-')[1].split('.txt')[0]
        plt.plot(*zip(*data), '-o', label=label)
    plt.legend()
    plt.savefig("./results/blockdiv_" + str(blockdiv) + "/" + "runtime_plot.png")


if __name__ == "__main__":
    performance_plot(2)
    plt.clf()
    runtime_plot(2)
    plt.clf()

    performance_plot(4)
    plt.clf()
    runtime_plot(4)
    plt.clf()

    performance_plot(8)
    plt.clf()
    runtime_plot(8)
    plt.clf()