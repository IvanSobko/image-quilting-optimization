import numpy as np
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


    data = []
    file = open("./gallery/results.txt", "r")
    for line in file:
        x_value = int(line.split('n=')[1].split(',')[0])
        y_value = float(line.split('performance=')[1])
        data.append((x_value, y_value))

    data.sort(key=lambda x:x[0])
    print(data)

    plt.xscale('log', base=2)
    plt.xticks([2 ** np.floor((np.log2(i[0]))) for i in data])
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%g'))
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%g'))
    # plt.minorticks_off()

    plt.plot(*zip(*data), '-o', label='mincut')
    plt.legend()
    plt.savefig("./gallery/performance_plot.png")
    plt.show()


if __name__ == "__main__":
    performance_plot()