import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def performance_plot(filenames_and_labels, cpu, compiler, output_filename):
    # Plot parameters
    ax = plt.axes()
    ax.set_facecolor('#f0f0f0')
    plt.xlabel(xlabel="Pixel count")
    plt.ylabel(ylabel="Performance [F/C]", loc="top", rotation=0, labelpad=-126)
    plt.grid(axis='y', color='white')
    plt.title(cpu + '; ' + compiler, loc='left', pad=15)

    # Log scale on the x axis
    plt.xscale('log', base=2)

    # Iterate over the filename, label pairs
    for (filename, label) in filenames_and_labels:
        data = []
        file = open(filename, "r")
        for i, line in enumerate(file):
            # Skip the header line
            if i == 0:
                pass
            # Process the following lines
            else:
                x_value = int(line.split('n=')[1].split(',')[0])
                y_value = float(line.split('performance=')[1].split(',')[0])
                data.append((x_value, y_value))

        # Sort the data by the x value and plot it
        data.sort(key=lambda x: x[0])
        plt.xticks([2 ** np.floor((np.log2(i[0]))) for i in data])
        plt.plot(*zip(*data), '-o', label=label)

    # Display the legend and save the plot
    plt.legend()
    plt.savefig(output_filename)
    plt.clf()


def runtime_plot(filenames_and_labels, cpu, compiler, output_filename):
    # Plot parameters
    ax = plt.axes()
    ax.set_facecolor('#f0f0f0')
    plt.xlabel(xlabel="Pixel count")
    plt.ylabel(ylabel="Runtime [Cycles]", loc="top", rotation=0, labelpad=-126)
    plt.grid(axis='y', color='white')
    plt.title(cpu + '; ' + compiler, loc='left', pad=15)

    # Log scale on the x and y axis
    plt.xscale('log', base=2)
    plt.yscale('log', base=2)

    # Iterate over the filename, label pairs
    for (filename, label) in filenames_and_labels:
        data = []
        file = open(filename, "r")
        for i, line in enumerate(file):
            # Skip the header line
            if i == 0:
                pass
            # Process the following lines
            else:
                x_value = int(line.split('n=')[1].split(',')[0])
                y_value = float(line.split('cycles=')[1])
                data.append((x_value, y_value))

        # Sort the data by the x value
        data.sort(key=lambda x: x[0])
        plt.xticks([2 ** np.floor((np.log2(i[0]))) for i in data])
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2g'))
        plt.plot(*zip(*data), '-o', label=label)

    # Display the legend and save the plot
    plt.legend()
    plt.savefig(output_filename)
    plt.clf()


if __name__ == "__main__":

    # Ivan's cpu and compiler
    cpu_ivan = ""
    cpu = "AMD Ryzen 7 5800H @3.2GHz"
    compiler = "GCC 12.2.0"

    # Filenames and labels
    filenames_and_labels = [
        ("default_-O1_06-06-19-26-46.txt", "Default Low"),
        ("default_-O3-fno-tree-vectorize_06-06-19-26-47.txt", "Default Mid"),
        ("default_-O3-ffast-math-march=native_06-06-19-26-48.txt", "Default High"),
    ]

    # Generate the desired plots
    performance_plot(filenames_and_labels, cpu, compiler, "my_performance.png")
    runtime_plot(filenames_and_labels, cpu, compiler, "my_runtime.png")
