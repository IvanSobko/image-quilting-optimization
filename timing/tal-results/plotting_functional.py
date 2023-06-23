import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import glob


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


def runtime_plot(filenames_and_labels, cpu, compiler, output_filename, ylogscale):
    # Plot parameters
    ax = plt.axes()
    ax.set_facecolor('#f0f0f0')
    plt.xlabel(xlabel="Pixel count")
    plt.ylabel(ylabel="Runtime [Cycles]", loc="top", rotation=0, labelpad=-126)
    plt.grid(axis='y', color='white')
    plt.title(cpu + '; ' + compiler, loc='left', pad=15)

    # Log scale on the x and y axis
    plt.xscale('log', base=2)
    if ylogscale:
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


def get_best_performances_and_runtimes(filenames):
    # Find the most performant and fastest runtime variants
    data_performances = {}
    best_performances = set()
    data_runtimes = {}
    best_runtimes = set()

    for filename in filenames:
        file = open(filename, "r")
        for i, line in enumerate(file):
            # Skip the header line
            if i == 0:
                pass
            # Process the following lines
            else:
                n = int(line.split('n=')[1].split(',')[0])
                cycles = float(line.split('cycles=')[1])
                performance = float(line.split('performance=')[1].split(',')[0])

                # Save the performance
                if n not in data_performances:
                    data_performances[n] = [(performance, filename)]
                else:
                    data_performances[n].append((performance, filename))

                # Save the runtime
                if n not in data_runtimes:
                    data_runtimes[n] = [(cycles, filename)]
                else:
                    data_runtimes[n].append((cycles, filename))

    # Sort the performances
    for n, performances_filenames in data_performances.items():
        performances_filenames.sort(key=lambda x: x[0])
        best_performances.add(performances_filenames[len(performances_filenames)-1][1])

    # Sort the runtimes
    for n, runtimes_filenames in data_runtimes.items():
        runtimes_filenames.sort(key=lambda x: x[0])
        best_runtimes.add(runtimes_filenames[0][1])

    return best_performances, best_runtimes

def get_best_performance_and_runtime_speedups(filenames):
    # Find the most performant and fastest runtime variants
    data_performances = {}
    best_performances = set()
    data_runtimes = {}
    best_runtimes = set()

    for filename in f:
        file = open(filename, "r")
        for i, line in enumerate(file):
            # Skip the header line
            if i == 0:
                pass
            # Process the following lines
            else:
                n = int(line.split('n=')[1].split(',')[0])
                cycles = float(line.split('cycles=')[1])
                performance = float(line.split('performance=')[1].split(',')[0])

                # Save the performance
                if n not in data_performances:
                    data_performances[n] = [(performance, filename)]
                else:
                    data_performances[n].append((performance, filename))

                # Save the runtime
                if n not in data_runtimes:
                    data_runtimes[n] = [(cycles, filename)]
                else:
                    data_runtimes[n].append((cycles, filename))

    # Sort the performances
    for n, performances_filenames in data_performances.items():
        performances_filenames.sort(key=lambda x: x[0])
        best_performances.add(performances_filenames[len(performances_filenames)-1][1])

    # Sort the runtimes
    for n, runtimes_filenames in data_runtimes.items():
        runtimes_filenames.sort(key=lambda x: x[0])
        best_runtimes.add(runtimes_filenames[0][1])

    return best_performances, best_runtimes


if __name__ == "__main__":

    # Tal's cpu and compiler
    cpu = "Intel(R) Core i7-1068NG7 @2.3GHz"
    compiler = "GCC 12.2.0"

    
    filenames_and_labels = [
        ("StdC_Algorithm_-O3-ffast-math-march=native_06-15-22-24-16.txt", "1/2"),
        ("StdC_Algorithm_-O3-ffast-math-march=native_06-15-22-30-03.txt", "1/4"),
        ("StdC_Algorithm_-O3-ffast-math-march=native_06-15-22-42-43.txt", "1/8"),
    ]

    performance_plot(
        filenames_and_labels, cpu, compiler, "../plots/248Performances.pdf")
    runtime_plot(
        filenames_and_labels,cpu, compiler, "../plots/248Runtimes.pdf", True)

