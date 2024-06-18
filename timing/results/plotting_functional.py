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

    for filename in glob.glob(filenames):
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

    # Helpers
    best_performances, best_runtimes = get_best_performances_and_runtimes(glob.glob("*.txt"))

    # Ivan's cpu and compiler
    cpu_ivan = ""
    cpu = "AMD Ryzen 7 5800H @3.2GHz"
    compiler = "GCC 12.2.0"

    # Bounds refactor, loop reorder, and blocking
    filenames_and_labels_1 = [
        ("StdC_KUnroll_BoundsRefactor_-O1_06-06-20-29-14.txt", "Base Low"),
        ("StdC_KUnroll_BoundsRefactor_-O3-ffast-math-march=native_06-06-19-47-11.txt", "Base High"),
        ("StdC_KUnroll_BoundsRefactor_-O3-fno-tree-vectorize_06-06-20-05-17.txt", "Base Mid"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_-O1_06-06-20-30-28.txt", "Loop Reorder Low"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_-O3-ffast-math-march=native_06-06-19-47-59.txt", "Loop Reorder High"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_-O3-fno-tree-vectorize_06-06-20-06-31.txt", "Loop Reorder Mid"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128_-O1_06-06-20-36-53.txt", "128 Low"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128_-O3-ffast-math-march=native_06-06-19-51-51.txt", "128 High"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking128_-O3-fno-tree-vectorize_06-06-20-12-26.txt", "128 Mid"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32_-O1_06-06-20-31-42.txt", "32 Low"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32_-O3-ffast-math-march=native_06-06-19-48-40.txt", "32 High"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking32_-O3-fno-tree-vectorize_06-06-20-07-44.txt", "32 Mid"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48_-O1_06-06-20-32-56.txt", "48 Low"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48_-O3-ffast-math-march=native_06-06-19-49-29.txt", "48 High"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking48_-O3-fno-tree-vectorize_06-06-20-08-54.txt", "48 Mid"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64_-O3-ffast-math-march=native_06-06-19-50-20.txt", "64 High"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking64_-O3-fno-tree-vectorize_06-06-20-10-03.txt", "64 Mid"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96_-O1_06-06-20-35-34.txt", "96 Low"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96_-O3-ffast-math-march=native_06-06-19-51-06.txt", "96 High"),
        ("StdC_KUnroll_BoundsRefactor_LoopReorder_Blocking96_-O3-fno-tree-vectorize_06-06-20-11-14.txt", "96 Mid"),
    ]

    performance_plot(
       [(filename, label.split(" ")[0]) for (filename, label) in filenames_and_labels_1 if "Mid" in label],
        cpu, compiler, "../plots/bounds_reorder_blocking_performance.pdf")

    filenames_and_labels_2 = [
        ("default_-O3-ffast-math-march=native_06-06-19-43-23.txt", "Default"),
        ("StdC_Algorithm_-O3-ffast-math-march=native_06-06-19-45-36.txt", "+StdC+Alg"),
        ("StdC_Algorithm_Vectorize_-O3-ffast-math-march=native_06-06-19-46-48.txt", "++Vec"),
        ("StdC_KSrc2Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32_-O3-ffast-math-march=native_06-06-19-54-46.txt", "++Vec+32+L2"),
        ("StdC_KSrc4Unroll_Vector_BoundsRefactor_LoopReorder_Blocking32_-O3-ffast-math-march=native_06-06-19-55-13.txt", "++Vec+32+L4"),
    ]

    performance_plot(
        filenames_and_labels_2, cpu, compiler, "../plots/best_performances.pdf")
    runtime_plot(
        filenames_and_labels_2,cpu, compiler, "../plots/best_runtimes.pdf", True)

    # Best performances and runtimes for high compiler optimizations
    default_performance_runtime = (1.63626, 1.42794e+11, "Default")
    best_performances_runtimes = [
        (9.26256, 2.34679e+10, "+StdC+Alg"),
        (8.63746, 2.51663e+10, "++Vec"),
        (6.3838,  7.25578e+10, "++Vec+32+L2"),
        (7.75122, 2.30573e+10, "++Vec+32+L4"),
    ]
    best_performance_speedups = [
        (performance / default_performance_runtime[0], label)
        for (performance, runtime, label) in
        sorted(best_performances_runtimes, key=lambda x: -x[0])
    ]
    # print(best_performance_speedups)
    # [(5.6608118514172565, '+StdC+Alg'), (5.278782100644152, '++Vec'), (4.737156686590151, '++Vec+32+L4'),
    best_runtime_speedups = [
        (default_performance_runtime[1] / runtime, label)
        for (performance, runtime, label) in
        sorted(best_performances_runtimes, key=lambda x: x[1])
    ]
    # print(best_runtime_speedups)
    # [(6.19300611953698, '++Vec+32+L4'), (6.084651800970688, '+StdC+Alg'), (git 40164426236674, '++Vec'), (1.9680034400160975, '++Vec+32+L2')]
