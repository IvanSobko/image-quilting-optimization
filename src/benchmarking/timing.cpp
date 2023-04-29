#include "timing.h"

#include <dirent.h>
#include <cstdio>
#include "tsc_x86.h"

#define CYCLES_REQUIRED 1e8
#define RDTSC_LATENCY 26
#define REP 50

void timing::run_timing() {
    // TODO: run rdtsc() on images of different sizes. Then write results to file and create performance plot
    //  with python

    char** files_list = NULL;
    int n = read_files(&files_list, "./gallery", "input");

    for (int i = 0; i < n; i++) {
        printf("%s\n", files_list[i]);
    }


    for (int i = 0; i < n; i++) {
        free(files_list[i]);
    }
    free(files_list);
}

double timing::rdtsc(ImageQuilting* quilting) {
    double cycles = 0;
    int64_t num_runs = 100;
    double multiplier = 1;
    myInt64 start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    do {
        num_runs = num_runs * multiplier;
        start = start_tsc();
        for (size_t i = 0; i < num_runs; i++) {
            quilting->synthesis();
        }
        end = stop_tsc(start);

        cycles = (double)end;
        multiplier = (CYCLES_REQUIRED) / (cycles);

    } while (multiplier > 2);

    // Actual performance measurements repeated REP times.
    // We simply store all results and compute medians during post-processing.
    double total_cycles = 0;
    for (size_t j = 0; j < REP; j++) {
        start = start_tsc();
        for (size_t i = 0; i < num_runs; ++i) {
            quilting->synthesis();
        }
        end = stop_tsc(start) - RDTSC_LATENCY;

        cycles = (double)end / (double)num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    return total_cycles;
}
int timing::read_files(char*** files, const char* directory, const char* filename_filter) {
    DIR* folder = opendir(directory);
    if (folder == NULL) {
        perror("Unable to read directory");
        return NULL;
    }

    if (filename_filter == NULL) {
        filename_filter = "";  // is this even legal? p.s. modifying 'const char*'
    }

    struct dirent* entry;
    int n = 0;
    while ((entry = readdir(folder))) {
        if (strstr(entry->d_name, filename_filter) != NULL) {
            *files = (char**)realloc(*files, sizeof(**files) * (n + 1));
            (*files)[n] = strdup(entry->d_name);
            n++;
        }
    }
    closedir(folder);
    return n;
}
