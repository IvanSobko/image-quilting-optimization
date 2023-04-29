#include "timing.h"

#include <dirent.h>

#include "PngReader.h"
#include "tsc_x86.h"

//TODO: should we modify cycles required? 1e8 is value from homeworks
#define CYCLES_REQUIRED 1e8
#define RDTSC_LATENCY 26
#define REP 1

void timing::run_timing() {
    std::string directory = "./gallery";
    std::vector<std::string> files = read_files(directory, "input0_");
    if (files.empty()) {
        perror("Failed to collect files for timing");
        return;
    }

    FILE* results_txt = NULL;
    std::string results_path = directory + '/' + "results.txt";
    results_txt = fopen(results_path.c_str(), "w");

    ImgData img_data;
    for (int i = 0; i < files.size(); i++) {
        printf("Timing for %s\n", files[i].c_str());
        file::read_png_file(files[i].c_str(), img_data);
        //TODO: specify block size implicitly
        img_data.output_w = img_data.width * 2;
        img_data.output_h = img_data.height * 2;
        img_data.block_w = img_data.width / 4;
        img_data.block_h = img_data.height / 4;

        ImageQuilting quilting(img_data);
        double cycles = rdtsc(&quilting);
        printf("Done. Quilting for %s took %f cycles\n\n", files[i].c_str(), cycles);

        int64_t data_size = img_data.width * img_data.height;
        //TODO: write performance=flops/cycles after implementation is finished and we can count flops
        fprintf(results_txt, "size=%ix%i, n=%lli, performance=%f\n", img_data.width, img_data.height, data_size,
                (double)data_size/cycles);

        for (int j = 0; j < img_data.height; j++) {
            free(img_data.data[j]);
        }
        free(img_data.data);
        img_data.data = NULL;

        for (int j = 0; j < img_data.output_h; j++) {
            free(img_data.output_d[j]);
        }
        free(img_data.output_d);
        img_data.output_d = NULL;
    }
    fclose(results_txt);
}

double timing::rdtsc(ImageQuilting* quilting) {
    double cycles = 0;
    int64_t num_runs = 100;
    double multiplier = 1;
    myInt64 start, end;

    // Warm-up phase: we determine a number of executions that allows
    // the code to be executed for at least CYCLES_REQUIRED cycles.
    // This helps excluding timing overhead when measuring small runtimes.
    printf("Doing warmup phase...");
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
    printf("actually measuring performance (%lli times).\n", num_runs);
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

std::vector<std::string> timing::read_files(const std::string& directory,
                                            const std::string& filename_filter) {
    DIR* folder = opendir(directory.c_str());
    if (folder == NULL) {
        perror("Unable to read directory");
        return {};
    }

    std::vector<std::string> result;
    struct dirent* entry;
    while ((entry = readdir(folder))) {
        if (strstr(entry->d_name, filename_filter.c_str()) != NULL) {
            std::string filepath = directory + '/' + entry->d_name;
            result.push_back(filepath);
        }
    }
    closedir(folder);
    return result;
}
