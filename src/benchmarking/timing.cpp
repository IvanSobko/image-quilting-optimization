#include "timing.h"
#include "tsc_x86.h"

#define CYCLES_REQUIRED 1e8
#define RDTSC_LATENCY 26
#define REP 50

void timing::run_timing()
{
 // TODO: run rdtsc() on images of different sizes. Then write results to file and create performance plot
 //  with python
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
            quilting->Synthesis();
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
            quilting->Synthesis();
        }
        end = stop_tsc(start) - RDTSC_LATENCY;

        cycles = (double)end / (double)num_runs;
        total_cycles += cycles;
    }
    total_cycles /= REP;
    return  total_cycles;
}