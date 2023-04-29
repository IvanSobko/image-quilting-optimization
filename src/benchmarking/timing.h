#ifndef TEAM19_TIMING_H
#define TEAM19_TIMING_H
#include "ImageQuilting.h"

namespace timing {

void run_timing();
double rdtsc(ImageQuilting* quilting);

int read_files(char ***files, const char* directory, const char* filename_filter = NULL);

}  // namespace timing

#endif  //TEAM19_TIMING_H
