#ifndef TEAM19_TIMING_H
#define TEAM19_TIMING_H
#include "ImageQuilting.h"

#include <string>
#include <vector>

namespace timing {

void run_timing();
double rdtsc(ImageQuilting* quilting);

std::vector<std::string> read_files(const std::string &directory, const std::string &filename_filter);

}  // namespace timing

#endif  //TEAM19_TIMING_H
