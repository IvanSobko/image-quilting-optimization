#include "AdvisorHelper.h"
#include "ImageQuilting.h"
#include "AdvanceAlgOptimiz.h"
#include "CompOverlapOptimiz.h"
#include "Blocking.h"

void Advisor::baseline(ImgData* imgData, int seed) {
    ImageQuilting imageQuilting(imgData);
    imageQuilting.Synthesis(0);
}

void Advisor::basicOpt(ImgData* imgData, int seed) {
    CompOverlapOptimiz::BasicOpt(imgData, seed);
}

void Advisor::unrollChnls(ImgData* imgData, int seed){
    CompOverlapOptimiz::UnrollChnls(imgData, seed);
}

void Advisor::unrollMemory(ImgData* imgData, int seed){
    AdvanceAlgOptimiz::DividedFuncOpt(imgData, seed);
}

void Advisor::unroll(ImgData* imgData, int seed) {
    CompOverlapOptimiz::UnrollOpt(imgData, seed);
}   

void Advisor::vectorize(ImgData* imgData, int seed) {
    CompOverlapOptimiz::VectorizeOpt(imgData, seed);
}

void Advisor::block(ImgData* imgData, int seed) {
    Blocking::Base(imgData, seed);
}