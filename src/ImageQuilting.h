#ifndef TEAM19_IMAGEQUILTING_H
#define TEAM19_IMAGEQUILTING_H

#include "ImgData.h"

class ImageQuilting {
   public:
    ImageQuilting() = delete;
    ImageQuilting(const ImgData& data) { mData = data; }

    ImgData Synthesis();
    void Transfer() {}

   private:
    ImgData mData;

    // Write a block from the source data to the output data specified their upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);

    // Synthesize a new texture sample by randomly choosing blocks
    ImgData RandomBlockPlacement();

    // Compute the left edge overlap between two left and right blocks specified by their upper-left corners
    double ComputeLeftEdgeOverlap(int leftY, int leftX, int rightY, int rightX);

    // Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
    ImgData OverlapConstraints();

};

#endif  //TEAM19_IMAGEQUILTING_H
