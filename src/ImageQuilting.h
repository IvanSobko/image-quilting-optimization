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

    // Write a block from the source data to the output data specified by the given upper-left corners
    void WriteBlock(int dstY, int dstX, int srcY, int srcX);

    // Synthesize a new texture sample by randomly choosing blocks
    ImgData RandomBlockPlacement();

    // Synthesize a new texture sample by randomly choosing blocks satisfying overlap constraints
    ImgData OverlapConstraints();

};

#endif  //TEAM19_IMAGEQUILTING_H
