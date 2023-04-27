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

    // Synthesize a new texture sample by randomly choosing blocks
    ImgData RandomBlockPlacement();
};

#endif  //TEAM19_IMAGEQUILTING_H
