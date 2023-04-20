#ifndef TEAM19_IMAGEQUILTING_H
#define TEAM19_IMAGEQUILTING_H

#include "ImgData.h"

class ImageQuilting {
public:
    ImageQuilting() = delete;
    ImageQuilting(const ImgData &data) { mData = data; }

    ImgData synthesis();
    void transfer() {}

private:
    ImgData mData;
};


#endif //TEAM19_IMAGEQUILTING_H
