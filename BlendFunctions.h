#ifndef BlendFunctions_h_
#define BlendFunctions_h_
#include "include/GBlendMode.h"
#include "include/GPixel.h"



GPixel clear(GPixel src, GPixel dst);

GPixel blendSrc(GPixel src, GPixel dst);
GPixel blendDst(GPixel src, GPixel dst);

GPixel blendSrcOver(GPixel src, GPixel dst);
GPixel blendDstOver(GPixel src, GPixel dst);

GPixel blendSrcIn(GPixel src, GPixel dst);

GPixel blendDstTop(GPixel src, GPixel dst);
GPixel blendXor(GPixel src, GPixel dst);
GPixel blendDstIn(GPixel src, GPixel dst);

GPixel blendSrcOut(GPixel src, GPixel dst);
GPixel blendDstOut(GPixel src, GPixel dst);

GPixel blendSrcTop(GPixel src, GPixel dst);

using BlendFunc = GPixel (*)(GPixel src, GPixel dst);
BlendFunc getBlendFunction(GBlendMode blendMode);

#endif