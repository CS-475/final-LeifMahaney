#include "BlendFunctions.h"


using BlendFunc = GPixel (*)(GPixel src, GPixel dst);

GPixel clear(GPixel src, GPixel dst) {
    return GPixel_PackARGB(0, 0, 0, 0);
}


GPixel blendSrcOver(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    if (srcA == 0) {
        return dst;
    }
    if (srcA == 255) {
        return src;
    }
    uint8_t dstA = GPixel_GetA(dst);
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t srcAlphaComp = 255 - srcA;

    uint8_t resultA = srcA + ((dstA * srcAlphaComp) >> 8);
    uint8_t resultR = srcR + ((dstR * srcAlphaComp) >> 8);
    uint8_t resultG = srcG + ((dstG * srcAlphaComp) >> 8);
    uint8_t resultB = srcB + ((dstB * srcAlphaComp) >> 8);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}



GPixel blendDstOver(GPixel src, GPixel dst) {
    uint8_t dstA = GPixel_GetA(dst);
    if (dstA == 0){
        return src;
    }
    if (dstA == 255){
        return dst;
    }

    uint8_t srcA = GPixel_GetA(src);
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    uint8_t dstAlphaComp = 255 - dstA;

    uint8_t resultA = dstA + ((srcA * dstAlphaComp) >> 8);
    uint8_t resultR = dstR + ((srcR * dstAlphaComp) >> 8);
    uint8_t resultG = dstG + ((srcG * dstAlphaComp) >> 8);
    uint8_t resultB = dstB + ((srcB * dstAlphaComp) >> 8);
    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}



GPixel blendSrcIn(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t dstA = GPixel_GetA(dst);

    if (srcA == 0 || dstA == 0) {
        return GPixel_PackARGB(0, 0, 0, 0);
    }
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);
    uint8_t resultA = (srcA * dstA) >> 8;
    uint8_t resultR = (srcR * dstA) >> 8;
    uint8_t resultG = (srcG * dstA) >> 8;
    uint8_t resultB = (srcB * dstA) >> 8;

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}



GPixel blendDstIn(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t dstA = GPixel_GetA(dst);

    if (srcA == 0 || dstA == 0) {
        return GPixel_PackARGB(0, 0, 0, 0);
    }
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);
    
    uint8_t resultA = (dstA * srcA) >> 8;
    uint8_t resultR = (dstR * srcA) >> 8;
    uint8_t resultG = (dstG * srcA) >> 8;
    uint8_t resultB = (dstB * srcA) >> 8;

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}




GPixel blendSrcOut(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    if (srcA == 0){
        return GPixel_PackARGB(0, 0, 0, 0);
    }
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);
    uint8_t dstA = GPixel_GetA(dst);

    uint8_t dstAlphaComp = 255 - dstA;

    uint8_t resultA = (srcA * dstAlphaComp) >> 8;
    uint8_t resultR = (srcR * dstAlphaComp) >> 8;
    uint8_t resultG = (srcG * dstAlphaComp) >> 8;
    uint8_t resultB = (srcB * dstAlphaComp) >> 8;

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}



GPixel blendDstOut(GPixel src, GPixel dst) {
    uint8_t dstA = GPixel_GetA(dst);
    if (dstA == 0){
        return GPixel_PackARGB(0, 0, 0, 0);
    }
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);
    uint8_t srcA = GPixel_GetA(src);

    uint8_t srcAlphaComp = 255 - srcA;

    uint8_t resultA = (dstA * srcAlphaComp) >> 8;
    uint8_t resultR = (dstR * srcAlphaComp) >> 8;
    uint8_t resultG = (dstG * srcAlphaComp) >> 8;
    uint8_t resultB = (dstB * srcAlphaComp) >> 8;

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}

GPixel blendSrc(GPixel src, GPixel dst) {
    return src;
}

GPixel blendDst(GPixel src, GPixel dst) {
    return dst;
}


GPixel blendXor(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstA = GPixel_GetA(dst);
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    float srcAlphaComp = 1.0f - (srcA / 255.0f);
    float dstAlphaComp = 1.0f - (dstA / 255.0f);

    uint8_t resultA = static_cast<uint8_t>(srcA * dstAlphaComp + dstA * srcAlphaComp + 0.5f);
    uint8_t resultR = static_cast<uint8_t>(srcR * dstAlphaComp + dstR * srcAlphaComp + 0.5f);
    uint8_t resultG = static_cast<uint8_t>(srcG * dstAlphaComp + dstG * srcAlphaComp + 0.5f);
    uint8_t resultB = static_cast<uint8_t>(srcB * dstAlphaComp + dstB * srcAlphaComp + 0.5f);

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}



GPixel blendSrcTop(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstA = GPixel_GetA(dst);
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    float dstAlpha = dstA / 255.0f;
    float srcAlphaComp = 1.0f - (srcA / 255.0f);

    uint8_t resultA = static_cast<uint8_t>(std::min(255.0f, (srcA * dstAlpha) + (dstA * srcAlphaComp) + 0.5f));
    uint8_t resultR = static_cast<uint8_t>((srcR * dstAlpha) + (dstR * srcAlphaComp) + 0.5f);
    uint8_t resultG = static_cast<uint8_t>((srcG * dstAlpha) + (dstG * srcAlphaComp) + 0.5f);
    uint8_t resultB = static_cast<uint8_t>((srcB * dstAlpha) + (dstB * srcAlphaComp) + 0.5f);

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}

GPixel blendDstTop(GPixel src, GPixel dst) {
    uint8_t srcA = GPixel_GetA(src);
    uint8_t srcR = GPixel_GetR(src);
    uint8_t srcG = GPixel_GetG(src);
    uint8_t srcB = GPixel_GetB(src);

    uint8_t dstA = GPixel_GetA(dst);
    uint8_t dstR = GPixel_GetR(dst);
    uint8_t dstG = GPixel_GetG(dst);
    uint8_t dstB = GPixel_GetB(dst);

    float srcAlpha = srcA / 255.0f;
    float dstAlphaComp = 1.0f - (dstA / 255.0f);

    uint8_t resultA = static_cast<uint8_t>(std::min(255.0f, (dstA * srcAlpha) + (srcA * dstAlphaComp) + 0.5f));
    uint8_t resultR = static_cast<uint8_t>((dstR * srcAlpha) + (srcR * dstAlphaComp) + 0.5f);
    uint8_t resultG = static_cast<uint8_t>((dstG * srcAlpha) + (srcG * dstAlphaComp) + 0.5f);
    uint8_t resultB = static_cast<uint8_t>((dstB * srcAlpha) + (srcB * dstAlphaComp) + 0.5f);

    return GPixel_PackARGB(resultA, resultR, resultG, resultB);
}

BlendFunc getBlendFunction(GBlendMode blendMode) {
    switch (blendMode) {
        case GBlendMode::kClear: return clear;
        case GBlendMode::kSrcATop: return blendSrcTop;
        case GBlendMode::kDstATop: return blendDstTop;
        case GBlendMode::kSrc: return blendSrc;
        case GBlendMode::kDst: return blendDst;
        case GBlendMode::kSrcOut: return blendSrcOut;
        case GBlendMode::kDstOut: return blendDstOut;
        case GBlendMode::kSrcOver: return blendSrcOver;
        case GBlendMode::kDstOver: return blendDstOver;
        case GBlendMode::kSrcIn: return blendSrcIn;
        case GBlendMode::kDstIn: return blendDstIn;
        default: return blendSrcOver; 
    }
}