#ifndef MYCANVAS_H
#define MYCANVAS_H
#include "include/GMatrix.h"
#include "include/GShader.h"
#include "include/GCanvas.h"
#include "include/GPaint.h"
#include "include/GBitmap.h"
#include "include/GRect.h"
#include <stack>

typedef GPixel(*BlendProc)(GPixel src, GPixel dst);
BlendProc find_blend(GBlendMode blendMode);


class MyCanvas : public GCanvas {
public:
    MyCanvas(GBitmap& bitmap);
    ~MyCanvas() override;

    void clear(const GColor& color) override;
    void concat(const GMatrix& matrix) override;
    void drawRect(const GRect& rect, const GPaint& paint) override;
    void save() override;
    void restore() override;
    void drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) override;
    void drawPath(const GPath& path, const GPaint& paint) override;
    void drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint& paint) override;
    void drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4], int level, const GPaint& paint) override;




private:
    GBitmap& fBitmap; 
    std::stack<GMatrix> fMatrixStack;
    GMatrix fMatrix;

    GPixel blendPixels(GPixel dstPixel, GPixel srcPixel, GBlendMode blendMode) const;
    GPixel colorsToPixels(const GColor& color) const;
    void computeBarycentric(int x, int y, const GPoint& p0, const GPoint& p1, const GPoint& p2, float& alpha, float& beta, float& gamma);
    void drawTriangleInline(const GPoint verts[3], const GColor colors[3], const GPoint texs[3],const GPaint& paint, const GBitmap& device);

};

#endif