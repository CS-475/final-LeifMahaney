#include "include/GMatrix.h"
#include <cmath>
#include "include/GPathBuilder.h"


inline void GPathBuilder::reset() {
    fPts.clear(); 
    fVbs.clear(); 
}
inline void GPathBuilder::transform(const GMatrix& matrix) {
    for (size_t i = 0; i < fPts.size(); ++i) {
        fPts[i] = matrix * fPts[i];
    }
}

inline std::shared_ptr<GPath> GPathBuilder::detach() {
    auto path = std::make_shared<GPath>(std::move(fPts), std::move(fVbs));
    reset();
    return path;
}


void GPathBuilder::addPolygon(const GPoint pts[], int count) {
    assert(count > 0); 
    moveTo(pts[0]); 
    for (int i = 1; i < count; ++i) {
        lineTo(pts[i]);
    }
}

constexpr float kCircleConstant = 0.551915;  

void GPathBuilder::addCircle(GPoint center, float radius, GPathDirection dir) {
    float off = kCircleConstant * radius;

    GPoint srcPoints[12] = {
        {1, 0},                
        {1, off},           
        {off, 1},           
        {0, 1},                
        {-off, 1},          
        {-1, off},          
        {-1, 0},               
        {-1, -off},         
        {-off, -1},         
        {0, -1},               
        {off, -1},          
        {1, -off}          
    };

    GPoint dstPoints[12];

    GMatrix matrix = GMatrix::Translate(center.x, center.y) * GMatrix::Scale(radius, radius);
    matrix.mapPoints(dstPoints, srcPoints, 12);

    int indices[12];
    if (dir == GPathDirection::kCW) {
        indices[0] = 0; indices[1] = 1; indices[2] = 2;
        indices[3] = 3; indices[4] = 4; indices[5] = 5;
        indices[6] = 6; indices[7] = 7; indices[8] = 8;
        indices[9] = 9; indices[10] = 10; indices[11] = 11;
    } else {
        indices[0] = 0; indices[1] = 11; indices[2] = 10;
        indices[3] = 9; indices[4] = 8; indices[5] = 7;
        indices[6] = 6; indices[7] = 5; indices[8] = 4;
        indices[9] = 3; indices[10] = 2; indices[11] = 1;
    }

    moveTo(dstPoints[indices[0]]);

    for (int i = 0; i < 12; i += 3) {
        quadTo(dstPoints[indices[i + 1]], dstPoints[indices[i + 2]]);
    }
}
void GPathBuilder::addRect(const GRect& rect, GPathDirection dir) {
    GPoint topL = {rect.left, rect.top};
    GPoint topR = {rect.right, rect.top};
    GPoint bottomR = {rect.right, rect.bottom};
    GPoint bottomL = {rect.left, rect.bottom};

    if (dir == GPathDirection::kCW) {
        moveTo(topL);
        lineTo(topR);
        lineTo(bottomR);
        lineTo(bottomL);
        lineTo(topL); 
    } else { 
        moveTo(topL);
        lineTo(bottomL);
        lineTo(bottomR);
        lineTo(topR);
        lineTo(topL); 
    }
}

