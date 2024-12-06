
#include <algorithm> 
#include "include/GRect.h"
#include "include/GPath.h"
#include "include/GMatrix.h"


inline std::shared_ptr<GPath> GPath::transform(const GMatrix& matrix) const {
    std::vector<GPoint> transformedPts(fPts.size());
    matrix.mapPoints(transformedPts.data(), fPts.data(), fPts.size());
    return std::make_shared<GPath>(std::move(transformedPts), fVbs);
}

inline GPath::Iter::Iter(const GPath& path)
    : fCurrPt(path.fPts.data())
    , fCurrVb(path.fVbs.data())
    , fStopVb(path.fVbs.data() + path.fVbs.size())
{}


void updateQuadBounds(const GPoint pts[3], float& minX, float& minY, float& maxX, float& maxY) {
    const float epsilon = 1e-6f;

    minX = std::min({minX, pts[0].x, pts[2].x});
    minY = std::min({minY, pts[0].y, pts[2].y});
    maxX = std::max({maxX, pts[0].x, pts[2].x});
    maxY = std::max({maxY, pts[0].y, pts[2].y});

    float d = pts[0].x - 2 * pts[1].x + pts[2].x;
    if (std::abs(d) > epsilon) {
        float tx = (pts[0].x - pts[1].x) / d;
        if (tx >= 0 && tx <= 1) {
            float x = (1 - tx) * (1 - tx) * pts[0].x + 2 * (1 - tx) * tx * pts[1].x + tx * tx * pts[2].x;
            minX = std::min(minX, x);
            maxX = std::max(maxX, x);
        }
    }

    d = pts[0].y - 2 * pts[1].y + pts[2].y;
    if (std::abs(d) > epsilon) {
        float ty = (pts[0].y - pts[1].y) / d;
        if (ty >= 0 && ty <= 1) {
            float y = (1 - ty) * (1 - ty) * pts[0].y + 2 * (1 - ty) * ty * pts[1].y + ty * ty * pts[2].y;
            minY = std::min(minY, y);
            maxY = std::max(maxY, y);
        }
    }
}


void updateCubicBounds(const GPoint pts[4], float& minX, float& minY, float& maxX, float& maxY) {
    const float epsilon = 1e-6f;

    minX = std::min({minX, pts[0].x, pts[3].x});
    minY = std::min({minY, pts[0].y, pts[3].y});
    maxX = std::max({maxX, pts[0].x, pts[3].x});
    maxY = std::max({maxY, pts[0].y, pts[3].y});

    float ax = -pts[0].x + 3 * pts[1].x - 3 * pts[2].x + pts[3].x;
    float bx = 3 * pts[0].x - 6 * pts[1].x + 3 * pts[2].x;
    float cx = -3 * pts[0].x + 3 * pts[1].x;

    std::vector<float> tValues;
    float discrim = bx * bx - 4 * ax * cx;
    if (discrim >= 0 && std::abs(ax) > epsilon) {
        float sqrtD = sqrt(discrim);
        tValues.push_back((-bx + sqrtD) / (2 * ax));
        tValues.push_back((-bx - sqrtD) / (2 * ax));
    }

    for (float t : tValues) {
        if (t >= 0 && t <= 1) {
            float x = (1 - t) * (1 - t) * (1 - t) * pts[0].x +
                      3 * (1 - t) * (1 - t) * t * pts[1].x +
                      3 * (1 - t) * t * t * pts[2].x +
                      t * t * t * pts[3].x;
            minX = std::min(minX, x);
            maxX = std::max(maxX, x);
        }
    }

    float ay = -pts[0].y + 3 * pts[1].y - 3 * pts[2].y + pts[3].y;
    float by = 3 * pts[0].y - 6 * pts[1].y + 3 * pts[2].y;
    float cy = -3 * pts[0].y + 3 * pts[1].y;

    discrim = by * by - 4 * ay * cy;
    if (discrim >= 0 && std::abs(ay) > epsilon) {
        float sqrtD = sqrt(discrim);
        tValues.clear();
        tValues.push_back((-by + sqrtD) / (2 * ay));
        tValues.push_back((-by - sqrtD) / (2 * ay));
    }

    for (float t : tValues) {
        if (t >= 0 && t <= 1) {
            float y = (1 - t) * (1 - t) * (1 - t) * pts[0].y +
                      3 * (1 - t) * (1 - t) * t * pts[1].y +
                      3 * (1 - t) * t * t * pts[2].y +
                      t * t * t * pts[3].y;
            minY = std::min(minY, y);
            maxY = std::max(maxY, y);
        }
    }
}

inline GPath::Edger::Edger(const GPath& path)
    : fPrevMove(nullptr)
    , fCurrPt(path.fPts.data())
    , fCurrVb(path.fVbs.data())
    , fStopVb(path.fVbs.data() + path.fVbs.size())
    , fPrevVerb(-1)
{}

GRect GPath::bounds() const {
    if (fPts.empty()) return GRect::XYWH(0, 0, 0, 0);

    float minX = fPts[0].x;
    float minY = fPts[0].y;
    float maxX = fPts[0].x;
    float maxY = fPts[0].y;

    GPath::Edger edger(*this);
    GPoint pts[4];  // Holds points for each segment

    while (auto verb = edger.next(pts)) {
        switch (*verb) {
            case GPathVerb::kMove:
            case GPathVerb::kLine:
                minX = std::min({minX, pts[0].x, pts[1].x});
                minY = std::min({minY, pts[0].y, pts[1].y});
                maxX = std::max({maxX, pts[0].x, pts[1].x});
                maxY = std::max({maxY, pts[0].y, pts[1].y});
                break;
            case GPathVerb::kQuad:
                minX = std::min({minX, pts[0].x, pts[1].x, pts[2].x});
                minY = std::min({minY, pts[0].y, pts[1].y, pts[2].y});
                maxX = std::max({maxX, pts[0].x, pts[1].x, pts[2].x});
                maxY = std::max({maxY, pts[0].y, pts[1].y, pts[2].y});
                break;
            case GPathVerb::kCubic:
                minX = std::min({minX, pts[0].x, pts[1].x, pts[2].x, pts[3].x});
                minY = std::min({minY, pts[0].y, pts[1].y, pts[2].y, pts[3].y});
                maxX = std::max({maxX, pts[0].x, pts[1].x, pts[2].x, pts[3].x});
                maxY = std::max({maxY, pts[0].y, pts[1].y, pts[2].y, pts[3].y});
                break;
            default:
                break;
        }
    }

    return GRect::LTRB(minX, minY, maxX, maxY);
}


void GPath::ChopCubicAt(const GPoint src[4], GPoint dst[7], float t) {
    GPoint p01 = (1 - t) * src[0] + t * src[1]; 
    GPoint p14 = (1 - t) * src[1] + t * src[2];
    GPoint p235 = (1 - t) * src[2] + t * src[3];
    
    GPoint p012 = (1 - t) * p01 + t * p14;      
    GPoint p1234 = (1 - t) * p14 + t * p235;      
    
    GPoint p01234 = (1 - t) * p012 + t * p1234;  
    
    dst[0] = src[0];
    dst[1] = p01;
    dst[2] = p012;
    dst[3] = p01234;
    
    dst[4] = p1234;
    dst[5] = p235;
    dst[6] = src[3];
}
inline nonstd::optional<GPathVerb> GPath::Edger::next(GPoint pts[]) {
    while (fCurrVb < fStopVb) {
        GPathVerb verb = *fCurrVb++;
        fPrevVerb = static_cast<int>(verb);
        
        switch (verb) {
            case kMove:
                fPrevMove = fCurrPt; 
                pts[0] = *fCurrPt++;
                return verb;
            case kLine: {
                if (fPrevMove) { 
                    pts[0] = *fPrevMove;
                    pts[1] = *fCurrPt++;
                    return verb;
                }
                break;
            }
            default:
                break;
        }
    }
    return nonstd::nullopt; 
}

void GPath::ChopQuadAt(const GPoint src[3], GPoint dst[5], float t) {
    GPoint p01 = (1 - t) * src[0] + t * src[1]; 
    GPoint p1234 = (1 - t) * src[1] + t * src[2];
    GPoint p01234 = (1 - t) * p01 + t * p1234;     
    
    dst[0] = src[0];
    dst[1] = p01;
    dst[2] = p01234;
    
    dst[3] = p1234;
    dst[4] = src[2];
}


inline nonstd::optional<GPathVerb> GPath::Iter::next(GPoint pts[]) {
    if (fCurrVb < fStopVb) {
        GPathVerb verb = *fCurrVb++;
        switch (verb) {
            case kMove:
                pts[0] = *fCurrPt++;
                return verb;
            case kLine:
                pts[0] = *(fCurrPt++);
                pts[1] = *(fCurrPt++);
                return verb;
            case kQuad:
                pts[0] = *(fCurrPt++);
                pts[1] = *(fCurrPt++);
                pts[2] = *(fCurrPt++);
                return verb;
            case kCubic:
                pts[0] = *(fCurrPt++);
                pts[1] = *(fCurrPt++);
                pts[2] = *(fCurrPt++);
                pts[3] = *(fCurrPt++);
                return verb;
            default:
                break;
        }
    }
    return nonstd::nullopt;
}



