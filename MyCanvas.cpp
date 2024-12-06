#include "MyCanvas.h"
#include "include/GPath.h"
#include "include/GPathBuilder.h"
#include "BlendFunctions.h"
#include "include/GMatrix.h"
#include "include/GShader.h"
#include "include/GPoint.h"
#include "include/GRect.h"
#include "include/GPixel.h"
#include "include/GBitmap.h"
#include "include/GBlendMode.h"
#include "include/GPaint.h"
#include "include/GCanvas.h"


MyCanvas::MyCanvas(GBitmap& bitmap) : fBitmap(bitmap) {
}

MyCanvas::~MyCanvas() {
}

GPixel MyCanvas::colorsToPixels(const GColor& color) const {
    uint8_t a = static_cast<uint8_t>(color.a * 255 + 0.5f);
    uint8_t r = static_cast<uint8_t>(color.r * color.a * 255 + 0.5f);
    uint8_t g = static_cast<uint8_t>(color.g * color.a * 255 + 0.5f);
    uint8_t b = static_cast<uint8_t>(color.b * color.a * 255 + 0.5f);
    return GPixel_PackARGB(a, r, g, b);
}


GPixel MyCanvas::blendPixels(GPixel dstPixel, GPixel srcPixel, GBlendMode blendMode) const {
    BlendFunc blendFunc = getBlendFunction(blendMode);
    return blendFunc(srcPixel, dstPixel);
}

void MyCanvas::clear(const GColor& color) {
    int height = fBitmap.height();
    int width = fBitmap.width();
    GPixel p = colorsToPixels(color);
    GPixel *row_addr = nullptr;
    for (int y = 0; y < height; y++) {
        row_addr = fBitmap.getAddr(0, y);
        for (int x = 0; x < width; x++) {
            row_addr[x] = p;
        }
    }
}
static inline int GDiv255(int value) {
    return (value + 128) * 257 >> 16;
}



void MyCanvas::drawRect(const GRect& rect, const GPaint& paint){
    GPoint points[4] = {
        {rect.left, rect.top},
        {rect.right, rect.top},
        {rect.right, rect.bottom},
        {rect.left, rect.bottom}
    };
    GPoint transformedPoints[4];
    for(int i = 0; i < 4; ++i){
        float fx = points[i].x;
        float fy = points[i].y;
        transformedPoints[i].x = fMatrix[0] * fx + fMatrix[2] * fy + fMatrix[4];
        transformedPoints[i].y = fMatrix[1] * fx + fMatrix[3] * fy + fMatrix[5];
    }

    GIRect irs = GIRect::LTRB(
        static_cast<int32_t>(std::min({transformedPoints[0].x, transformedPoints[1].x, transformedPoints[2].x, transformedPoints[3].x})),
        static_cast<int32_t>(std::min({transformedPoints[0].y, transformedPoints[1].y, transformedPoints[2].y, transformedPoints[3].y})),
        static_cast<int32_t>(std::max({transformedPoints[0].x, transformedPoints[1].x, transformedPoints[2].x, transformedPoints[3].x})),
        static_cast<int32_t>(std::max({transformedPoints[0].y, transformedPoints[1].y, transformedPoints[2].y, transformedPoints[3].y}))
    );
    int left = std::max(0, irs.left);
    int top = std::max(0, irs.top);
    int right = std::min(fBitmap.width(), irs.right);
    int bottom = std::min(fBitmap.height(), irs.bottom);

    if(left >= right || top >= bottom){
        return;
    }

    GShader* shade = paint.peekShader();
    GPixel srPixel;
    bool useShade = (shade != nullptr && shade-> setContext(fMatrix));
    for(int y = top; y < bottom; ++y){
        GPixel* rowAddr = fBitmap.getAddr(0, y);
        for(int x = left ; x < right; ++x){
            if(useShade){
                shade->shadeRow(x, y, 1, &srPixel);
            }else{
                srPixel = colorsToPixels(paint.getColor());
            }
            if(paint.getBlendMode() == GBlendMode::kSrc){
                rowAddr[x] = srPixel;
            }else{
                rowAddr[x] = blendPixels(rowAddr[x], srPixel, paint.getBlendMode());
            }
        }
    }
}


struct Edge {
    float m;    
    float b;   
    int top;  
    int bottom;
    int xL; 
    int xR; 
    int currentX; 
    int windings;

    bool isValid(int y) const {
        return (y >= top && y < bottom);  
    }

    float computeX(int y) const {
        return m * y + b;  
    }

    bool isUseful() const {
        return top < bottom;
    }
};

Edge clipE(const Edge& edge, int canvasWidth, int canvasHeight) {
    Edge clippedEdge = edge;


    if (clippedEdge.top < 0) {
        clippedEdge.xL = static_cast<int>(std::round(clippedEdge.m * -clippedEdge.top + clippedEdge.b));
        clippedEdge.top = 0;
    }


    if (clippedEdge.bottom > canvasHeight) {
        clippedEdge.bottom = canvasHeight;
    }


    if (clippedEdge.xL < 0) {
        clippedEdge.xL = 0;
        clippedEdge.top = static_cast<int>(std::round((0 - clippedEdge.b) / clippedEdge.m));
    }


    if (clippedEdge.xR > canvasWidth) {
        clippedEdge.xR = canvasWidth;
        clippedEdge.bottom = static_cast<int>(std::round((canvasWidth - clippedEdge.b) / clippedEdge.m));
    }

    return clippedEdge;
}

Edge makeE(const GPoint& p0, const GPoint& p1) {
    GPoint top = p0, bottom = p1;
    int windings = 1;
    if (p0.y > p1.y) {
        std::swap(top, bottom);
        windings = -1;
    }

    Edge edge;
    edge.m = (bottom.x - top.x) / (bottom.y - top.y);  
    edge.b = top.x - edge.m * top.y; 
    edge.top = GRoundToInt(top.y);
    edge.bottom = GRoundToInt(bottom.y);
    edge.currentX = GRoundToInt(top.x); 
    edge.windings = windings;
    return edge;
}

void MyCanvas::drawConvexPolygon(const GPoint points[], int count, const GPaint& paint) {
    GPoint transformedPoints[count];
    for (int i = 0; i < count; ++i) {
        float sx = points[i].x;
        float sy = points[i].y;
        transformedPoints[i].x = fMatrix[0] * sx + fMatrix[2] * sy + fMatrix[4];
        transformedPoints[i].y = fMatrix[1] * sx + fMatrix[3] * sy + fMatrix[5];
    }

    std::vector<Edge> edges;
    for (int i = 0; i < count; ++i) {
        int next = (i + 1) % count;
        Edge edge = makeE(transformedPoints[i], transformedPoints[next]);
        if (edge.isUseful()) {
            edges.push_back(edge);
        }
    }

    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        return a.top < b.top;
    });

    int canvasHeight = fBitmap.height();
    GShader* shader = paint.peekShader();
    bool useShader = (shader != nullptr && shader->setContext(fMatrix));
    GPixel srcPixel;

    for (int y = 0; y < canvasHeight; ++y) {
        std::vector<int> intersections;

        for (const auto& edge : edges) {
            if (y >= edge.top && y <= edge.bottom) {
                int x = static_cast<int>(std::round(edge.m * y + edge.b));
                intersections.push_back(x);
            }
        }

        std::sort(intersections.begin(), intersections.end());

        if (useShader) {
            shader->shadeRow(0, y, 1, &srcPixel); 
        } else {
            srcPixel = colorsToPixels(paint.getColor());
        }

        GPixel* rowAddr = fBitmap.getAddr(0, y);
        for (size_t i = 0; i < intersections.size(); i += 2) {
            if (i + 1 < intersections.size()) {
                int xStart = std::max(0, intersections[i]);
                int xEnd = std::min(fBitmap.width(), intersections[i + 1]);

                if (xStart < xEnd) {
                    for (int x = xStart; x < xEnd; ++x) {
                        if (useShader) {
                            shader->shadeRow(x, y, 1, &srcPixel);
                        }
                        rowAddr[x] = blendPixels(rowAddr[x], srcPixel, paint.getBlendMode());
                    }
                }
            }
        }
    }
}


void MyCanvas::save() {
    fMatrixStack.push(fMatrix);
}

void MyCanvas::restore() {
    if (!fMatrixStack.empty()) {
        fMatrix = fMatrixStack.top();
        fMatrixStack.pop();
    }
}

void MyCanvas::concat(const GMatrix& matrix) {
    fMatrix = GMatrix::Concat(fMatrix, matrix);
}

int computeQuadSegments(const GPoint pts[3], float tolerance) {
    float ax = pts[0].x - 2 * pts[1].x + pts[2].x;
    float ay = pts[0].y - 2 * pts[1].y + pts[2].y;
    float maxDist = std::sqrt(ax * ax + ay * ay);
    return static_cast<int>(std::ceil(std::sqrt(maxDist / tolerance)));

}

int computeCubicSegments(const GPoint pts[4], float tolerance) {
    float ax = -pts[0].x + 3 * (pts[1].x - pts[2].x) + pts[3].x;
    float ay = -pts[0].y + 3 * (pts[1].y - pts[2].y) + pts[3].y;
    float maxDist = std::sqrt(ax * ax + ay * ay);
    return static_cast<int>(std::ceil(std::sqrt(maxDist / tolerance)));
}

GPoint evalQuad(const GPoint pts[3], float t) {
    float u = 1 - t;
    return {
        u * u * pts[0].x + 2 * u * t * pts[1].x + t * t * pts[2].x,
        u * u * pts[0].y + 2 * u * t * pts[1].y + t * t * pts[2].y
    };
}

GPoint evalCubic(const GPoint pts[4], float t) {
    float u = 1 - t;
    return {
        u * u * u * pts[0].x + 3 * u * u * t * pts[1].x + 3 * u * t * t * pts[2].x + t * t * t * pts[3].x,
        u * u * u * pts[0].y + 3 * u * u * t * pts[1].y + 3 * u * t * t * pts[2].y + t * t * t * pts[3].y
    };
}

void MyCanvas::drawPath(const GPath& path, const GPaint& paint) {
    GPathBuilder builder;
    builder.reset();
    GPath::Iter iter(path);
    GPoint pts[GPath::kMaxNextPoints];
    const float tolerance = 0.25f;

    while (auto verb = iter.next(pts)) {
        switch (verb.value()) {
            case GPathVerb::kMove:
                builder.moveTo(pts[0]);
                break;
            case GPathVerb::kLine:
                builder.lineTo(pts[1]);
                break;
            case GPathVerb::kQuad: {
                int seg = computeQuadSegments(pts, tolerance);
                for (int i = 0; i < seg; ++i) {
                    float t1 = static_cast<float>(i + 1) / seg;
                    builder.lineTo(evalQuad(pts, t1));
                }
                break;
            }
            case GPathVerb::kCubic: {
                int seg = computeCubicSegments(pts, tolerance);
                for (int i = 0; i < seg; ++i) {
                    float t1 = static_cast<float>(i + 1) / seg;
                    builder.lineTo(evalCubic(pts, t1));
                }
                break;
            }
            default:
                break;
        }
    }

    builder.transform(fMatrix);
    auto transformedPath = builder.detach();
    
    GPath::Edger edger(*transformedPath);
    std::vector<Edge> edges;
    while (auto verb = edger.next(pts)) {
        if (verb.value() == GPathVerb::kLine) {
            Edge edge = makeE(pts[0], pts[1]);
            if (edge.isUseful()) {
                edges.push_back(edge);
            }
        }
    }
    
    std::sort(edges.begin(), edges.end(), [](const Edge& a, const Edge& b) {
        return a.top < b.top || (a.top == b.top && a.currentX < b.currentX);
    });

    int canvasHeight = fBitmap.height();
    int canvasWidth = fBitmap.width();
    GShader* shader = paint.peekShader();
    bool useShader = (shader != nullptr && shader->setContext(fMatrix));
    GPixel srcPixel;

    for (int y = 0; y < canvasHeight; ++y) {
        std::vector<int> xIntervals;
        int winding = 0;

        for (auto& edge : edges) {
            if (edge.isValid(y)) {
                int x = static_cast<int>(std::round(edge.computeX(y)));
                if (x < 0 || x >= canvasWidth) continue;
                xIntervals.push_back(x);
                winding += edge.windings;
            }
        }

        std::sort(xIntervals.begin(), xIntervals.end());

        for (size_t i = 0; i + 1 < xIntervals.size(); i += 2) {
            int left = xIntervals[i];
            int right = xIntervals[i + 1];
            left = std::max(0, left);
            right = std::min(canvasWidth, right);

            for (int xSpan = left; xSpan < right; ++xSpan) {
                if (useShader) {
                    shader->shadeRow(xSpan, y, 1, &srcPixel);
                } else {
                    srcPixel = colorsToPixels(paint.getColor());
                }

                GPixel* dst = fBitmap.getAddr(xSpan, y);
                *dst = blendPixels(*dst, srcPixel, paint.getBlendMode());
            }
        }
    }
}


std::unique_ptr<GCanvas> GCreateCanvas(const GBitmap& bitmap) {
    return std::make_unique<MyCanvas>(const_cast<GBitmap&>(bitmap));
}


std::string GDrawSomething(GCanvas* canvas, GISize dim) {
    canvas->clear(GColor::RGBA(0.7f, 0.9f, 0.5f, 1.0f));

    float CenterX = dim.width / 2.0f;
    float CenterY = dim.height / 2.0f;

    GBitmap bitmap;
    GMatrix localMatrix;

    std::shared_ptr<GShader> shader = GCreateBitmapShader(bitmap, localMatrix);

    GPoint hexagon[] = {
        {CenterX - 30, CenterY - 60},
        {CenterX + 30, CenterY - 50},
        {CenterX + 60, CenterY},
        {CenterX + 30, CenterY + 50},
        {CenterX - 30, CenterY + 50},
        {CenterX - 60, CenterY}
    };

    int hexCount = sizeof(hexagon) / sizeof(hexagon[0]);
    GPaint redPaint(GColor::RGBA(0.0f, 1.0f, 0.0f, 1.0f));
    redPaint.setShader(shader);
    canvas->drawConvexPolygon(hexagon, hexCount, redPaint);

    float octScale = 50.0f;

    GPoint octagon[] = {
        {CenterX - octScale, CenterY - 2 * octScale},
        {CenterX + octScale, CenterY - 2 * octScale},
        {CenterX + 2 * octScale, CenterY - octScale},
        {CenterX + 2 * octScale, CenterY + octScale},
        {CenterX + octScale, CenterY + 2 * octScale},
        {CenterX - octScale, CenterY + 2 * octScale},
        {CenterX - 2 * octScale, CenterY + octScale},
        {CenterX - 2 * octScale, CenterY - octScale},
    };
    int octCount = sizeof(octagon) / sizeof(octagon[0]);
    GPaint bluePaint(GColor::RGBA(0.5f, 0.0f, 0.5f, 1.0f));
    bluePaint.setShader(shader);
    canvas->drawConvexPolygon(octagon, octCount, bluePaint);

    GPoint parallelogram[] = {
        {CenterX - 60, CenterY + 80},
        {CenterX + 20, CenterY + 80},
        {CenterX + 60, CenterY + 140},
        {CenterX - 20, CenterY + 140}
    };
    int parCount = sizeof(parallelogram) / sizeof(parallelogram[0]);
    GPaint cyanPaint(GColor::RGBA(1.0f, 0.0f, 1.0f, 1.0f));
    cyanPaint.setShader(shader);
    canvas->drawConvexPolygon(parallelogram, parCount, cyanPaint);

    float pentagonScale = 0.9f;
    GPoint diamond[] = {
        {CenterX, CenterY - 40},
        {CenterX + 40, CenterY},
        {CenterX, CenterY + 40},
        {CenterX - 40 , CenterY},
    };

    int pentCount = sizeof(diamond) / sizeof(diamond[0]);
    GPaint orangePaint(GColor::RGBA(1.0f, 1.0f, 0.0f, 0.8f));
    orangePaint.setShader(shader);
    canvas->drawConvexPolygon(diamond, pentCount, orangePaint);

    return "Draw Something New";
}


void MyCanvas::drawMesh(const GPoint verts[], const GColor colors[], const GPoint texs[], int count, const int indices[], const GPaint& paint) {
    if (!paint.peekShader()) {
        texs = nullptr;
    }
    if (!colors && !texs) {
        return;
    }

    for (int i = 0; i < count; ++i) {
        int p0 = indices[3 * i];
        int p1 = indices[3 * i + 1];
        int p2 = indices[3 * i + 2];

        GPoint transformedVerts[3];
        GPoint originalVerts[3] = {verts[p0], verts[p1], verts[p2]};
        fMatrix.mapPoints(transformedVerts, originalVerts, 3);

        GColor triColors[3];
        if (colors) {
            triColors[0] = colors[p0];
            triColors[1] = colors[p1];
            triColors[2] = colors[p2];
        }

        GPoint triTexs[3];
        if (texs) {
            triTexs[0] = texs[p0];
            triTexs[1] = texs[p1];
            triTexs[2] = texs[p2];
        }

        drawTriangleInline(transformedVerts, 
                          colors ? triColors : nullptr,
                          texs ? triTexs : nullptr,
                          paint, fBitmap);
    }
}

void MyCanvas::drawQuad(const GPoint verts[4], const GColor colors[4], const GPoint texs[4],
                        int level, const GPaint& paint) {
    
    if (level < 1 || !verts) {
        return;
    }
    if (!colors && !texs) {
        return;
    }

    try {
        
        std::vector<GColor> generatedColors;
        std::vector<GPoint> generatedTexs;
        std::vector<GPoint> generatedVerts;
        std::vector<int> indices;

        int numVerts = (level + 1) * (level + 1);

        for (int i = 0; i <= level; ++i) {
            float v = static_cast<float>(i) / level;
            for (int j = 0; j <= level; ++j) {
                float u = static_cast<float>(j) / level;

                float x = verts[0].x * (1 - u) * (1 - v) +
                          verts[1].x * u * (1 - v) +
                          verts[3].x * (1 - u) * v +
                          verts[2].x * u * v;
                float y = verts[0].y * (1 - u) * (1 - v) +
                          verts[1].y * u * (1 - v) +
                          verts[3].y * (1 - u) * v +
                          verts[2].y * u * v;
                generatedVerts.push_back({x, y});

                if (colors) {
                    float c00 = (1 - u) * (1 - v);
                    float c10 = u * (1 - v);
                    float c01 = (1 - u) * v;
                    float c11 = u * v;
                    
                    GColor interpolatedColor;
                    interpolatedColor.a = c00 * colors[0].a + c10 * colors[1].a + c01 * colors[3].a + c11 * colors[2].a;
                    interpolatedColor.r = c00 * colors[0].r + c10 * colors[1].r + c01 * colors[3].r + c11 * colors[2].r;
                    interpolatedColor.g = c00 * colors[0].g + c10 * colors[1].g + c01 * colors[3].g + c11 * colors[2].g;
                    interpolatedColor.b = c00 * colors[0].b + c10 * colors[1].b + c01 * colors[3].b + c11 * colors[2].b;
                    generatedColors.push_back(interpolatedColor);
                }
                if (texs) {
                    float tx = (1 - u) * (1 - v) * texs[0].x +
                               u * (1 - v) * texs[1].x +
                               (1 - u) * v * texs[3].x +
                               u * v * texs[2].x;
                    float ty = (1 - u) * (1 - v) * texs[0].y +
                               u * (1 - v) * texs[1].y +
                               (1 - u) * v * texs[3].y +
                               u * v * texs[2].y;
                    generatedTexs.push_back({tx, ty});
                }
            }
        }

        for (int i = 0; i < level; ++i) {
            for (int j = 0; j < level; ++j) {
                int idx = i * (level + 1) + j;

                indices.push_back(idx);
                indices.push_back(idx + 1);
                indices.push_back(idx + level + 1);

                indices.push_back(idx + 1);
                indices.push_back(idx + level + 2);
                indices.push_back(idx + level + 1);
            }
        }

       
        drawMesh(generatedVerts.data(),
                colors ? generatedColors.data() : nullptr,
                texs ? generatedTexs.data() : nullptr,
                indices.size() / 3, indices.data(), paint);
        

    } catch (const std::exception& e) {    } catch (...) {
    }
}



void MyCanvas::drawTriangleInline(const GPoint verts[3], const GColor colors[3], const GPoint texs[3], 
                        const GPaint& paint, const GBitmap& device) {

    int left = std::max(0, static_cast<int>(std::floor(std::min({verts[0].x, verts[1].x, verts[2].x}))));
    int right = std::min(device.width(), static_cast<int>(std::ceil(std::max({verts[0].x, verts[1].x, verts[2].x}))));
    int top = std::max(0, static_cast<int>(std::floor(std::min({verts[0].y, verts[1].y, verts[2].y}))));
    int bottom = std::min(device.height(), static_cast<int>(std::ceil(std::max({verts[0].y, verts[1].y, verts[2].y}))));

    float ax1 = verts[1].x - verts[0].x;
    float ay1 = verts[1].y - verts[0].y;
    float ax2 = verts[2].x - verts[1].x;
    float ay2 = verts[2].y - verts[1].y;
    float ax3 = verts[0].x - verts[2].x;
    float ay3 = verts[0].y - verts[2].y;

    GShader* shader = paint.peekShader();
    bool useShader = shader && shader->setContext(GMatrix());

    float px = 0.5f;
    float py = 0.5f;

    for (int y = top; y < bottom; ++y) {
        GPixel* row = device.getAddr(0, y);
        float fy = y + py;

        for (int x = left; x < right; ++x) {
            float fx = x + px;

            float a0 = (ax2 * (fy - verts[1].y) - ay2 * (fx - verts[1].x));
            float a1 = (ax3 * (fy - verts[2].y) - ay3 * (fx - verts[2].x));
            float a2 = (ax1 * (fy - verts[0].y) - ay1 * (fx - verts[0].x));

            float area = ax2 * ay3 - ax3 * ay2;
            if (area == 0) continue;

            float iArea = 1.0f / area;
            a0 *= iArea;
            a1 *= iArea;
            a2 *= iArea;

            if (a0 >= 0 && a1 >= 0 && a2 >= 0) {
                GPixel srcPixel;

                if (useShader && texs) {
                    float tx = a0 * texs[0].x + a1 * texs[1].x + a2 * texs[2].x;
                    float ty = a0 * texs[0].y + a1 * texs[1].y + a2 * texs[2].y;
                    shader->shadeRow(tx, ty, 1, &srcPixel);
                } else if (colors) {
                    float a = a0 * colors[0].a + a1 * colors[1].a + a2 * colors[2].a;
                    float r = a0 * colors[0].r + a1 * colors[1].r + a2 * colors[2].r;
                    float g = a0 * colors[0].g + a1 * colors[1].g + a2 * colors[2].g;
                    float b = a0 * colors[0].b + a1 * colors[1].b + a2 * colors[2].b;

                    a = std::max(0.0f, std::min(1.0f, a));
                    r = std::max(0.0f, std::min(1.0f, r)) * a;
                    g = std::max(0.0f, std::min(1.0f, g)) * a;
                    b = std::max(0.0f, std::min(1.0f, b)) * a;

                    srcPixel = GPixel_PackARGB(
                        static_cast<uint8_t>(a * 255 + 0.5f),
                        static_cast<uint8_t>(r * 255 + 0.5f),
                        static_cast<uint8_t>(g * 255 + 0.5f),
                        static_cast<uint8_t>(b * 255 + 0.5f)
                    );
                }

                row[x] = blendPixels(row[x], srcPixel, paint.getBlendMode());
            }
        }
    }
}

void computeBarycentric(float x, float y, const GPoint& v0, const GPoint& v1, const GPoint& v2,
                        float& alpha, float& beta, float& gamma) {
    float denom = (v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y);
    alpha = ((v1.y - v2.y) * (x - v2.x) + (v2.x - v1.x) * (y - v2.y)) / denom;
    beta  = ((v2.y - v0.y) * (x - v2.x) + (v0.x - v2.x) * (y - v2.y)) / denom;
    gamma = 1.0f - alpha - beta;
}