#ifndef MY_G_FINAL_H
#define MY_G_FINAL_H
#include <memory>   
#include <vector>   
#include "include/GFinal.h" 
#include "include/GShader.h" 
#include "include/GPath.h"  

class MyGFinal : public GFinal {
public:
    MyGFinal();
    std::shared_ptr<GShader> createLinearPosGradient(GPoint p0, GPoint p1, const GColor colors[], const float pos[], int count) override;
    std::shared_ptr<GShader> createVoronoiShader(const GPoint points[], const GColor colors[], int count) override;
    std::shared_ptr<GPath> strokePolygon(const GPoint points[], int count, float width, bool isClosed) override;
};
std::unique_ptr<GFinal> GCreateFinal();
#endif
