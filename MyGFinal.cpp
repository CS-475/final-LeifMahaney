#include "MyGFinal.h"
#include "include/GPathBuilder.h"
#include <cmath>



MyGFinal::MyGFinal() = default;
std::unique_ptr<GFinal> GCreateFinal() {
    return std::make_unique<MyGFinal>();
}

std::shared_ptr<GShader> MyGFinal::createLinearPosGradient(GPoint p0, GPoint p1, const GColor colors[], const float pos[], int count) {
    if (count < 1 || !colors || !pos) {
        return nullptr;
    }

    std::vector<GColor> interpolatedColors;
    
    for (float t = 0; t <= 1.0f; t += 0.01f) { 
        int index = 0;
        while (index < count - 1 && t > pos[index + 1]) {
            index++;
        }

        GColor color;
        if (t <= pos[0]) {
            color = colors[0];
        } else if (t >= pos[count - 1]) {
            color = colors[count - 1];
        } else {
            float t00 = pos[index];
            float t11 = pos[index + 1];
            float weighted = (t - t00) / (t11 - t00);
            const GColor& c0 = colors[index];
            const GColor& c1 = colors[index + 1];
            color = {
                c0.r + (c1.r - c0.r) * weighted,
                c0.g + (c1.g - c0.g) * weighted,
                c0.b + (c1.b - c0.b) * weighted,
                c0.a + (c1.a - c0.a) * weighted
            };
        }
        interpolatedColors.push_back(color);
    }
    return GCreateLinearGradient(p0, p1, interpolatedColors.data(), interpolatedColors.size(), GTileMode::kClamp);
}


std::shared_ptr<GShader> MyGFinal::createVoronoiShader(const GPoint points[], const GColor colors[], int count) {
    return nullptr;
}

std::shared_ptr<GPath> MyGFinal::strokePolygon(const GPoint points[], int count, float width, bool isClosed) {
    if (count < 2) return nullptr;
    GPathBuilder builder;
    float halfWidth = width/2;
    auto getNorm = [](GPoint p0, GPoint p1, float scale) -> std::pair<float, float> {
        float dx = p1.x - p0.x;
        float dy = p1.y - p0.y;
        float length = sqrt(dx*dx + dy*dy);
        if (length < 0.0001f) {
            return {0, scale};
        }
        return {(-dy / length) * scale, (dx / length) * scale};
    };
    auto addRJ = [&builder, halfWidth](GPoint center, float nx1, float ny1, float nx2, float ny2) {
        float angle1 = atan2(ny1, nx1);
        float angle2 = atan2(ny2, nx2);
        if (fabs(angle1 - angle2) < 0.0001f) {
            builder.lineTo({center.x + nx2, center.y + ny2});
            return;
        }
        if (angle1 - angle2 > M_PI) angle1 -= 2 * M_PI;
        if (angle2 - angle1 > M_PI) angle2 -= 2 * M_PI;
        float angleDifference = fabs(angle2 - angle1);
        int totalsteps = std::max(2, int(8 * angleDifference / M_PI));
        for (float t = 0; t <= 1; t += 1.0f/totalsteps) {
            float angle = angle1 + (angle2 - angle1) * t;
            float x = center.x + cos(angle) * halfWidth;
            float y = center.y + sin(angle) * halfWidth;
            builder.lineTo({x, y});
        }
    };
    auto validateNormal = [](float nx, float ny, float scale) -> std::pair<float, float> {
        float length = sqrt(nx*nx + ny*ny);
        if (length < 0.0001f) return {0, scale};
        float factor = scale / length;
        return {nx * factor, ny * factor};
    };

    auto [nx0, ny0] = getNorm(points[0], points[1], halfWidth);

    builder.moveTo({points[0].x + nx0, points[0].y + ny0});
    for (int i = 1; i < count; i++) {
        auto [nx_curr, ny_curr] = getNorm(points[i-1], points[i], halfWidth);
        
        if (i < count-1) {
            auto [nx_next, ny_next] = getNorm(points[i], points[i+1], halfWidth);
            addRJ(points[i], nx_curr, ny_curr, nx_next, ny_next);
        } else {
            builder.lineTo({points[i].x + nx_curr, points[i].y + ny_curr});
        }
    }


    if (isClosed) {
        auto [nx_last, ny_last] = getNorm(points[count-1], points[0], halfWidth);
        auto [nx_first, ny_first] = getNorm(points[0], points[1], halfWidth);
        addRJ(points[0], nx_last, ny_last, nx_first, ny_first);
        builder.lineTo({points[0].x + nx_first, points[0].y + ny_first});
    } else {

        auto [nx_last, ny_last] = getNorm(points[count-2], points[count-1], halfWidth);
        builder.lineTo({points[count-1].x + nx_last, points[count-1].y + ny_last});
    }

    for (int i = count-1; i >= 0; i--) {
        auto [nx_curr, ny_curr] = getNorm(points[i], points[(i + 1) % count], halfWidth);
        
        if (i > 0) {
            auto [nx_prev, ny_prev] = getNorm(points[i-1], points[i], halfWidth);
            addRJ(points[i], -nx_curr, -ny_curr, -nx_prev, -ny_prev);
        } else {
            builder.lineTo({points[i].x - nx_curr, points[i].y - ny_curr});
        }
    }


    builder.lineTo({points[0].x + nx0, points[0].y + ny0});

    return builder.detach();
}