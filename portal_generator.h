#pragma once

#include "geometry.h"
#include <vector>

class PortalGenerator {
public:
    // Generates k equally spaced portals on each of the 4 edges of the square.
    // The order is canonical: Top, Right, Bottom, Left.
    std::vector<Point> GeneratePortals(const BoundingSquare& square, int k) const {
        std::vector<Point> portals;
        if (k == 0) return portals;
        portals.reserve(4 * k);

        double dx = (square.maxX - square.minX) / (k + 1);
        double dy = (square.maxY - square.minY) / (k + 1);

        // Top edge (from left to right)
        for (int i = 1; i <= k; ++i) portals.push_back({square.minX + i * dx, square.maxY});
        // Right edge (from top to bottom)
        for (int i = 1; i <= k; ++i) portals.push_back({square.maxX, square.maxY - i * dy});
        // Bottom edge (from right to left)
        for (int i = 1; i <= k; ++i) portals.push_back({square.maxX - i * dx, square.minY});
        // Left edge (from bottom to top)
        for (int i = 1; i <= k; ++i) portals.push_back({square.minY + i * dy, square.minX});
        
        return portals;
    }
};