#include "perturbation.h"
#include <algorithm>

Perturbation::Perturbation(std::vector<Point>& points, double grid_size)
    : points(points), grid_size(grid_size) {
    // Calculate the bounding box for the points
    if (!points.empty()) {
        double minX = points[0].x, minY = points[0].y;
        double maxX = points[0].x, maxY = points[0].y;

        for (const auto& p : points) {
            minX = std::min(minX, p.x);
            minY = std::min(minY, p.y);
            maxX = std::max(maxX, p.x);
            maxY = std::max(maxY, p.y);
        }

        boundingSquare = BoundingSquare(minX, minY, maxX, maxY);
    }
}

void Perturbation::Apply() {
    for (auto& p : points) {
        // Snapping the points to the nearest grid point
        p.x = grid_size * std::round(p.x / grid_size);
        p.y = grid_size * std::round(p.y / grid_size);
    }
}

BoundingSquare Perturbation::GetBoundingSquare() const {
    return boundingSquare;
}

const std::vector<Point>& Perturbation::GetPoints() const {
    return points;
}
