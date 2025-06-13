#include "portals.h"

PortalGenerator::PortalGenerator(int numPortalsPerSide) : k(numPortalsPerSide) {}
std::vector<Point> PortalGenerator::GeneratePortals(const BoundingSquare& square) const {
    std::vector<Point> portals;

    double dx = (square.maxX - square.minX) / (k + 1);  // Use this for spacing
    double dy = (square.maxY - square.minY) / (k + 1);

    // Top edge (minY)
    for (int i = 1; i <= k; ++i) {
        double x = square.minX + i * dx;  // Position portals correctly along X-axis
        portals.emplace_back(x, square.minY);
    }

    // Bottom edge (maxY)
    for (int i = 1; i <= k; ++i) {
        double x = square.minX + i * dx;  // Position portals correctly along X-axis
        portals.emplace_back(x, square.maxY);
    }

    // Left edge (minX)
    for (int i = 1; i <= k; ++i) {
        double y = square.minY + i * dy;  // Position portals correctly along Y-axis
        portals.emplace_back(square.minX, y);
    }

    // Right edge (maxX)
    for (int i = 1; i <= k; ++i) {
        double y = square.minY + i * dy;  // Position portals correctly along Y-axis
        portals.emplace_back(square.maxX, y);
    }

    return portals;
}
