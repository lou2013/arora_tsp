#ifndef PORTALS_H
#define PORTALS_H

#include "geometry.h"
#include <vector>

class PortalGenerator {
public:
    PortalGenerator(int numPortalsPerSide);

    // Generates portals for a given square
    std::vector<Point> GeneratePortals(const BoundingSquare& square) const;

private:
    int k; // Number of portals per edge
};

#endif // PORTALS_H
