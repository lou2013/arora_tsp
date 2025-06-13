#pragma once
#ifndef PERTURBATION_H
#define PERTURBATION_H

#include <vector>
#include "geometry.h"

class Perturbation {
public:
    Perturbation(std::vector<Point>& points, double grid_size = 1.0);

    // Applies the perturbation (snapping points to grid)
    void Apply();

    // Returns the bounding square of the perturbed points
    BoundingSquare GetBoundingSquare() const;

    // Returns the perturbed points
    const std::vector<Point>& GetPoints() const;

private:
    std::vector<Point>& points;  // Reference to the original points
    double grid_size;            // Size of the grid cells
    BoundingSquare boundingSquare;  // Bounding square for the points
};

#endif // PERTURBATION_H
