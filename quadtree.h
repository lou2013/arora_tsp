#ifndef QUADTREE_H
#define QUADTREE_H

#include <vector>
#include "geometry.h"

class Quadtree {
public:
    Quadtree(BoundingSquare bounds, int capacity = 4);
    void Insert(const Point& point);
    std::vector<Point> QueryRange(const BoundingSquare& range) const;
    bool IsSubdivided() const;
    void Subdivide();

    // Get the points in this quadtree node
    const std::vector<Point>& GetPoints() const;

private:
    BoundingSquare bounds;   // Bounding square for this node
    int capacity;            // Maximum number of points per node
    std::vector<Point> points;  // Points stored in this node
    Quadtree* NW;            // North-West quadrant
    Quadtree* NE;            // North-East quadrant
    Quadtree* SW;            // South-West quadrant
    Quadtree* SE;            // South-East quadrant

    void Split();  // Helper method to split a node

};

#endif // QUADTREE_H
