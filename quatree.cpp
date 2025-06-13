#include "quadtree.h"
#include <algorithm>

Quadtree::Quadtree(BoundingSquare bounds, int capacity)
    : bounds(bounds), capacity(capacity), NW(nullptr), NE(nullptr), SW(nullptr), SE(nullptr) {}
void Quadtree::Insert(const Point& point) {
    // If the point is outside the bounds of the current quadtree node, do nothing
    if (!bounds.Contains(point)) {
        return;
    }

    // If the current node has space, just add the point
    if (points.size() < capacity) {
        points.push_back(point);
    }
    else {
        // Otherwise, subdivide the node if it hasn't been subdivided already
        if (NW == nullptr) {
            Subdivide();
        }

        // Insert the point into the correct child quadrant
        if (NW->bounds.Contains(point)) {
            NW->Insert(point);
        }
        else if (NE->bounds.Contains(point)) {
            NE->Insert(point);
        }
        else if (SW->bounds.Contains(point)) {
            SW->Insert(point);
        }
        else if (SE->bounds.Contains(point)) {
            SE->Insert(point);
        }
    }
}

void Quadtree::Subdivide() {
    // Split the bounding box into four quadrants and create the child nodes
    double midX = (bounds.minX + bounds.maxX) / 2;
    double midY = (bounds.minY + bounds.maxY) / 2;

    NW = new Quadtree(BoundingSquare(bounds.minX, bounds.minY, midX, midY), capacity);
    NE = new Quadtree(BoundingSquare(midX, bounds.minY, bounds.maxX, midY), capacity);
    SW = new Quadtree(BoundingSquare(bounds.minX, midY, midX, bounds.maxY), capacity);
    SE = new Quadtree(BoundingSquare(midX, midY, bounds.maxX, bounds.maxY), capacity);
}

bool Quadtree::IsSubdivided() const {
    return NW != nullptr;
}

std::vector<Point> Quadtree::QueryRange(const BoundingSquare& range) const {
    std::vector<Point> result;

    // If the range does not overlap with this quadtree node's bounds, return empty
    if (!bounds.Intersects(range)) {
        return result;
    }
    auto p = Point(2, 3);
    // Add points in this quadtree node's region to the result
    for (const auto& point : points) {
        if (range.Contains(point)) {
            result.push_back(point);
        }
    }

    // If the node is subdivided, query the child quadrants
    if (IsSubdivided()) {
        // Collect points from all child quadrants
        auto nw_points = NW->QueryRange(range);
        result.insert(result.end(), nw_points.begin(), nw_points.end());

        auto ne_points = NE->QueryRange(range);
        result.insert(result.end(), ne_points.begin(), ne_points.end());

        auto sw_points = SW->QueryRange(range);
        result.insert(result.end(), sw_points.begin(), sw_points.end());

        auto se_points = SE->QueryRange(range);
        result.insert(result.end(), se_points.begin(), se_points.end());
    }

    return result;
}

const std::vector<Point>& Quadtree::GetPoints() const {
    return points;
}
