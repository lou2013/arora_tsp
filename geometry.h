#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <algorithm>
#include <limits>
#include <cmath>
#include <vector>
#include <tuple>
struct Point {
    double x, y;

    Point(double x = 0, double y = 0) : x(x), y(y) {}

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }

};

inline double Distance(const Point& a, const Point& b) {
    return std::sqrt(std::pow(a.x - b.x, 2) + std::pow(a.y - b.y, 2));
}

class BoundingSquare {
public:
    double minX, minY, maxX, maxY;
    BoundingSquare(){}
    BoundingSquare(double minX, double minY, double maxX, double maxY)
        : minX(minX), minY(minY), maxX(maxX), maxY(maxY) {}

    // Check if a point is within this bounding box
    bool Contains(const Point& point) const {
        return point.x >= minX && point.x <= maxX && point.y >= minY && point.y <= maxY;
    }
    // Check if this bounding box intersects with another
    bool Intersects(const BoundingSquare& other) const {
        return !(other.maxX < minX || other.minX > maxX || other.maxY < minY || other.minY > maxY);
    }
    double SideLength() const {
        return std::max(maxX - minX, maxY - minY);
    }
    std::vector<Point> FilterPoints(const std::vector<Point>& points) const {
        std::vector<Point> filtered;
        for (const auto& p : points) {
            if (p.x >= minX && p.x <= maxX && p.y >= minY && p.y <= maxY) {
                filtered.push_back(p);
            }
        }
        return filtered;
    }

    static BoundingSquare ComputeBoundingSquare(const std::vector<Point>& points) {
        if (points.empty()) {
            return {0, 0, 0, 0};
        }
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double max_y = std::numeric_limits<double>::lowest();

        for (const auto& p : points) {
            min_x = std::min(min_x, p.x);
            max_x = std::max(max_x, p.x);
            min_y = std::min(min_y, p.y);
            max_y = std::max(max_y, p.y);
        }
        return {min_x, min_y, max_x, max_y};
    }

    BoundingSquare ExpandToPowerOfTwo() const {
        double side = std::max(maxX - minX, maxY - minY);
        double powerOfTwoSide = std::pow(2, std::ceil(std::log2(side)));
        return {minX, minY, minX + powerOfTwoSide, minY + powerOfTwoSide};
    }
    double Size() const;
    bool operator<(const BoundingSquare& other) const {
        return std::tie(minX, minY, maxX, maxY) < std::tie(other.minX, other.minY, other.maxX, other.maxY);
    }
};

#endif // GEOMETRY_H
