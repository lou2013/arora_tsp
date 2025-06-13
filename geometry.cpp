#include "geometry.h"

//BoundingSquare BoundingSquare::ComputeBoundingSquare(const std::vector<Point>& points) {
    //double minX = points[0].x;
    //double maxX = points[0].x;
    //double minY = points[0].y;
    //double maxY = points[0].y;
 
    //for (const auto& p : points) {
        //minX = std::min(minX, p.x);
        //maxX = std::max(maxX, p.x);
        //minY = std::min(minY, p.y);
        //maxY = std::max(maxY, p.y);
    //}

    //return BoundingSquare(minX, minY, maxX, maxY);
//}

// BoundingSquare BoundingSquare::ExpandToPowerOfTwo() const {
//     double width = maxX - minX;
//     double height = maxY - minY;
//     double maxDim = std::max(width, height);
//     int powerOfTwo = 1;
//     while (powerOfTwo < maxDim) {
//         powerOfTwo *= 2;
//     }
//     double centerX = (minX + maxX) / 2.0;
//     double centerY = (minY + maxY) / 2.0;

//     double halfSide = powerOfTwo / 2.0;

//     return BoundingSquare(
//         centerX - halfSide,
//         centerY - halfSide,
//         centerX + halfSide,
//         centerY + halfSide
//     );
// }

double BoundingSquare::Size() const {
    return std::max(maxX - minX, maxY - minY);
}

// std::vector<Point> BoundingSquare::FilterPoints(const std::vector<Point>& points) const {
//     std::vector<Point> result;
//     for (const auto& p : points) {
//         if (Contains(p)) {
//             result.push_back(p);
//         }
//     }
//     return result;
// }
