#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <string>
#include "geometry.h"

bool ReadPointsFromFile(const std::string& filename, std::vector<Point>& points);
void ShufflePoints(std::vector<Point>& points);

#endif // UTILITIES_H
