#include "utilities.h"
#include <fstream>
#include <sstream>
#include <random>
#include <algorithm>

bool ReadPointsFromFile(const std::string& filename, std::vector<Point>& points) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    double x, y;
    while (file >> x >> y) {
        points.emplace_back(x, y);
    }
    return true;
}

void ShufflePoints(std::vector<Point>& points) {
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(points.begin(), points.end(), g);
}
