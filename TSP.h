#ifndef TSP_H
#define TSP_H

#include "geometry.h"
#include <string>
#include <vector>

class TSP {
public:
    TSP(const std::string& filename, int iterations = 1);
    void Run();

private:
    std::vector<Point> points;
    int numIterations;
};

#endif // TSP_H
