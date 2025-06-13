#ifndef FAST_ARORA_H
#define FAST_ARORA_H

#include <vector>
#include "geometry.h"

namespace FastArora {
    // Runs the algorithm and returns the length of the best tour found
    double Run(std::vector<Point>& points);
}

#endif // FAST_ARORA_H
