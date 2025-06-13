//#include "TSP.h"
////#include "fast_arora.h"
//#include "utilities.h"
//#include <iostream>
//
//TSP::TSP(const std::string& filename, int iterations)
//    : numIterations(iterations) {
//    // Read points from file
//    if (!ReadPointsFromFile(filename, points)) {
//        std::cerr << "Failed to read points from file: " << filename << std::endl;
//        exit(1);
//    }
//}
//
//void TSP::Run() {
//    std::cout << "Running Arora's Algorithm on " << points.size() << " points..." << std::endl;
//
//    double bestLength = 1e9;
//    std::vector<Point> bestTour;
//
//    for (int i = 0; i < numIterations; ++i) {
//        std::vector<Point> tour = points;
//
//        // Shuffle points to get different starting conditions (optional)
//        ShufflePoints(tour);
//
//        double length = FastArora::Run(tour);
//
//        if (length < bestLength) {
//            bestLength = length;
//            bestTour = tour;
//        }
//
//        std::cout << "Iteration " << i + 1 << ": Tour length = " << length << std::endl;
//    }
//
//    std::cout << "Best tour length: " << bestLength << std::endl;
//}
