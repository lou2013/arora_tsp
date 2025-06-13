//#include "fast_arora.h"
//#include "geometry.h"
//#include "perturbation.h"
//#include "quad_tree.h"
//#include "portalization.h"
//#include "dynamic_programming.h"
//
//#include <iostream>
//
//namespace FastArora {
//
//    double Run(std::vector<Point>& points) {
//        // Step 1: Perturb the points (snap to grid)
//        Perturbation perturbation(points);
//        perturbation.Apply();
//
//        // Step 2: Build quadtree structure
//        QuadTree qt(perturbation.GetBoundingSquare(), perturbation.GetPoints());
//
//        // Step 3: Place portals
//        Portalization portalization(qt);
//        portalization.PlacePortals();
//
//        // Step 4: Run dynamic programming to find minimal tour
//        DynamicProgramming dp(qt, portalization);
//        std::vector<Point> tour = dp.Solve();
//
//        // Compute length of final tour
//        double length = 0.0;
//        for (size_t i = 0; i < tour.size(); ++i) {
//            length += Distance(tour[i], tour[(i + 1) % tour.size()]);
//        }
//
//        points = tour; // overwrite with computed tour
//        return length;
//    }
//
//}
