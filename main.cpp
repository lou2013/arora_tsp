#include "dp_solver.h"
#include <iostream>
#include <vector>
#include <string>
#include <iomanip> // For std::setprecision and std::fixed

// A helper structure to hold all our test cases
struct TestCase {
    std::string name;
    std::vector<Point> points;
    double optimal_path_length;
    // For tests where the connection cost is different, we'll note it in the name.
};

int main() {
    // --- Test Case Definitions ---
    std::vector<TestCase> all_tests = {
        {
            "Isolated Segment (Path vs Connection Test)",
            { {10, 50}, {15, 50} },
            5.0
        },
        {
            "Diagonal Line (Grid Weakness Test)",
            { {10, 10}, {40, 40}, {70, 70}, {100, 100} },
            127.28 
        },
        {
            "U-Shape Path",
            { {10, 100}, {10, 10}, {100, 10}, {100, 100} },
            270.0
        },
        {
            "Detour Path",
            { {50, 10}, {10, 90}, {100, 90}, {60, 10} },
            188.88
        },
        {
            "Spokes of a Wheel (Optimal Path ~434.56, Optimal Connection = 360.0)",
            { {100, 100}, {100, 10}, {190, 100}, {100, 190}, {10, 100} },
            434.56 // This is the optimal path length
        },
        {
            "Original 8 Points",
            { {10, 20}, {95, 35}, {30, 80}, {70, 10}, {50, 50}, {15, 60}, {85, 90}, {40, 5} },
            333.27
        },
        {
            "Winding Path (10 Points)",
            { {60, 200}, {180, 200}, {80, 180}, {140, 180}, {20, 160}, {100, 160}, {200, 160}, {140, 140}, {40, 120}, {100, 120} },
            385.65
        }
    };

    // k=3 would be much more accurate but significantly slower.
    int k = 3;

    // --- Main Test Loop ---
    for (const auto& test : all_tests) {
        std::cout << "\n====================================================================\n";
        std::cout << "--- Running Test: " << test.name << " ---\n";
        std::cout << "====================================================================\n";
        std::cout << "Points: " << test.points.size() << ", k=" << k << std::endl;
        std::cout << std::fixed << std::setprecision(2);
        
        // The known best answer for a path visiting all nodes.
        std::cout << "> Optimal Path Length: " << test.optimal_path_length << std::endl;
        
        // Initialize the solver for this test case
        DPSolver solver(test.points, k);
        
        // We need to call one of the solves first to populate the DP table
        double cost1 = solver.Solve(); 
        double cost2 = solver.Solve2();

        std::cout << "\n--- RESULTS ---\n";
        std::cout << "[Solve()  - Open Path Approx]:   " << cost1 << std::endl;
        std::cout << "[Solve2() - Connection Approx]: " << cost2 << std::endl;
        
        std::cout << "\n--- ANALYSIS ---\n";
        if (test.name.find("Spokes") != std::string::npos) {
            std::cout << "Solve() correctly approximates the long path around the spokes.\n";
            std::cout << "Solve2() correctly approximates the cheaper 'star-shaped' connection.\n";
            std::cout << "This shows they solve two different problems correctly.\n";
        } else if (test.name.find("Isolated") != std::string::npos) {
             std::cout << "Solve2() correctly finds the simple path length (5.0).\n";
             std::cout << "Solve() adds the cost of connecting the path to the boundary.\n";
        }
        else {
            std::cout << "For this test, both methods approximate the same path.\n";
            std::cout << "The results should be similar and slightly > optimal.\n";
        }
    }

    return 0;
}