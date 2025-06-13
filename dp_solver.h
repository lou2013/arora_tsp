#pragma once

#include "geometry.h"
#include "dpstate.h"
#include "portal_generator.h"

class DPSolver {
public:
    DPSolver(const std::vector<Point>& points, int portals_per_edge);
    double Solve();
    double Solve2();

    // --- FIX IS HERE ---
    // Moved ComputeMSTCost to public so our helper function can access it.
    double ComputeMSTCost(const std::vector<Point>& points) const;


private:
    // Core state
    std::vector<Point> inputPoints;
    int k; // Portals per edge
    PortalGenerator portalGen;
    DPCache memo; // Memoization table

    // Main recursive functions
    DPCostTable SolveRecursive(const BoundingSquare& square, int depth);
    DPCostTable ComputeBaseCase(const BoundingSquare& square, const std::vector<Point>& local_points, int depth);
    DPCostTable Combine(
        const DPCostTable& nw_costs, const DPCostTable& ne_costs,
        const DPCostTable& sw_costs, const DPCostTable& se_costs,
        const BoundingSquare& parent_square, int depth);
};