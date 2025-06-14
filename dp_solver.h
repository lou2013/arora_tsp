#pragma once

#include "geometry.h"
#include "dpstate.h" // Contains the new PathInfo and Tour structs
#include "portal_generator.h"

class DPSolver {
public:
    DPSolver(const std::vector<Point>& points, int portals_per_edge);

    // ================== MODIFICATION START ==================
    // Changed return types from double to the new Tour struct
    // which contains both the cost and the final path.
    Tour Solve();
    Tour Solve2();
    // =================== MODIFICATION END ===================


    // --- FIX IS HERE ---
    // Moved ComputeMSTCost to public so our helper function can access it.
    // MODIFIED: This now returns PathInfo (cost + edges).
    PathInfo ComputeMSTCost(const std::vector<Point>& points) const;


private:
    // Core state
    std::vector<Point> inputPoints;
    int k; // Portals per edge
    PortalGenerator portalGen;
    DPCache memo; // Memoization table

    // ================== MODIFICATION START ==================
    // Main recursive functions now return DPCostTable with PathInfo
    DPCostTable SolveRecursive(const BoundingSquare& square, int depth);
    DPCostTable ComputeBaseCase(const BoundingSquare& square, const std::vector<Point>& local_points, int depth);
    DPCostTable Combine(
        const DPCostTable& nw_costs, const DPCostTable& ne_costs,
        const DPCostTable& sw_costs, const DPCostTable& se_costs,
        const BoundingSquare& parent_square, int depth);
    // =================== MODIFICATION END ===================
};