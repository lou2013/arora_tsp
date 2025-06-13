#include "dp_solver.h"
#include <iostream>
#include <numeric>
#include <queue>
#include <map>
#include <vector>
#include <algorithm>

#pragma region Helper Functions
void generate_pairings_recursive(
    int portal_idx,
    std::vector<int>& open_portals,
    PortalPairing& current_pairing,
    std::vector<PortalPairing>& all_pairings,
    int num_portals)
{
    if (portal_idx == num_portals) {
        if (open_portals.empty()) {
            for(auto& p : current_pairing.pairs) {
                if (p.first > p.second) std::swap(p.first, p.second);
            }
            std::sort(current_pairing.pairs.begin(), current_pairing.pairs.end());
            all_pairings.push_back(current_pairing);
        }
        return;
    }
    open_portals.push_back(portal_idx);
    generate_pairings_recursive(portal_idx + 1, open_portals, current_pairing, all_pairings, num_portals);
    open_portals.pop_back();
    if (!open_portals.empty()) {
        int last_open = open_portals.back();
        open_portals.pop_back();
        current_pairing.pairs.push_back({last_open, portal_idx});
        generate_pairings_recursive(portal_idx + 1, open_portals, current_pairing, all_pairings, num_portals);
        current_pairing.pairs.pop_back();
        open_portals.push_back(last_open);
    }
}
// A map to memoize results for performance
std::map<int, std::vector<PortalPairing>> memo_all_pairings;

// Generates ALL non-crossing pairings, including subsets where some portals are unpaired.
std::vector<PortalPairing> GenerateAllPossibleNonCrossingPairings(int num_portals) {
    if (num_portals < 0) return {};
    if (num_portals == 0) return { {} }; // Base case: one option, the empty pairing
    if (memo_all_pairings.count(num_portals)) {
        return memo_all_pairings[num_portals];
    }

    std::vector<PortalPairing> result;

    // Case 1: Portal (num_portals - 1) is left unpaired.
    // Get all pairings of the first (n-1) portals.
    auto pairings_without_last = GenerateAllPossibleNonCrossingPairings(num_portals - 1);
    result.insert(result.end(), pairings_without_last.begin(), pairings_without_last.end());

    // Case 2: Portal (num_portals - 1) is paired with some portal j.
    // j must be such that the number of points between them is even.
    for (int j = 0; j < num_portals - 1; ++j) {
        int points_between = (num_portals - 1) - j - 1;
        if (points_between % 2 != 0) continue; // Must be even

        // Recursively find pairings for the two independent subproblems
        auto inner_pairings = GenerateAllPossibleNonCrossingPairings(points_between);
        auto outer_pairings = GenerateAllPossibleNonCrossingPairings(j);

        // Combine the results
        for (const auto& outer_p : outer_pairings) {
            for (const auto& inner_p : inner_pairings) {
                PortalPairing new_pairing = outer_p; // Start with the outer

                // Add the new pair {j, num_portals - 1}
                new_pairing.pairs.push_back({ j, num_portals - 1 });

                // Add the inner pairings, shifting their indices
                for (auto inner_pair : inner_p.pairs) {
                    new_pairing.pairs.push_back({ inner_pair.first + j + 1, inner_pair.second + j + 1 });
                }

                // Sort to maintain a canonical representation for the pairing
                for (auto& p : new_pairing.pairs) { if (p.first > p.second) std::swap(p.first, p.second); }
                std::sort(new_pairing.pairs.begin(), new_pairing.pairs.end());
                result.push_back(new_pairing);
            }
        }
    }

    // Remove duplicates and memoize
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return memo_all_pairings[num_portals] = result;
}
#pragma endregion

DPSolver::DPSolver(const std::vector<Point>& points, int portals_per_edge)
    : inputPoints(points), k(portals_per_edge) {}

double DPSolver::Solve() {
    if (inputPoints.empty()) return 0.0;
    std::cout << "--- Starting Solver ---" << std::endl;
    BoundingSquare root_square = BoundingSquare::ComputeBoundingSquare(inputPoints);
    root_square = root_square.ExpandToPowerOfTwo();

    DPCostTable final_costs = SolveRecursive(root_square, 0);

    std::cout << "--- Solver Finished ---" << std::endl;


    double min_path_cost = std::numeric_limits<double>::infinity();
    
    for (const auto& map_entry : final_costs) {
        const PortalPairing& pairing = map_entry.first;
        double cost = map_entry.second;

        // We are looking for the cheapest path that enters and exits the bounding box once.
        if (pairing.pairs.size() == 1) {
            if (cost < min_path_cost) {
                min_path_cost = cost;
            }
        }
    }

    return min_path_cost;
    // --- MODIFICATION END ---
}

 double DPSolver::Solve2() {
     if (inputPoints.empty()) return 0.0;
     std::cout << "--- Starting Solver ---" << std::endl;
     BoundingSquare root_square = BoundingSquare::ComputeBoundingSquare(inputPoints);
     root_square = root_square.ExpandToPowerOfTwo();
     DPCostTable final_costs = SolveRecursive(root_square, 0);
     PortalPairing empty_pairing;
     std::cout << "--- Solver Finished ---" << std::endl;
     if (final_costs.count(empty_pairing)) {
         return final_costs.at(empty_pairing);
     }
     return std::numeric_limits<double>::infinity();
 }

DPCostTable DPSolver::SolveRecursive(const BoundingSquare& square, int depth) {
    std::cout << std::string(depth * 2, ' ') << "Solving for square [" << square.minX << ", " << square.minY << "] at depth " << depth << std::endl;
    if (memo.count(square)) { return memo.at(square); }

    std::vector<Point> localPoints = square.FilterPoints(inputPoints);
    if (localPoints.size() <= (size_t)4 * k || depth > 5) {
        return memo[square] = ComputeBaseCase(square, localPoints, depth);
    }

    double midX = (square.minX + square.maxX) / 2;
    double midY = (square.minY + square.maxY) / 2;

    BoundingSquare nw{square.minX, midY, midX, square.maxY};
    BoundingSquare ne{midX, midY, square.maxX, square.maxY};
    BoundingSquare sw{square.minX, square.minY, midX, midY};
    // FIX: Corrected the Y coordinates for the South-East square
    BoundingSquare se{midX, square.minY, square.maxX, midY};

    DPCostTable nw_costs = SolveRecursive(nw, depth + 1);
    DPCostTable ne_costs = SolveRecursive(ne, depth + 1);
    DPCostTable sw_costs = SolveRecursive(sw, depth + 1);
    DPCostTable se_costs = SolveRecursive(se, depth + 1);

    return memo[square] = Combine(nw_costs, ne_costs, sw_costs, se_costs, square, depth);
}

DPCostTable DPSolver::ComputeBaseCase(const BoundingSquare& square, const std::vector<Point>& local_points, int depth) {
    DPCostTable cost_table;
    auto portals = portalGen.GeneratePortals(square, k);
    int m = 4 * k;
    std::vector<PortalPairing> all_possible_pairings;
    auto full_pairings = GenerateAllPossibleNonCrossingPairings(m);
    all_possible_pairings.insert(all_possible_pairings.end(), full_pairings.begin(), full_pairings.end());
    if (m > 1) {
        for (int i = 0; i < m; ++i) {
            for (int j = i + 1; j < m; ++j) {
                all_possible_pairings.push_back({{{i, j}}});
            }
        }
    }
    all_possible_pairings.push_back({});

    for(const auto& pairing : all_possible_pairings) {
        std::vector<Point> mst_points = local_points;
        for (const auto& p : pairing.pairs) {
            if (p.first < m && p.second < m) {
                mst_points.push_back(portals[p.first]);
                mst_points.push_back(portals[p.second]);
            }
        }
        cost_table[pairing] = ComputeMSTCost(mst_points);
    }
    
    std::cout << std::string(depth * 2, ' ') << "--> Base Case! Points: " << local_points.size() << ". Cost for empty pairing: " << cost_table[PortalPairing{}] << std::endl;
    return cost_table;
}

PortalPairing get_required_pairing(int child_idx, const std::map<int, int>& connections) {
    PortalPairing req;
    std::vector<int> unpaired;
    for (int i = 0; i < 4; ++i) {
        int global_portal_id = child_idx * 4 + i;
        if (connections.find(global_portal_id) == connections.end()) {
            unpaired.push_back(i);
        } else {
            int partner_global_id = connections.at(global_portal_id);
            if (partner_global_id / 4 == child_idx) {
                if (global_portal_id < partner_global_id) {
                    req.pairs.push_back({i, partner_global_id % 4});
                }
            }
        }
    }

    if (unpaired.size() == 2) {
        req.pairs.push_back({unpaired[0], unpaired[1]});
    } else if (unpaired.size() == 4) {
        req.pairs.push_back({0, 3});
        req.pairs.push_back({1, 2});
    }

    for(auto& p : req.pairs) { if(p.first > p.second) std::swap(p.first, p.second); }
    std::sort(req.pairs.begin(), req.pairs.end());
    return req;
}

// Helper struct for the edges in our portal graph
struct PortalEdge {
    int u, v;
    double weight;

    bool operator<(const PortalEdge& other) const {
        return weight < other.weight;
    }
};

// Helper struct for Kruskal's Algorithm (to find MST)
struct DSU {
    std::vector<int> parent;
    DSU(int n) {
        parent.resize(n);
        std::iota(parent.begin(), parent.end(), 0);
    }
    int find(int i) {
        if (parent[i] == i) return i;
        return parent[i] = find(parent[i]);
    }
    void unite(int i, int j) {
        int root_i = find(i);
        int root_j = find(j);
        if (root_i != root_j) {
            parent[root_i] = root_j;
        }
    }
};

DPCostTable DPSolver::Combine(
    const DPCostTable& nw_costs, const DPCostTable& ne_costs,
    const DPCostTable& sw_costs, const DPCostTable& se_costs,
    const BoundingSquare& parent_square, int depth)
{
    std::cout << std::string(depth * 2, ' ') << "Combining for square [" << parent_square.minX << ", " << parent_square.minY << "]" << std::endl;
    DPCostTable parent_cost_table;
    if (k != 1) return parent_cost_table;

    // Step 1: Get portal locations and base costs for each child
    double midX = (parent_square.minX + parent_square.maxX) / 2;
    double midY = (parent_square.minY + parent_square.maxY) / 2;
    BoundingSquare nw_sq{parent_square.minX, midY, midX, parent_square.maxY},     ne_sq{midX, midY, parent_square.maxX, parent_square.maxY};
    BoundingSquare sw_sq{parent_square.minX, parent_square.minY, midX, midY}, se_sq{midX, parent_square.minY, parent_square.maxX, midY};
    
    auto nw_p = portalGen.GeneratePortals(nw_sq, k), ne_p = portalGen.GeneratePortals(ne_sq, k);
    auto sw_p = portalGen.GeneratePortals(sw_sq, k), se_p = portalGen.GeneratePortals(se_sq, k);
    
    const std::vector<const DPCostTable*> costs = {&nw_costs, &ne_costs, &sw_costs, &se_costs};
    const std::vector<std::vector<Point>*> portals = {&nw_p, &ne_p, &sw_p, &se_p};

    auto get_cost = [&](int child_idx, const PortalPairing& p) {
        return costs[child_idx]->count(p) ? costs[child_idx]->at(p) : std::numeric_limits<double>::infinity();
    };
    double nw_base_cost = get_cost(0, {});
    double ne_base_cost = get_cost(1, {});
    double sw_base_cost = get_cost(2, {});
    double se_base_cost = get_cost(3, {});
    double children_base_total_cost = nw_base_cost + ne_base_cost + sw_base_cost + se_base_cost;


    // Step 2: Define the graph for the connection problem.
    // We have 12 nodes: 4 parent portals (0-3) and 8 center-facing child portals (4-11).
    // Node mapping:
    // 0-3: Parent portals {Top, Right, Bottom, Left}
    // 4: NW-right, 5: NW-bottom
    // 6: NE-bottom, 7: NE-left
    // 8: SW-top,   9: SW-right
    // 10: SE-left, 11: SE-top
    std::vector<PortalEdge> edges;

    // Add "stitch" edges between adjacent child portals
    edges.push_back({4, 7, Distance(nw_p[1], ne_p[3])}); // NW-right to NE-left
    edges.push_back({5, 8, Distance(nw_p[2], sw_p[0])}); // NW-bottom to SW-top
    edges.push_back({6, 11, Distance(ne_p[2], se_p[0])});// NE-bottom to SE-top
    edges.push_back({9, 10, Distance(sw_p[1], se_p[3])});// SW-right to SE-left

    // Add edges for paths *through* the children.
    // The cost of a path is the EXTRA cost over the child's base cost.
    auto add_child_path_edges = [&](int child_idx, int p1_local, int p2_local, int n1_global, int n2_global){
        double path_cost = get_cost(child_idx, {{{p1_local, p2_local}}});
        double base_cost = get_cost(child_idx, {});
        if (path_cost != std::numeric_limits<double>::infinity()) {
            edges.push_back({n1_global, n2_global, path_cost - base_cost});
        }
    };

    // Paths connecting parent portals to center-facing portals
    add_child_path_edges(0, 0, 1, 0, 4); add_child_path_edges(0, 0, 2, 0, 5); add_child_path_edges(0, 3, 1, 3, 4); add_child_path_edges(0, 3, 2, 3, 5); // NW
    add_child_path_edges(1, 0, 2, 0, 6); add_child_path_edges(1, 0, 3, 0, 7); add_child_path_edges(1, 1, 2, 1, 6); add_child_path_edges(1, 1, 3, 1, 7); // NE
    add_child_path_edges(2, 0, 1, 2, 9); add_child_path_edges(2, 2, 0, 2, 8); add_child_path_edges(2, 2, 1, 2, 9); add_child_path_edges(2, 3, 0, 3, 8); add_child_path_edges(2, 3, 1, 3, 9); // SW (Note: Parent portal 2 is bottom, 3 is left)
    add_child_path_edges(3, 0, 3, 2, 10);add_child_path_edges(3, 0, 1, 2, 11);add_child_path_edges(3, 1, 3, 1, 10);add_child_path_edges(3, 1, 0, 1, 11); // SE
    
    // Paths between center-facing portals of the same child
    add_child_path_edges(0, 1, 2, 4, 5);   // NW
    add_child_path_edges(1, 2, 3, 6, 7);   // NE
    add_child_path_edges(2, 0, 1, 8, 9);   // SW
    add_child_path_edges(3, 0, 3, 11, 10); // SE -- check portal indices, should be 11(top) and 10(left) vs local 0 and 3
    
    std::sort(edges.begin(), edges.end());

    // Step 3: Iterate through all possible parent pairings and solve the MST problem
    auto all_parent_pairings = GenerateAllPossibleNonCrossingPairings(4);

    for (const auto& p_pairing : all_parent_pairings) {
        DSU dsu(12);
        double connection_mst_cost = 0;
        int edges_count = 0;

        // Force connections for the current parent pairing by pre-uniting their nodes
        for (const auto& pair : p_pairing.pairs) {
            dsu.unite(pair.first, pair.second);
        }

        // Run Kruskal's algorithm to find the MST cost for the remaining connections
        for (const auto& edge : edges) {
            if (dsu.find(edge.u) != dsu.find(edge.v)) {
                dsu.unite(edge.u, edge.v);
                connection_mst_cost += edge.weight;
                edges_count++;
            }
        }
        
        parent_cost_table[p_pairing] = children_base_total_cost + connection_mst_cost;
    }
    
    std::cout << std::string(depth * 2, ' ') << "--> Combined. Populated DP table. Cost for empty pairing: " << parent_cost_table[PortalPairing{}] << std::endl;
    return parent_cost_table;
}
double DPSolver::ComputeMSTCost(const std::vector<Point>& points) const {
    size_t n = points.size();
    if (n <= 1) return 0.0;
    std::vector<double> min_cost(n, std::numeric_limits<double>::infinity());
    std::vector<bool> visited(n, false);
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
    min_cost[0] = 0.0;
    pq.push({0.0, 0});
    double total_cost = 0.0;
    int edges = 0;
    while (!pq.empty() && edges < n) {
        double cost = pq.top().first;
        int u = pq.top().second;
        pq.pop();
        if (visited[u]) continue;
        visited[u] = true;
        total_cost += cost;
        edges++;
        for (int v = 0; v < n; ++v) {
            if (!visited[v]) {
                double weight = Distance(points[u], points[v]);
                if (weight < min_cost[v]) {
                    min_cost[v] = weight;
                    pq.push({weight, v});
                }
            }
        }
    }
    return total_cost;
}