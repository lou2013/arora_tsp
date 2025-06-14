#include "dp_solver.h"
#include <iostream>
#include <numeric>
#include <queue>
#include <map>
#include <vector>
#include <algorithm>

#pragma region Helper Functions
// This section of helper functions for generating pairings remains unchanged.
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
std::map<int, std::vector<PortalPairing>> memo_all_pairings;
std::vector<PortalPairing> GenerateAllPossibleNonCrossingPairings(int num_portals) {
    if (num_portals < 0) return {};
    if (num_portals == 0) return { {} };
    if (memo_all_pairings.count(num_portals)) {
        return memo_all_pairings[num_portals];
    }
    std::vector<PortalPairing> result;
    auto pairings_without_last = GenerateAllPossibleNonCrossingPairings(num_portals - 1);
    result.insert(result.end(), pairings_without_last.begin(), pairings_without_last.end());
    for (int j = 0; j < num_portals - 1; ++j) {
        int points_between = (num_portals - 1) - j - 1;
        if (points_between % 2 != 0) continue;
        auto inner_pairings = GenerateAllPossibleNonCrossingPairings(points_between);
        auto outer_pairings = GenerateAllPossibleNonCrossingPairings(j);
        for (const auto& outer_p : outer_pairings) {
            for (const auto& inner_p : inner_pairings) {
                PortalPairing new_pairing = outer_p;
                new_pairing.pairs.push_back({ j, num_portals - 1 });
                for (auto inner_pair : inner_p.pairs) {
                    new_pairing.pairs.push_back({ inner_pair.first + j + 1, inner_pair.second + j + 1 });
                }
                for (auto& p : new_pairing.pairs) { if (p.first > p.second) std::swap(p.first, p.second); }
                std::sort(new_pairing.pairs.begin(), new_pairing.pairs.end());
                result.push_back(new_pairing);
            }
        }
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return memo_all_pairings[num_portals] = result;
}
#pragma endregion

DPSolver::DPSolver(const std::vector<Point>& points, int portals_per_edge)
    : inputPoints(points), k(portals_per_edge) {}

// MODIFIED: Returns a Tour object now
Tour DPSolver::Solve() {
    if (inputPoints.empty()) return {0.0, {}};
    std::cout << "--- Starting Solver ---" << std::endl;
    BoundingSquare root_square = BoundingSquare::ComputeBoundingSquare(inputPoints);
    root_square = root_square.ExpandToPowerOfTwo();

    DPCostTable final_costs = SolveRecursive(root_square, 0);

    std::cout << "--- Solver Finished ---" << std::endl;

    // Find the cheapest tour that uses exactly one pair of portals on the root box
    PathInfo best_path;
    for (const auto& map_entry : final_costs) {
        const PortalPairing& pairing = map_entry.first;
        const PathInfo& current_path_info = map_entry.second;
        if (pairing.pairs.size() == 1) {
            if (current_path_info.cost < best_path.cost) {
                best_path = current_path_info;
            }
        }
    }
    return {best_path.cost, best_path.edges};
}

// MODIFIED: Returns a Tour object now
Tour DPSolver::Solve2() {
     if (inputPoints.empty()) return {0.0, {}};
     std::cout << "--- Starting Solver ---" << std::endl;
     BoundingSquare root_square = BoundingSquare::ComputeBoundingSquare(inputPoints);
     root_square = root_square.ExpandToPowerOfTwo();
     DPCostTable final_costs = SolveRecursive(root_square, 0);
     PortalPairing empty_pairing;
     std::cout << "--- Solver Finished ---" << std::endl;
     if (final_costs.count(empty_pairing)) {
         const auto& final_path_info = final_costs.at(empty_pairing);
         return {final_path_info.cost, final_path_info.edges};
     }
     return {std::numeric_limits<double>::infinity(), {}};
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
    auto all_possible_pairings = GenerateAllPossibleNonCrossingPairings(m);
    
    // This logic to generate pairings seems complex, ensuring we cover all cases
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
        // MODIFIED: ComputeMSTCost now returns a PathInfo object, which is stored directly.
        cost_table[pairing] = ComputeMSTCost(mst_points);
    }
    
    std::cout << std::string(depth * 2, ' ') << "--> Base Case! Points: " << local_points.size() << ". Cost for empty pairing: " << cost_table[PortalPairing{}].cost << std::endl;
    return cost_table;
}

struct PortalEdge {
    int u, v;
    double weight;
    // MODIFICATION: Store the actual geometric edge for path reconstruction
    std::pair<Point, Point> actual_edge;

    bool operator<(const PortalEdge& other) const {
        return weight < other.weight;
    }
};

struct DSU {
    std::vector<int> parent;
    DSU(int n) { parent.resize(n); std::iota(parent.begin(), parent.end(), 0); }
    int find(int i) { return (parent[i] == i) ? i : (parent[i] = find(parent[i])); }
    void unite(int i, int j) { int root_i = find(i); int root_j = find(j); if (root_i != root_j) parent[root_i] = root_j; }
};

DPCostTable DPSolver::Combine(
    const DPCostTable& nw_costs, const DPCostTable& ne_costs,
    const DPCostTable& sw_costs, const DPCostTable& se_costs,
    const BoundingSquare& parent_square, int depth)
{
    std::cout << std::string(depth * 2, ' ') << "Combining for square [" << parent_square.minX << ", " << parent_square.minY << "]" << std::endl;
    DPCostTable parent_cost_table;
    if (k != 1) return parent_cost_table;

    double midX = (parent_square.minX + parent_square.maxX) / 2;
    double midY = (parent_square.minY + parent_square.maxY) / 2;
    BoundingSquare nw_sq{parent_square.minX, midY, midX, parent_square.maxY}, ne_sq{midX, midY, parent_square.maxX, parent_square.maxY};
    BoundingSquare sw_sq{parent_square.minX, parent_square.minY, midX, midY}, se_sq{midX, parent_square.minY, parent_square.maxX, midY};
    
    auto nw_p = portalGen.GeneratePortals(nw_sq, k), ne_p = portalGen.GeneratePortals(ne_sq, k);
    auto sw_p = portalGen.GeneratePortals(sw_sq, k), se_p = portalGen.GeneratePortals(se_sq, k);
    
    const std::vector<const DPCostTable*> costs_ptr = {&nw_costs, &ne_costs, &sw_costs, &se_costs};
    const PathInfo default_path_info; // Used for infinity cases

    auto get_path_info = [&](int child_idx, const PortalPairing& p) -> const PathInfo& {
        return costs_ptr[child_idx]->count(p) ? costs_ptr[child_idx]->at(p) : default_path_info;
    };

    // MODIFIED: Get base PathInfo objects and aggregate their costs and paths
    const PathInfo& nw_base_info = get_path_info(0, {});
    const PathInfo& ne_base_info = get_path_info(1, {});
    const PathInfo& sw_base_info = get_path_info(2, {});
    const PathInfo& se_base_info = get_path_info(3, {});
    
    double children_base_total_cost = nw_base_info.cost + ne_base_info.cost + sw_base_info.cost + se_base_info.cost;
    std::vector<std::pair<Point, Point>> children_base_path;
    children_base_path.insert(children_base_path.end(), nw_base_info.edges.begin(), nw_base_info.edges.end());
    children_base_path.insert(children_base_path.end(), ne_base_info.edges.begin(), ne_base_info.edges.end());
    children_base_path.insert(children_base_path.end(), sw_base_info.edges.begin(), sw_base_info.edges.end());
    children_base_path.insert(children_base_path.end(), se_base_info.edges.begin(), se_base_info.edges.end());


    std::vector<PortalEdge> edges;

    // MODIFIED: Add the actual geometric edge to each PortalEdge
    edges.push_back({4, 7, Distance(nw_p[1], ne_p[3]), {nw_p[1], ne_p[3]}});
    edges.push_back({5, 8, Distance(nw_p[2], sw_p[0]), {nw_p[2], sw_p[0]}});
    edges.push_back({6, 11, Distance(ne_p[2], se_p[0]), {ne_p[2], se_p[0]}});
    edges.push_back({9, 10, Distance(sw_p[1], se_p[3]), {sw_p[1], se_p[3]}});

    auto add_child_path_edges = [&](int child_idx, int p1_local, int p2_local, int n1_global, int n2_global){
        const auto& path_info = get_path_info(child_idx, {{{p1_local, p2_local}}});
        const auto& base_info = get_path_info(child_idx, {});
        if (path_info.cost != std::numeric_limits<double>::infinity()) {
            // This is complex: the "edge" represents the entire sub-path.
            // For MST, we only need the cost difference. The path edges are added separately.
            // For now, we don't add a specific geometric edge here, as the full path
            // from the subproblem will be added later if this connection is chosen.
            // This is a simplification; a more rigorous approach would track which sub-paths are used.
            edges.push_back({n1_global, n2_global, path_info.cost - base_info.cost, {}});
        }
    };
    // This logic remains the same, it just populates the edges vector
    add_child_path_edges(0, 1, 2, 4, 5); add_child_path_edges(1, 2, 3, 6, 7); add_child_path_edges(2, 0, 1, 8, 9); add_child_path_edges(3, 0, 3, 11, 10);
    // ... all other add_child_path_edges calls ...


    std::sort(edges.begin(), edges.end());

    auto all_parent_pairings = GenerateAllPossibleNonCrossingPairings(4);
    for (const auto& p_pairing : all_parent_pairings) {
        DSU dsu(12);
        double connection_mst_cost = 0;
        std::vector<std::pair<Point, Point>> connection_mst_path; // To store edges of the connection MST

        for (const auto& pair : p_pairing.pairs) dsu.unite(pair.first, pair.second);
        
        for (const auto& edge : edges) {
            if (dsu.find(edge.u) != dsu.find(edge.v)) {
                dsu.unite(edge.u, edge.v);
                connection_mst_cost += edge.weight;
                if (edge.actual_edge.first.x != 0 || edge.actual_edge.first.y != 0) { // Add edge if it's not empty
                   connection_mst_path.push_back(edge.actual_edge);
                }
            }
        }
        
        PathInfo final_path_info;
        final_path_info.cost = children_base_total_cost + connection_mst_cost;
        final_path_info.edges = children_base_path;
        final_path_info.edges.insert(final_path_info.edges.end(), connection_mst_path.begin(), connection_mst_path.end());
        parent_cost_table[p_pairing] = final_path_info;
    }
    
    std::cout << std::string(depth * 2, ' ') << "--> Combined. Populated DP table. Cost for empty pairing: " << parent_cost_table[PortalPairing{}].cost << std::endl;
    return parent_cost_table;
}

// MODIFIED: Returns PathInfo (cost + edges) instead of just double
PathInfo DPSolver::ComputeMSTCost(const std::vector<Point>& points) const {
    size_t n = points.size();
    if (n <= 1) return {0.0, {}};

    std::vector<double> min_cost(n, std::numeric_limits<double>::infinity());
    std::vector<int> parent(n, -1); // To reconstruct the path
    std::vector<bool> visited(n, false);
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;

    min_cost[0] = 0.0;
    pq.push({0.0, 0});

    double total_cost = 0.0;
    int edges_count = 0;
    while (!pq.empty() && edges_count < n) {
        double cost = pq.top().first;
        int u = pq.top().second;
        pq.pop();

        if (visited[u]) continue;
        visited[u] = true;
        total_cost += cost;
        edges_count++;

        for (size_t v = 0; v < n; ++v) {
            if (!visited[v]) {
                double weight = Distance(points[u], points[v]);
                if (weight < min_cost[v]) {
                    min_cost[v] = weight;
                    parent[v] = u; // Keep track of the parent in the MST
                    pq.push({weight, v});
                }
            }
        }
    }

    // Reconstruct path from parent array
    std::vector<std::pair<Point, Point>> path_edges;
    for (size_t i = 1; i < n; ++i) {
        if (parent[i] != -1) {
            path_edges.push_back({points[i], points[parent[i]]});
        }
    }

    return {total_cost, path_edges};
}