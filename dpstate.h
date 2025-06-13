#pragma once

#include "geometry.h"
#include <vector>
#include <map>
#include <utility> // For std::pair

// Represents a non-crossing pairing of portals on a square's boundary.
// Portals are indexed 0 to m-1. A pairing is a set of pairs (i, j).
// We sort the vector of pairs to make it a canonical key for a map.
struct PortalPairing {
    std::vector<std::pair<int, int>> pairs;
    
    bool operator<(const PortalPairing& other) const {
        return pairs < other.pairs;
    }

    bool operator==(const PortalPairing& other) const {
        return pairs == other.pairs;
    }
};
// The complete DP information for a single square.
// It maps each possible portal pairing to its minimum cost.
using DPCostTable = std::map<PortalPairing, double>;

// The main cache for our memoization (the DP table).
// It maps a square's identifier to its DPCostTable.
using DPCache = std::map<BoundingSquare, DPCostTable>;