#include <iostream>

#include <fstream>

#include <limits>

#include "Heap.h"

struct Node {
    uint32_t label;
    uint32_t dist;
};
bool operator==(const Node& lhs, const Node& rhs) {
    return lhs.dist == rhs.dist;
}
bool operator<(const Node& lhs, const Node& rhs) { return lhs.dist < rhs.dist; }
bool operator>(const Node& lhs, const Node& rhs) { return lhs.dist > rhs.dist; }
bool operator<=(const Node& lhs, const Node& rhs) {
    return lhs.dist <= rhs.dist;
}
bool operator>=(const Node& lhs, const Node& rhs) {
    return lhs.dist >= rhs.dist;
}

int main() {
    uint32_t node_count;

    std::ifstream fs("dj20.txt");
    if (!fs.is_open()) {
        std::cerr << "bleh\n";
        std::exit(EXIT_FAILURE);
    }
    fs >> node_count;

    std::vector<uint32_t> adj_matrix;
    adj_matrix.reserve(node_count * node_count);

    for (uint32_t j = 0; j < node_count; ++j) {
        for (uint32_t i = 0; i < node_count; ++i) {
            if (i < j + 1) {
                if (i < j) {
                    uint32_t index = i * node_count + j;
                    adj_matrix.push_back(adj_matrix[index]);
                } else {
                    adj_matrix.push_back(0);
                }
            } else {
                uint32_t v;
                fs >> v;
                adj_matrix.push_back(v);
            }
        }
    }

    for (uint32_t i = 0; i < node_count; ++i) {
        for (uint32_t j = 0; j < node_count; ++j) {
            uint32_t index = i * node_count + j;
            std::cout << adj_matrix[index] << "\t";
        }
        std::cout << "\n";
    }

    std::vector<Node> S;
    S.reserve(node_count);

    std::vector<Node> nodes;
    for (uint32_t i = 1; i <= node_count; ++i)
        nodes.push_back({i, std::numeric_limits<uint32_t>::max()});

    nodes[0].dist = 0;

    Heap<Node, HeapType::Min> Q(nodes);

    while (!Q.is_empty()) {
        for (auto& i : Q.get_array())
            std::cout << "(Node: " << i.label << " Dist: " << i.dist << ");";
        std::cout << "\n ";

        Node v = Q.extract();

        S.push_back(v);

        for (uint32_t i = 0; i < Q.get_size(); ++i) {
            Node u = Q.at(i);

            uint32_t matrix_index = (v.label - 1) * node_count + (u.label - 1);
            uint32_t adj = adj_matrix[matrix_index];

            if ((adj) && (u.dist > v.dist + adj)) {
                u.dist = v.dist + adj;
                Q.set(i, u);
            }
        }
    }

    std::cout << "\n";

    for (auto& i : S) {
        std::cout << "(Node: " << i.label << " Dist: " << i.dist << "); ";
    }
    std::cout << "\n";

    return 0;
}
