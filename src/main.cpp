#include <iostream>

#include <algorithm>
#include <fstream>
#include <limits>

#include "Heap.h"

struct Node {
    uint32_t id;
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

int main(int argc, char** argv) {
    uint32_t node_count;

    if (argc != 2) {
        std::cout << "dijkstra <file_path>" << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    std::string file_path = argv[1];

    std::ifstream fs(file_path);
    if (!fs.is_open()) {
        std::cerr << "Failed to open file: " << file_path;
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

    std::vector<Node> S;
    S.reserve(node_count);

    std::vector<Node> nodes;
    for (uint32_t i = 0; i < node_count; ++i)
        nodes.push_back({i, std::numeric_limits<uint32_t>::max()});
    nodes[0].dist = 0;

    std::vector<uint32_t> P(nodes.size(), std::numeric_limits<uint32_t>::max());

    MinHeap<Node> Q(nodes);
    while (!Q.is_empty()) {
        Node v = Q.extract();
        S.push_back(v);

        for (uint32_t i = 0; i < Q.get_size(); ++i) {
            Node u = Q.at(i);

            uint32_t matrix_index = v.id * node_count + u.id;
            uint32_t adj = adj_matrix[matrix_index];

            if ((adj) && (u.dist > v.dist + adj)) {
                u.dist = v.dist + adj;
                nodes[u.id].dist = u.dist;
                P[u.id] = v.id;
                Q.set(i, u);
            }
        }
    }

    for (auto& u : S) {
        std::cout << "{Node = " << u.id << ", Dist = " << u.dist << "} ";
    }
    std::cout << "\n";

    S.clear();
    Node u = nodes.back();

    while (true) {
        S.push_back(u);
        if (P[u.id] < std::numeric_limits<uint32_t>::max()) {
            u = nodes[P[u.id]];
        } else {
            break;
        }
    }

    while (!S.empty()) {
        Node u = S.back();
        S.pop_back();
        std::cout << "{Node = " << u.id << ", Dist = " << u.dist << "} ";
    }
    std::cout << "\n";

    std::cout << "Solucao: " << nodes.back().dist << std::endl;

    return 0;
}
