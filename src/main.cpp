#include <iostream>

#include <algorithm>
#include <fstream>
#include <limits>
#include <list>

struct HeapNode {
    uint32_t id;
    uint32_t dist;
    std::vector<uint32_t> adjs;  // not optimal, has to move a vector
    uint32_t p;
};

class MinHeap {
   public:
    MinHeap(const std::vector<HeapNode>& nodes_) : nodes(nodes_) {
        uint32_t k = 0;
        for (auto& n : nodes) n.id = k++;

        pos.resize(nodes.size());
        k = 0;
        for (auto& p : pos) p = k++;

        size = nodes.size();
        build();
    }
    ~MinHeap() {}

    HeapNode extract() {
        HeapNode extracted = nodes[0];
        nodes[0] = nodes[size - 1];
        pos[0] = size - 1;
        --size;

        heapify(0);
        return extracted;
    }

    void heapify(uint32_t index) {
        while (true) {
            uint32_t t = index;
            uint32_t l = left_child(t);
            uint32_t r = right_child(t);

            if (l < size && (nodes[t].dist > nodes[l].dist)) t = l;
            if (r < size && (nodes[t].dist > nodes[r].dist)) t = r;

            if (t != index) {
                pos[nodes[t].id] = index;
                pos[nodes[index].id] = t;

                std::swap(nodes[index], nodes[t]);
            } else {
                break;
            }
        }
    }

    void decrement(HeapNode n, uint32_t dist) {
        uint32_t index = pos[n.id];

        if (dist == nodes[index].dist) {
            return;
        } else if (dist < nodes[index].dist) {
            nodes[index].dist = dist;
            nodes[index].p = n.p;

            while (index && (nodes[index].dist < nodes[parent(index)].dist)) {
                pos[nodes[index].id] = parent(index);
                pos[nodes[parent(index)].id] = index;

                std::swap(nodes[index], nodes[parent(index)]);

                index = parent(index);
            }
        }
    }

    const HeapNode& at(uint32_t index) const {
        if (index < size) {
            return nodes[index];
        } else {
            std::cerr << "Heap (at): out of bounds" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    uint32_t get_pos(uint32_t index) const { return pos[index]; }

    std::vector<HeapNode> get_nodes() const {
        return std::vector<HeapNode>(nodes.begin(), (nodes.begin() + size));
    }
    bool is_empty() const { return size == 0; }
    uint32_t get_size() const { return size; }

    // private:
    std::vector<HeapNode> nodes;
    std::vector<uint32_t> pos;
    uint32_t size;

    uint32_t parent(uint32_t index) const { return index / 2.0f; }
    uint32_t left_child(uint32_t index) const { return 2 * index + 1; }
    uint32_t right_child(uint32_t index) const { return 2 * index + 2; }

    void build() {
        uint32_t j = (size / 2.0f) - 1;
        for (uint32_t i = j + 1; i != 0; --i) heapify(i - 1);
    }
};

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
    fs.close();

    std::ofstream ofs("tst");

    ofs << "\t";
    for (uint32_t j = 0; j < node_count; ++j) ofs << j << "\t";
    ofs << "\n";

    // ofs << "\t";

    for (uint32_t j = 0; j < node_count; ++j) {
        ofs << j << "\t";
        for (uint32_t i = 0; i < node_count; ++i) {
            uint32_t matrix_index = j * node_count + i;
            ofs << adj_matrix[matrix_index] << "\t";
        }
        ofs << "\n";
    }
    ofs.close();

    std::vector<HeapNode> S;
    S.reserve(node_count);

    std::vector<HeapNode> nodes;
    for (uint32_t i = 0; i < node_count; ++i) {
        std::vector<uint32_t> adj;
        adj.reserve(node_count);

        for (uint32_t j = 0; j < node_count; ++j) {
            uint32_t matrix_index = i * node_count + j;
            if (adj_matrix[matrix_index]) adj.push_back(j);
        }

        nodes.push_back({i, std::numeric_limits<uint32_t>::max(), adj,
                         std::numeric_limits<uint32_t>::max()});
    }
    nodes[0].dist = 0;

    MinHeap Q(nodes);
    while (!Q.is_empty()) {
        HeapNode v = Q.extract();
        S.push_back(v);

        // std::cout << "Q: ";
        // for (auto& u : Q.get_nodes())
        //    std::cout << "{Node = " << u.id << ", Dist = " << u.dist << "} ";
        // std::cout << "\n";

        for (uint32_t i = 0; i < v.adjs.size(); ++i) {
            uint32_t p = v.adjs.at(i);
            HeapNode u = nodes[p];

            uint32_t matrix_index = v.id * node_count + u.id;
            uint32_t adj_dist = adj_matrix[matrix_index];

            if (u.dist > v.dist + adj_dist) {
                u.p = v.id;
                Q.decrement(u, v.dist + adj_dist);
            }
        }
    }

    for (auto& u : S)
        std::cout << "{Node = " << u.id << ", Dist = " << u.dist
                  << ", P = " << u.p << "} ";
    std::cout << "\n";

    /*
    S.clear();
    HeapNode u = S.back();

    while (true) {
        S.push_back(u);
        if (u.p < std::numeric_limits<uint32_t>::max()) {
            u = Q.nodes[u.p];
        } else {
            break;
        }
    }

    std::cout << "\n";

    while (!S.empty()) {
        HeapNode u = S.back();
        S.pop_back();
        std::cout << "{Node = " << u.id << ", Dist = " << u.dist << "} ";
    }
    std::cout << "\n";*/

    return 0;
}
