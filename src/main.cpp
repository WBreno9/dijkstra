#include <iostream>

#include <algorithm>
#include <fstream>
#include <limits>
#include <vector>

struct Node {
    uint32_t id;
    uint32_t dist;
};

class MinHeap {
   public:
    MinHeap() = default;
    MinHeap(const std::vector<Node>& array_) : array(array_) {
        uint32_t k = 0;
        for (auto& n : array) n.id = k++;

        pos.resize(array.size());
        k = 0;
        for (auto& p : pos) p = k++;

        size = array.size();
        build();
    }
    ~MinHeap() {}

    Node extract() {
        Node extracted = array[0];
        array[0] = array[(size--) - 1];
        pos[0] = size - 1;
        heapify(0);

        return extracted;
    }
    void heapify(uint32_t index) {
        while (true) {
            uint32_t t = index;
            uint32_t l = left_child(t);
            uint32_t r = right_child(t);

            if (l < size && (array[t].dist > array[l].dist)) t = l;
            if (r < size && (array[t].dist > array[r].dist)) t = r;

            if (t != index) {
                pos[array[t].id] = index;
                pos[array[index].id] = t;

                std::swap(array[index], array[t]);
            } else {
                break;
            }
        }
    }
    void decrease(Node& v, uint32_t new_dist) {
        uint32_t index = pos[v.id];

        v.dist = new_dist;
        array[index] = v;

        while ((index) && (array[index].dist < array[parent(index)].dist)) {
            pos[array[index].id] = parent(index);
            pos[array[parent(index)].id] = index;

            std::swap(array[index], array[parent(index)]);
            index = parent(index);
        }
    }

    std::vector<Node> get_array() const {
        return std::vector<Node>(array.begin(), (array.begin() + size));
    }

    bool is_empty() const { return size == 0; }

    std::vector<uint32_t> pos;
    std::vector<Node> array;
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

    std::vector<Node> S;
    S.reserve(node_count);

    std::vector<Node> nodes;
    for (uint32_t i = 0; i < node_count; ++i)
        nodes.push_back({i, std::numeric_limits<uint32_t>::max()});
    nodes[0].dist = 0;

    std::vector<uint32_t> P(nodes.size(), std::numeric_limits<uint32_t>::max());

    MinHeap Q(nodes);
    while (!Q.is_empty()) {
        Node v = Q.extract();
        S.push_back(v);

        for (uint32_t i = 0; i < node_count; ++i) {
            Node& u = nodes[i];

            uint32_t matrix_index = v.id * node_count + i;
            uint32_t adj = adj_matrix[matrix_index];

            if (u.dist > v.dist + adj) {
                uint32_t new_dist = v.dist + adj;
                P[u.id] = v.id;
                Q.decrease(u, new_dist);
            }
        }
    }

    for (auto& u : S)
        std::cout << "{Node = " << u.id << ", Dist = " << u.dist << "} ";
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
