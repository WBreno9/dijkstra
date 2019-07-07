#include <iostream>

#include <algorithm>
#include <fstream>
#include <limits>
#include <list>
#include <vector>

struct GraphNode {
    uint32_t id;
    uint32_t p;
    std::vector<uint32_t> adjs;
};

struct HeapNode {
    uint32_t id;
    uint32_t dist;
};

class MinHeap {
   public:
    MinHeap(const std::vector<HeapNode>& nodes_) : nodes(nodes_) {
        pos.resize(nodes.size());
        uint32_t k = 0;
        for (auto& p : pos) p = k++;

        size = nodes.size();

        uint32_t j = (size / 2.0f) - 1;
        for (uint32_t i = j + 1; i != 0; --i) heapify(i - 1);
    }

    HeapNode extract() {
        HeapNode extracted = nodes[0];

        // pos[nodes[0].id] = 0;
        // pos[extracted.id] = size - 1;

        nodes[0] = nodes[(size--) - 1];

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
                // pos[nodes[t].id] = index;
                // pos[nodes[index].id] = t;

                std::swap(nodes[index], nodes[t]);
            } else {
                break;
            }
        }
    }

    bool decrement(uint32_t id, uint32_t dist) {
        uint32_t index;

        for (uint32_t i = 0; i < nodes.size(); ++i) {
            HeapNode& n = nodes[i];
            if (n.id == id) index = i;
        }

        if (index >= size) return false;

        if (nodes[index].dist < dist) return false;

        nodes[index].dist = dist;

        while (index && (nodes[index].dist < nodes[parent(index)].dist)) {
            // pos[nodes[index].id] = parent(index);
            // pos[nodes[parent(index)].id] = index;

            uint32_t p_id = nodes[parent(index)].id;

            std::swap(nodes[index], nodes[parent(index)]);

            id = p_id;
            for (uint32_t i = 0; i < nodes.size(); ++i) {
                HeapNode& n = nodes[i];
                if (n.id == id) index = i;
            }
        }

        return true;
    }

    std::vector<HeapNode> get_nodes() const {
        return std::vector<HeapNode>(nodes.begin(), (nodes.begin() + size));
    }
    bool is_empty() const { return size == 0; }

    // private:
    std::vector<uint32_t> pos;
    std::vector<HeapNode> nodes;

    uint32_t size;

    uint32_t parent(uint32_t index) const { return index / 2.0f; }
    uint32_t left_child(uint32_t index) const { return 2 * index + 1; }
    uint32_t right_child(uint32_t index) const { return 2 * index + 2; }
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

    // std::ofstream ofs("tst");
    //
    // ofs << "\t";
    // for (uint32_t j = 0; j < node_count; ++j) ofs << j << "\t";
    // ofs << "\n";
    //
    // ofs << "\t";
    //
    // for (uint32_t j = 0; j < node_count; ++j) {
    //    ofs << j << "\t";
    //    for (uint32_t i = 0; i < node_count; ++i) {
    //        uint32_t matrix_index = j * node_count + i;
    //        ofs << adj_matrix[matrix_index] << "\t";
    //    }
    //    ofs << "\n";
    //}
    // ofs.close();

    std::vector<GraphNode> graph_nodes;
    for (uint32_t i = 0; i < node_count; ++i) {
        std::vector<uint32_t> adj;
        adj.reserve(node_count);

        for (uint32_t j = 0; j < node_count; ++j) {
            uint32_t matrix_index = i * node_count + j;
            if (adj_matrix[matrix_index]) adj.push_back(j);
        }

        GraphNode gnode;
        gnode.id = i;
        gnode.p = std::numeric_limits<uint32_t>::max();
        gnode.adjs = adj;

        graph_nodes.push_back(gnode);
    }

    std::vector<HeapNode> heap_nodes;
    heap_nodes.reserve(graph_nodes.size());

    for (uint32_t i = 0; i < graph_nodes.size(); ++i) {
        HeapNode n;
        n.id = i;
        n.dist = std::numeric_limits<uint32_t>::max();

        heap_nodes.push_back(n);
    }
    heap_nodes[0].dist = 0;

    std::vector<HeapNode> S;
    S.reserve(node_count);

    for (uint32_t i = 1; i < node_count; ++i) {
        heap_nodes[i].dist = adj_matrix[i];
    }

    for (auto& u : heap_nodes)
        std::cout << "{Node = " << u.id << ", Dist = " << u.dist << ", P = "
                  << "} ";
    std::cout << "\n";
    std::cout << "\n";

    MinHeap Q(heap_nodes);

    while (!Q.is_empty()) {
        HeapNode v = Q.extract();
        S.push_back(v);
    }

    for (auto& u : S)
        std::cout << "{Node = " << u.id << ", Dist = " << u.dist << ", P = "
                  << "} ";
    std::cout << "\n";

    /*

    while (!Q.is_empty()) {
        HeapNode v = Q.extract();
        S.push_back(v);

        GraphNode gp = graph_nodes[v.id];

        for (uint32_t i = 0; i < gp.adjs.size(); ++i) {
            uint32_t p = gp.adjs.at(i);

            uint32_t adj_dist = adj_matrix[gp.id * node_count + p];

            Q.decrement(p, v.dist + adj_dist);
        }
    }


    for (auto& u : S)
        std::cout << "{Node = " << u.id << ", Dist = " << u.dist << ", P = "
                  << "} ";
    std::cout << "\n";
    */

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
