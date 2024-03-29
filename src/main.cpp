#include <iostream>

#include <algorithm>
#include <fstream>
#include <limits>
#include <list>
#include <vector>

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

        pos[nodes[size - 1].id] = 0;
        pos[extracted.id] = size - 1;

        std::swap(nodes[0], nodes[size - 1]);
        --size;

        heapify(0);
        return extracted;
    }

    void heapify(uint32_t index) {
        uint32_t t = index;

        while (true) {
            uint32_t l = left_child(t);
            uint32_t r = right_child(t);

            if ((l < size) && (nodes[t].dist > nodes[l].dist)) t = l;
            if ((r < size) && (nodes[t].dist > nodes[r].dist)) t = r;

            if (t != index) {
                pos[nodes[t].id] = index;
                pos[nodes[index].id] = t;

                std::swap(nodes[index], nodes[t]);
                index = t;
            } else {
                break;
            }
        }
    }

    bool decrement(uint32_t id, uint32_t dist) {
        uint32_t index = pos[id];

        if ((index >= size) || (nodes[index].dist < dist)) return false;

        nodes[index].dist = dist;

        while (index && (nodes[index].dist < nodes[parent(index)].dist)) {
            pos[nodes[index].id] = parent(index);
            pos[nodes[parent(index)].id] = index;

            uint32_t p_id = nodes[parent(index)].id;

            std::swap(nodes[index], nodes[parent(index)]);

            index = pos[p_id];
        }

        return true;
    }

    std::vector<HeapNode> get_nodes() const {
        return std::vector<HeapNode>(nodes.begin(), (nodes.begin() + size));
    }
    bool is_empty() const { return size == 0; }

   private:
    std::vector<uint32_t> pos;
    std::vector<HeapNode> nodes;

    uint32_t size;

    uint32_t parent(uint32_t index) const { return index / 2.0f; }
    uint32_t left_child(uint32_t index) const { return 2 * index + 1; }
    uint32_t right_child(uint32_t index) const { return 2 * index + 2; }
};

struct GraphNode {
    uint32_t id;
    uint32_t p;
    std::vector<uint32_t> adjs;
};

class Graph {
   public:
    Graph() = default;

    void open_matrix(const std::string& file_path) {
        uint32_t node_count = 0;

        std::ifstream fs(file_path);
        if (!fs.is_open()) {
            std::cerr << "Failed to open file: " << file_path;
            std::exit(EXIT_FAILURE);
        }
        fs >> node_count;

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

        graph_nodes.reserve(node_count);
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
    }

    std::pair<std::vector<GraphNode>, uint32_t> dijkstra() {
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
        S.reserve(graph_nodes.size());

        MinHeap Q(heap_nodes);

        while (!Q.is_empty()) {
            HeapNode v = Q.extract();
            S.push_back(v);

            if (v.id == graph_nodes.size() - 1) break;

            GraphNode gp = graph_nodes[v.id];

            for (uint32_t i = 0; i < gp.adjs.size(); ++i) {
                uint32_t p = gp.adjs.at(i);

                uint32_t adj_dist = adj_matrix[gp.id * graph_nodes.size() + p];

                if (Q.decrement(p, v.dist + adj_dist)) graph_nodes[p].p = v.id;
            }
        }

        std::vector<GraphNode> st;
        GraphNode u = graph_nodes[S.back().id];
        while (true) {
            st.push_back(u);
            if (u.p < std::numeric_limits<uint32_t>::max()) {
                u = graph_nodes[u.p];
            } else {
                break;
            }
        }

        std::vector<GraphNode> path;
        path.reserve(st.size());
        while (!st.empty()) {
            GraphNode u = st.back();
            path.push_back(u);
            st.pop_back();
        }

        return {path, S.back().dist};
    }

   private:
    std::vector<GraphNode> graph_nodes;
    std::vector<uint32_t> adj_matrix;
};

int main(int argc, char** argv) {
    uint32_t node_count;

    if (argc != 2) {
        std::cout << "dijkstra <file_path>" << std::endl;
        std::exit(EXIT_SUCCESS);
    }

    std::string file_path = argv[1];

    Graph gp;
    gp.open_matrix(file_path);
    auto result = gp.dijkstra();

    std::cout << "Solucao: " << result.second << "\n";
    for (auto& n : result.first) std::cout << "{ Node: " << n.id << " } ";
    std::cout << std::endl;

    return 0;
}
