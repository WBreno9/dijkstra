#pragma once

#include <iostream>

#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

enum class HeapType { Min, Max };

template <typename T, HeapType H>
class Heap {
   public:
    Heap(const std::vector<T>& array) {
        nodes.reserve(array.size());
        for (uint32_t i = 0; i < array.size(); ++i)
            nodes.push_back({i, array[i]});

        pos.resize(nodes.size());
        for (uint32_t i = 0; i < nodes.size(); ++i) pos[i] = i;

        size = nodes.size();
        build();
    }
    ~Heap() {}

    struct HeapNode {
        uint32_t id;
        T d;
    };

    T extract() {
        HeapNode extracted = nodes[0];
        nodes[0] = nodes[size - 1];
        pos[0] = size - 1;
        --size;

        heapify(0);
        return extracted.d;
    }

    void heapify(uint32_t index) {
        while (true) {
            uint32_t t = index;
            uint32_t l = left_child(t);
            uint32_t r = right_child(t);

            if constexpr (H == HeapType::Min) {
                if (l < size && (nodes[t].d > nodes[l].d)) t = l;
                if (r < size && (nodes[t].d > nodes[r].d)) t = r;
            } else {
                if (l < size && nodes[t].d < nodes[l].d) t = l;
                if (r < size && nodes[t].d < nodes[r].d) t = r;
            }

            if (t != index) {
                pos[nodes[t].id] = index;
                pos[nodes[index].id] = t;

                std::swap(nodes[index], nodes[t]);
            } else {
                break;
            }
        }
    }

    /*
    void set(uint32_t index, T v) {
        if (index >= size) {
            std::cerr << "Heap (set): out of bounds - " << index << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if constexpr (H == HeapType::Min) {
            if (v == nodes[index]) {
                return;
            } else if (v < nodes[index]) {
                nodes[index] = v;
                while ((index != 0) && (nodes[index] < nodes[parent(index)])) {
                    uint32_t tmp = pos[index];
                    uint32_t pp = parent(index);

                    pos[index] = pp;
                    pos[pos[pp]] = tmp;

                    std::swap(nodes[index], nodes[parent(index)]);

                    // std::swap(pos[index], pos[parent(index)]);
                    index = parent(index);
                }
            } else {
                nodes[index] = v;
                heapify(index);
            }
        } else {
            if (v == nodes[index]) {
                return;
            } else if (v > nodes[index]) {
                nodes[index] = v;
                while ((index != 0) && (nodes[index] > nodes[parent(index)])) {
                    std::swap(nodes[index], nodes[parent(index)]);
                    std::swap(pos[index], pos[parent(index)]);
                    index = parent(index);
                }
            } else {
                nodes[index] = v;
                heapify(index);
            }
        }
    }*/

    const T& at(uint32_t index) const {
        if (index < size) {
            return nodes[index].d;
        } else {
            std::cerr << "Heap (at): out of bounds" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    uint32_t get_pos(uint32_t index) const { return pos[index]; }

    std::vector<T> get_nodes() const {
        return std::vector<T>(nodes.begin(), (nodes.begin() + size));
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

template <typename V, typename K, HeapType H>
using AssociativeHeap = Heap<std::pair<V, K>, H>;
template <typename T>
using MinHeap = Heap<T, HeapType::Min>;
template <typename T>
using MaxHeap = Heap<T, HeapType::Max>;