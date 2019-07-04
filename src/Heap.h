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
    Heap() = default;
    Heap(const std::vector<T>& array_) : array(array_) {
        size = array.size();
        build();
    }
    ~Heap() {}

    T extract() {
        T extracted = array[0];
        array[0] = array[(size--) - 1];
        heapify(0);

        return extracted;
    }
    void heapify(uint32_t index) {
        while (true) {
            uint32_t t = index;
            uint32_t l = left_child(t);
            uint32_t r = right_child(t);

            if constexpr (H == HeapType::Min) {
                if (l < size && (array[t] > array[l])) t = l;
                if (r < size && (array[t] > array[r])) t = r;
            } else {
                if (l < size && array[t] < array[l]) t = l;
                if (r < size && array[t] < array[r]) t = r;
            }

            if (t != index) {
                std::swap(array[index], array[t]);
            } else {
                break;
            }
        }
    }
    void set(uint32_t index, T v) {
        if (index >= size) {
            std::cerr << "Heap: out of bounds" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if constexpr (H == HeapType::Min) {
            if (v == array[index]) {
                return;
            } else if (v < array[index]) {
                array[index] = v;
                while ((index != 0) && (array[index] < array[parent(index)])) {
                    std::swap(array[index], array[parent(index)]);
                    index = parent(index);
                }
            } else {
                array[index] = v;
                heapify(index);
            }
        } else {
            if (v == array[index]) {
                return;
            } else if (v > array[index]) {
                array[index] = v;
                while ((index != 0) && (array[index] > array[parent(index)])) {
                    std::swap(array[index], array[parent(index)]);
                    index = parent(index);
                }
            } else {
                array[index] = v;
                heapify(index);
            }
        }
    }

    const T& at(uint32_t index) const {
        if (index < size) {
            return array[index];
        } else {
            std::cerr << "Heap: out of bounds" << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    std::vector<T> get_array() const {
        return std::vector<T>(array.begin(), (array.begin() + size));
    }
    T top() const { return array[0]; }
    bool is_empty() const { return size == 0; }
    uint32_t get_size() const { return size; }

   private:
    std::vector<T> array;
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