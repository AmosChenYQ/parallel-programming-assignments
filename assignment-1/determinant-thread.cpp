#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <typeinfo>
#include <vector>

#include "./hpc_helpers.hpp"

#define MAXN 8

template <typename value_t>
void print_vector(const std::vector<value_t>& v) {
  for (const auto& i : v) {
    std::cout << " " << i;
  }
  std::cout << std::endl;
}

// for given a permutaion sequence range from [0, size of sequence) with index
// [0, size of sequence), return the count of reverse pair in this sequence
uint64_t count_reverse_pair(const std::vector<uint64_t>& seq) {
  // fenwick tree is a array from [0, size of tree), denotes si in sequence
  // fenwick_tree[si] is the count of si appearances when construcing the tree
  std::vector<uint64_t> fenwick_tree(seq.size(), 0);
  // query prefix sum of fenwick tree from [0, index]
  auto query = [&](int64_t index) -> uint64_t {
    uint64_t sum = 0;

    while (index >= 0) {
      sum += fenwick_tree[index];
      index = (index & (index + 1)) - 1;
    }
    return sum;
  };
  // update element at index of fenwick tree with delta
  auto update = [&](int64_t index, const int64_t delta) -> void {
    while (index < static_cast<int64_t>(seq.size())) {
      fenwick_tree[index] += delta;
      index = index | (index + 1);
    }
  };

  // iterate sequence list reversly and create fenwick tree at the same time
  // according to defination, fenwick_tree[seq[idx] - 1] stands for the count of
  // numbers smaller than seq[idx] and at the right side of seq[idx]
  uint64_t reverse_pair_count = 0;
  for (auto reverse_iter = seq.rbegin(); reverse_iter != seq.rend();
       ++reverse_iter) {
    reverse_pair_count += query(static_cast<int64_t>(*reverse_iter) - 1);
    update(static_cast<int64_t>(*reverse_iter), 1);
  }
  return reverse_pair_count;
}

void init(std::vector<uint64_t>& A, std::vector<int64_t>& minors,
          uint64_t n) {
  // todo: add input to init matrix A here
  for (uint64_t row = 0; row < n; ++row) {
    for (uint64_t col = 0; col < n; ++col) {
      A[row * n + col] = 1;
    }
    minors[row] = 0;
  }
}

// for given a sequence indexing from [0, size) ranging from [0, size) like "3,
// 2, 1, 0" return the pai(3, 2, 1, 0) * A[0, 3] * A[1, 2] * A[2, 1] * A[3, 0]
// pai() is +1 when the count of reverser pair in sequence is even, otherwise
// pai() should return -1.
int64_t calculate_product(const std::vector<uint64_t>& A, uint64_t n,
                           const std::vector<uint64_t>& seq) {
  int64_t product = count_reverse_pair(seq) % 2 ? -1 : 1;
  for (uint64_t row = 0; row < n; ++row) {
    product *= A[row * n + seq[row]];
  }
  return product;
}

void calculate_determinant(const std::vector<uint64_t>& A,
                           std::vector<int64_t>& minors, uint64_t n,
                           uint64_t num_threads = 8) {
  const uint64_t chunk_size = SDIV(n, num_threads);

  auto break_down_by_col = [&](const uint64_t& id) -> void {
    const uint64_t lower = id * chunk_size;
    const uint64_t upper = std::min(n, lower + chunk_size);
    std::vector<int64_t> minors_segment(chunk_size, 0);

    for (uint64_t col = lower; col < upper; ++col) {
      uint64_t i = 0;
      std::vector<uint64_t> short_sequence(n - 1);
      for (uint64_t index = 0; index < n; ++index) {
        if (index != col) {
          short_sequence[i++] = index;
        }
      }
      do {
        std::vector<uint64_t> sequence(short_sequence);
        // construct generator vector standing for calcuate sequence
        sequence.insert(sequence.begin(), col);
        minors_segment[col - lower] += calculate_product(A, n, sequence);
      } while (
          std::next_permutation(short_sequence.begin(), short_sequence.end()));
    }

    // to prevent false sharing, move minor segment value into minor at the end
    // of each thread
    for (uint64_t col = lower; col < upper; ++col) {
      minors[col] = minors_segment[col - lower];
    }
  };

  std::vector<std::thread> threads;
  for (uint64_t id = 0; id < num_threads; id++) {
    threads.emplace_back(break_down_by_col, id);
  }
  for (auto&& thread : threads) {
    thread.join();
  }
}

int main(int argc, char* argv[]) {

  std::string s_dimensions(argv[1]);
  int dimensions = std::stoi(s_dimensions);
  uint64_t n = static_cast<uint64_t>(dimensions);

  std::string s_num_threads(argv[2]);
  int num_threads = std::stoi(s_num_threads);

  std::cout << "dimensions " << n << " threads " << num_threads << std::endl;

  TIMERSTART(overall)

  TIMERSTART(alloc)
  std::vector<uint64_t> A(n * n);
  std::vector<int64_t> minors(n);
  TIMERSTOP(alloc)

  TIMERSTART(init)
  init(A, minors, n);
  TIMERSTOP(init)

  TIMERSTART(calc)
  calculate_determinant(A, minors, n, num_threads);
  TIMERSTOP(calc)

  TIMERSTOP(overall)

  // check minors if they are correct
  for (uint64_t index = 0; index < minors.size(); ++index) {
    if(minors[index] != 0) {
      std::cout << "error at positon " << index << std::endl;
    }
  }

  return 0;
}
