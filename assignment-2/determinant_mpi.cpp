#include <inttypes.h>

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "include/hpc_helpers.hpp"
#include "mpi.h"

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
int64_t count_reverse_pair(const std::vector<int64_t>& seq) {
  // fenwick tree is a array from [0, size of tree), denotes si in sequence
  // fenwick_tree[si] is the count of si appearances when construcing the tree
  std::vector<int64_t> fenwick_tree(seq.size(), 0);
  // query prefix sum of fenwick tree from [0, index]
  auto query = [&](int64_t index) -> int64_t {
    int64_t sum = 0;

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
  int64_t reverse_pair_count = 0;
  for (auto reverse_iter = seq.rbegin(); reverse_iter != seq.rend();
       ++reverse_iter) {
    reverse_pair_count += query(static_cast<int64_t>(*reverse_iter) - 1);
    update(static_cast<int64_t>(*reverse_iter), 1);
  }
  return reverse_pair_count;
}

void init(std::vector<int64_t>& A, std::vector<int64_t>& minors, int64_t n) {
  // todo: add input to init matrix A here
  for (int64_t row = 0; row < n; ++row) {
    for (int64_t col = 0; col < n; ++col) {
      A[row * n + col] = 1;
    }
    minors[row] = 0;
  }
}

// for given a sequence indexing from [0, size) ranging from [0, size) like "3,
// 2, 1, 0" return the pai(3, 2, 1, 0) * A[0, 3] * A[1, 2] * A[2, 1] * A[3, 0]
// pai() is +1 when the count of reverser pair in sequence is even, otherwise
// pai() should return -1.
int64_t calculate_product(const std::vector<int64_t>& A, int64_t n,
                          const std::vector<int64_t>& seq) {
  int64_t product = count_reverse_pair(seq) % 2 ? -1 : 1;
  for (int64_t row = 0; row < n; ++row) {
    product *= A[row * n + seq[row]];
  }
  return product;
}

void calculate_determinant(const std::vector<int64_t>& A,
                           std::vector<int64_t>& minors, int64_t n,
                           int num_processes, int my_id) {
  int64_t my_chunk_size = n / num_processes;

  if (my_id < n % num_processes) {
    my_chunk_size++;
  }

  std::vector<int64_t> minors_segment(my_chunk_size, 0);
  std::vector<int> send_counts(num_processes, 0);
  std::vector<int> displacements(num_processes, 0);

  if (!my_id) {
    for (int i = 0; i < num_processes; i++) {
      send_counts[i] =
          i < n % num_processes ? n / num_processes + 1 : n / num_processes;
    }
    for (int i = 1; i < num_processes; i++) {
      displacements[i] = displacements[i - 1] + send_counts[i - 1];
    }
    // for (int i = 0; i < num_processes; i++) {
    //   std::cout << send_counts[i] << " " << displacements[i] << std::endl;
    // }
  }

  int64_t lower = displacements[my_id];
  int64_t upper = lower + my_chunk_size;
  for (int64_t col = lower; col < upper; ++col) {
    int64_t i = 0;
    std::vector<int64_t> short_sequence(n - 1);
    for (int64_t index = 0; index < n; ++index) {
      if (index != col) {
        short_sequence[i++] = index;
      }
    }
    do {
      std::vector<int64_t> sequence(short_sequence);
      // construct generator vector standing for calcuate sequence
      sequence.insert(sequence.begin(), col);
      minors_segment[col - lower] += calculate_product(A, n, sequence);
    } while (
        std::next_permutation(short_sequence.begin(), short_sequence.end()));
  }

  MPI::COMM_WORLD.Gatherv(
      minors_segment.data(), my_chunk_size, MPI::LONG,
      minors.data(), send_counts.data(), displacements.data(), MPI::LONG, 0);
}

int main(int argc, char* argv[]) {
  MPI::Init(argc, argv);
  int num_processes = MPI::COMM_WORLD.Get_size();
  int my_id = MPI::COMM_WORLD.Get_rank();

  std::string s_dimensions(argv[1]);
  int dimensions = std::stoi(s_dimensions);
  int64_t n = static_cast<int64_t>(dimensions);

  MPI::COMM_WORLD.Barrier();
  double start;
  if (!my_id) {
    std::cout << "dimensions " << n << " processes " << num_processes
              << std::endl;
    start = MPI::Wtime();
  }

  std::vector<int64_t> A(n * n);
  std::vector<int64_t> minors(n);
  init(A, minors, n);

  calculate_determinant(A, minors, n, num_processes, my_id);

  MPI::COMM_WORLD.Barrier();
  double end;

  // check minors if they are correct
  if (!my_id) {
    for (auto index = 0u; index < minors.size(); ++index) {
      if (minors[index] != 0) {
        std::cout << "error at positon " << index << std::endl;
      }
    }
    end = MPI::Wtime();
    std::cout << "Used " << end - start << " seconds" << std::endl;
  }

  MPI::Finalize();
  return 0;
}
