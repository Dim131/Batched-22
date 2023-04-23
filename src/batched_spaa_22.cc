/* Code for producing the experiments in Section 8 of 
      "Balanced Allocations in Batches: Simplified and Generalized"
      by Dimitrios Los and Thomas Sauerwald (SPAA'22)
      [https://arxiv.org/abs/2203.13902]. */
#include <algorithm>
#include <functional>
#include <iostream>
#include <random>

/* Runs the b-Batched setting. This process was introduced in
*     "Multiple-choice balanced allocation in (almost) parallel",
         by Berenbrink, Czumaj, Englert, Friedetzky, and Nagel (2012)
         [https://arxiv.org/abs/1501.04822].

   It starts from an empty load vector and in each round:
     - Allocates b (potentially weighted) balls using the process provided, 
       using the load information at the beginning of the batch.

   This class keeps track of the load-vector, the maximum load and gap.
   */
template<typename LoadType, typename Generator>
class BatchedSetting {
public:

  /* Initializes the b-Batched setting, given:
      - Number of bins (n)
      - Batch size (b)
      - Function that samples a bin given the load vector
      - Function that samples the weight of a ball. */
  BatchedSetting(
    size_t num_bins,
    size_t batch_size,
    const std::function<size_t(const std::vector<LoadType>&, Generator&)> bin_selector,
    const std::function<LoadType(Generator&)> weight_generator)
    : load_vector_(num_bins, 0), buffer_vector_(num_bins, 0), batch_size_(batch_size), max_load_(0), total_weight_(0),
    bin_selector_(bin_selector), weight_generator_(weight_generator) {

  }

  /* Performs an allocation of a batch. */
  void nextRound(Generator& generator) {
    // Phase 1: Perform b allocations.
    size_t n = load_vector_.size();

    for (size_t i = 0; i < batch_size_; ++i) {
      size_t idx = bin_selector_(load_vector_, generator);
      LoadType ball_weight = weight_generator_(generator);
      buffer_vector_[idx] += ball_weight;
      total_weight_ += ball_weight;
    }

    // Phase 2: Update and sort the load vector.
    for (int i = 0; i < n; ++i) {
      load_vector_[i] += buffer_vector_[i];
      buffer_vector_[i] = 0;
      max_load_ = std::max(max_load_, load_vector_[i]);
    }
    std::sort(load_vector_.begin(), load_vector_.end(), std::greater<LoadType>());
  }

  /* Returns the current maximum load. */
  LoadType getMaxLoad() const {
    return max_load_;
  }

  /* Returns the current gap. */
  LoadType getGap() const {
    return max_load_ - total_weight_ / LoadType(load_vector_.size());
  }

  /* Returns the current load vector. */
  std::vector<size_t> getLoadVector() const {
    return load_vector_;
  }

private:

  /* Current load vector of the process. */
  std::vector<LoadType> load_vector_;

  /* Buffer vector for the balls allocated in the current batch. */
  std::vector<LoadType> buffer_vector_;

  /* Batch size used in the setting. */
  const size_t batch_size_;

  /* Current maximum load in the load vector. */
  LoadType max_load_;

  /* Total weight of balls in the load vector. */
  LoadType total_weight_;

  /* Function that selects the bin. */
  std::function<size_t(const std::vector<LoadType>&, Generator&)> bin_selector_;

  /* Function that generates the weight of the ball. */
  std::function<LoadType(Generator&)> weight_generator_;
};

/* Sample a bin via the d-Choice process. */
template<class LoadType, class Generator, int d>
size_t d_choice(const std::vector<LoadType>& load_vector, Generator& generator) {
  size_t n = load_vector.size();
  std::uniform_int<size_t> uar(0, n - 1);
  size_t sample = uar(generator);
  for (int i = 1; i < d; ++i) {
    sample = std::max(sample, uar(generator));
  }
  return sample;
}

/* Sample a bin via the Two-Choice process with random tie-breaking. */
template<typename LoadType, typename Generator>
size_t two_choice_with_random_tie_breaking(const std::vector<LoadType>& load_vector, Generator& generator) {
  size_t n = load_vector.size();
  std::uniform_int<size_t> uar(0, n - 1);
  size_t i1 = uar(generator), i2 = uar(generator);
  // Break ties randomly. 
  if (load_vector[i1] <= load_vector[i2]) return i1;
  return i2;
}

/* Sample a bin via the (1+beta)-process. */
template<typename LoadType, typename Generator>
std::function<size_t(const std::vector<LoadType>&, Generator&)> one_plus_beta(double beta) {
  return [beta](const std::vector<LoadType>& load_vector, Generator& generator) {
    size_t n = load_vector.size();
    std::bernoulli_distribution use_two_choices(beta);
    std::uniform_int<size_t> uar(0, n - 1);
    size_t i1 = uar(generator);
    if (!use_two_choices(generator)) {
      return i1;
    }
    size_t i2 = uar(generator);
    return i1 < i2 ? i2 : i1;
  };
}

/* Sample a bin via the (1+beta)-process with random tie-breaking. */
template<typename LoadType, typename Generator>
std::function<size_t(const std::vector<LoadType>&, Generator&)> one_plus_beta_with_random_tie_breaking(double beta) {
  return [beta](const std::vector<LoadType>& load_vector, Generator& generator) {
    size_t n = load_vector.size();
    std::bernoulli_distribution use_two_choices(beta);
    std::uniform_int<size_t> uar(0, n - 1);
    size_t i1 = uar(generator);
    if (!use_two_choices(generator)) {
      return i1;
    }
    size_t i2 = uar(generator);
    if (load_vector[i1] <= load_vector[i2]) return i1;
    return i2;
  };
}

/* Sample a bin selector for the hybrid of the Quantile and 
   the Two-Choice process. */
template<typename LoadType, typename Generator>
std::function<size_t(const std::vector<LoadType>&, Generator&)> hybrid_quantile(double delta) {
  return [delta](const std::vector<LoadType>& load_vector, Generator& generator) {
    size_t n = load_vector.size();
    std::uniform_int<size_t> uar(0, n - 1);
    size_t i1 = uar(generator);
    if (i1 >= delta * n) {
      return i1;
    }
    return uar(generator);
  };
}

/* Generate unit weight balls. */
template<typename Generator>
size_t unit_weights_generator(Generator& generator) {
  return 1;
}

/* Generate balls with weights from an Exponential distribution. */
template<typename Generator>
double exp_weights_generator(Generator& generator) {
  std::exponential_distribution<double> exp_distribution(1.0);
  return exp_distribution(generator);
}

/* Returns the average (maximum) gap for the given process. */
template<typename LoadType = size_t, typename Generator = std::mt19937>
double avg_gap(
  int num_bins,
  int normalised_batch_size,
  int num_batches,
  int repetitions,
  const std::function<size_t(const std::vector<LoadType>&, Generator&)> process,
  const std::function<LoadType(Generator&)> weights_generator = unit_weights_generator<Generator>) {

  Generator generator;
  double gap_total = 0.0;
  for (int i = 0; i < repetitions; ++i) {
    BatchedSetting<LoadType, std::mt19937> batched_setting(
      num_bins, normalised_batch_size * num_bins, process, weights_generator);
    LoadType cur_gap = 0;
    for (int round = 0; round < num_batches; ++round) {
      batched_setting.nextRound(generator);
      cur_gap = std::max(cur_gap, batched_setting.getGap());
    }
    gap_total += cur_gap;
  }
  return gap_total / double(repetitions);
}

/* Auxiliary function for generating Figures 8.1 and 8.3.*/
template<typename LoadType, typename Generator = std::mt19937>
void generate_several_processes_vs_batch_size(const std::function<LoadType(Generator&)> weights_generator) {
  int num_bins = 1'000;
  int repetitions = 50;
  std::cout << "Three-Choice:" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<LoadType, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions,
        d_choice<LoadType, Generator, 3>, weights_generator) << ")" << std::endl;
  }

  std::cout << "Two-Choice:" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<LoadType, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions, 
        two_choice_with_random_tie_breaking<LoadType, Generator>, weights_generator) << ")" << std::endl;
  }

  std::cout << "(1+beta) with beta=0.7:" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<LoadType, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions, 
        one_plus_beta_with_random_tie_breaking<LoadType, Generator>(0.7), weights_generator) << ")" << std::endl;
  }

  std::cout << "(1+beta) with beta=0.5:" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 50; ++normalised_batch_size) {
    int num_rounds = num_bins / normalised_batch_size;
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<LoadType, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions,
        one_plus_beta_with_random_tie_breaking<LoadType, Generator>(0.5), weights_generator) << ")" << std::endl;
  }
}

/* Figure 8.1: Average gap for various processes for b = {n, 2n, ... , 50n}
   at m = n^2 balls. */
void generate_several_processes_vs_batch_size_unit_weights() {
  generate_several_processes_vs_batch_size<size_t, std::mt19937>(unit_weights_generator<std::mt19937>);
}

/* Figure 8.3: Average gap for various processes for b = {n, 2n, ... , 50n}
   at m = n^2 balls with weights from Exp(1). */
void generate_several_processes_vs_batch_size_exp_weights() {
  generate_several_processes_vs_batch_size<double, std::mt19937>(exp_weights_generator<std::mt19937>);
}

template<typename Generator = std::mt19937>
void process_with_small_p_max_vs_batch_size() {
  int num_bins = 1'000;
  int repetitions = 50;

  std::cout << "(1+beta) with beta=0.5:" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 225; normalised_batch_size += 5) {
    int num_rounds = num_bins / normalised_batch_size;
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<size_t, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions, 
        one_plus_beta_with_random_tie_breaking<size_t, Generator>(0.5)) << ")" << std::endl;
  }

  std::cout << "(1+beta) with beta=" << 1 / std::log(num_bins) << ":" << std::endl;
  for (int normalised_batch_size = 1; normalised_batch_size <= 225; normalised_batch_size += 5) {
    int num_rounds = num_bins / normalised_batch_size;
    double beta = 1 / std::log(num_bins);
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<size_t, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions,
        one_plus_beta_with_random_tie_breaking<size_t, Generator>(beta)) << ")" << std::endl;
  }

  std::cout << "Quantile(delta) with delta=" << 1 / std::log(num_bins) << ":" << std::endl;
  for (int normalised_batch_size = 210; normalised_batch_size <= 225; normalised_batch_size += 5) {
    int num_rounds = num_bins / normalised_batch_size;
    double delta = 1 / std::log(num_bins);
    std::cout << "(" << normalised_batch_size << ", " <<
      avg_gap<size_t, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions,
        hybrid_quantile<size_t, Generator>(delta)) << ")" << std::endl;
  }
}

/* Returns the average gap for every batch of the process. */
template<typename LoadType = size_t, typename Generator = std::mt19937>
void avg_gap_every_batch(
  int num_bins,
  int normalised_batch_size,
  int num_batches,
  int repetitions,
  const std::function<size_t(const std::vector<LoadType>&, Generator&)> process,
  const std::function<LoadType(Generator&)> weights_generator = unit_weights_generator<Generator>) {

  Generator generator;
  double gap_total = 0.0;
  std::unordered_map<int, LoadType> avg_gap;
  for (int i = 0; i < repetitions; ++i) {
    BatchedSetting<LoadType, Generator> batched_setting(
      num_bins, normalised_batch_size * num_bins, process, weights_generator);
    LoadType cur_gap = 0.0;
    for (int round = 0; round < num_batches; ++round) {
      batched_setting.nextRound(generator);
      cur_gap = std::max(cur_gap, batched_setting.getGap());
      avg_gap[round] += cur_gap;
    }
  }
  for (int round = 0; round < num_batches; ++round) {
    std::cout << "(" << round << ", " << avg_gap[round] / double(repetitions) << ")" << std::endl;
  }
}

/* Figure 8.4: Average gap for Two-Choice with and without tie-breaking. */
template<typename Generator = std::mt19937>
void two_choice_with_and_without_random_tie_breaking() {
  int num_bins = 1'000;
  int normalised_batch_size = 25;
  int num_rounds = 40;
  int repetitions = 100;

  avg_gap_every_batch<size_t, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions,
    d_choice<size_t, Generator, 2>);

  avg_gap_every_batch<size_t, Generator>(num_bins, normalised_batch_size, num_rounds, repetitions,
    two_choice_with_random_tie_breaking<size_t, Generator>);
}

int main() {
  std::cout << "=== Figure 8.1 ===" << std::endl;
  generate_several_processes_vs_batch_size_unit_weights();
  std::cout << "=== Figure 8.2 ===" << std::endl;
  process_with_small_p_max_vs_batch_size();
  std::cout << "=== Figure 8.3 ===" << std::endl;
  generate_several_processes_vs_batch_size_exp_weights();
  std::cout << "=== Figure 8.4 ===" << std::endl;
  two_choice_with_and_without_random_tie_breaking();

  return 0;
}
