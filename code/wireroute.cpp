/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#include "wireroute.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>
#include <cmath>
#include <random>

#include <unistd.h>
#include <omp.h>

void print_stats(const std::vector<std::vector<int>>& occupancy) {
  int max_occupancy = 0;
  long long total_cost = 0;

  for (const auto& row : occupancy) {
    for (const int count : row) {
      max_occupancy = std::max(max_occupancy, count);
      total_cost += count * count;
    }
  }

  std::cout << "Max occupancy: " << max_occupancy << '\n';
  std::cout << "Total cost: " << total_cost << '\n';
}

void write_output(const std::vector<Wire>& wires, const int num_wires, const std::vector<std::vector<int>>& occupancy, const int dim_x, const int dim_y, const int num_threads, std::string input_filename) {
  if (std::size(input_filename) >= 4 && input_filename.substr(std::size(input_filename) - 4) == ".txt") {
    input_filename.resize(std::size(input_filename) - 4);
  }

  const std::string occupancy_filename = input_filename + "_occupancy_" + std::to_string(num_threads) + ".txt";
  const std::string wires_filename = input_filename + "_wires_" + std::to_string(num_threads) + ".txt";

  std::ofstream out_occupancy(occupancy_filename, std::fstream::out);
  if (!out_occupancy) {
    std::cerr << "Unable to open file: " << occupancy_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_occupancy << dim_x << ' ' << dim_y << '\n';
  for (const auto& row : occupancy) {
    for (const int count : row) {
      out_occupancy << count << ' ';
    }
    out_occupancy << '\n';
  }

  out_occupancy.close();

  std::ofstream out_wires(wires_filename, std::fstream:: out);
  if (!out_wires) {
    std::cerr << "Unable to open file: " << wires_filename << '\n';
    exit(EXIT_FAILURE);
  }

  out_wires << dim_x << ' ' << dim_y << '\n' << num_wires << '\n';

  for (const auto& [start_x, start_y, end_x, end_y, bend1_x, bend1_y] : wires) {
    out_wires << start_x << ' ' << start_y << ' ' << bend1_x << ' ' << bend1_y << ' ';

    if (start_y == bend1_y) {
    // first bend was horizontal

      if (end_x != bend1_x) {
        // two bends

        out_wires << bend1_x << ' ' << end_y << ' ';
      }
    } else if (start_x == bend1_x) {
      // first bend was vertical

      if(end_y != bend1_y) {
        // two bends

        out_wires << end_x << ' ' << bend1_y << ' ';
      }
    }
    out_wires << end_x << ' ' << end_y << '\n';
  }

  out_wires.close();
}

long long calculate_wire_cost(const Wire& wire, const std::vector<std::vector<int>>& occupancy) {
  long long cost = 0;
  int x1 = wire.start_x, y1 = wire.start_y;
  int x2 = wire.bend1_x, y2 = wire.bend1_y;
  int x3 = wire.end_x, y3 = wire.end_y;
  
  if (x1 == x2) {
    int min_y = std::min(y1, y2);
    int max_y = std::max(y1, y2);
    for (int y = min_y; y <= max_y; y++) {
      int occ = occupancy[y][x1];
      cost += occ * occ;
    }
  } else {
    int min_x = std::min(x1, x2);
    int max_x = std::max(x1, x2);
    for (int x = min_x; x <= max_x; x++) {
      int occ = occupancy[y1][x];
      cost += occ * occ;
    }
  }
  
  if (x2 == x3) {
    int min_y = std::min(y2, y3);
    int max_y = std::max(y2, y3);
    for (int y = min_y; y <= max_y; y++) {
      if (y != y2) {
        int occ = occupancy[y][x2];
        cost += occ * occ;
      }
    }
  } else {
    int min_x = std::min(x2, x3);
    int max_x = std::max(x2, x3);
    for (int x = min_x; x <= max_x; x++) {
      if (x != x2) {
        int occ = occupancy[y2][x];
        cost += occ * occ;
      }
    }
  }
  
  return cost;
}

void update_occupancy(const Wire& wire, std::vector<std::vector<int>>& occupancy, int delta) {
  int x1 = wire.start_x, y1 = wire.start_y;
  int x2 = wire.bend1_x, y2 = wire.bend1_y;
  int x3 = wire.end_x, y3 = wire.end_y;
  
  if (x1 == x2) {
    int min_y = std::min(y1, y2);
    int max_y = std::max(y1, y2);
    for (int y = min_y; y <= max_y; y++) {
      occupancy[y][x1] += delta;
    }
  } else {
    int min_x = std::min(x1, x2);
    int max_x = std::max(x1, x2);
    for (int x = min_x; x <= max_x; x++) {
      occupancy[y1][x] += delta;
    }
  }
  
  if (x2 == x3) {
    int min_y = std::min(y2, y3);
    int max_y = std::max(y2, y3);
    for (int y = min_y; y <= max_y; y++) {
      if (y != y2) {
        occupancy[y][x2] += delta;
      }
    }
  } else {
    int min_x = std::min(x2, x3);
    int max_x = std::max(x2, x3);
    for (int x = min_x; x <= max_x; x++) {
      if (x != x2) {
        occupancy[y2][x] += delta;
      }
    }
  }
}

void route_within_wires(std::vector<Wire>& wires, std::vector<std::vector<int>>& occupancy, 
                        int num_threads, double SA_prob, int SA_iters) {
  const int num_wires = wires.size();
  
  for (int i = 0; i < num_wires; i++) {
    Wire& wire = wires[i];
    wire.bend1_x = wire.start_x;
    wire.bend1_y = wire.start_y;
    update_occupancy(wire, occupancy, 1);
  }
  
  for (int iter = 0; iter < SA_iters; iter++) {
    for (int i = 0; i < num_wires; i++) {
      Wire& wire = wires[i];
      
      if ((wire.start_x == wire.end_x) || (wire.start_y == wire.end_y)) {
        continue;
      }
      
      update_occupancy(wire, occupancy, -1);
      
      int delta_x = std::abs(wire.end_x - wire.start_x);
      int delta_y = std::abs(wire.end_y - wire.start_y);
      int num_routes = delta_x + delta_y;
      
      long long best_cost = LLONG_MAX;
      int best_bend_x = wire.bend1_x;
      int best_bend_y = wire.bend1_y;
      
      #pragma omp parallel num_threads(num_threads)
      {
        long long local_best_cost = LLONG_MAX;
        int local_best_bend_x = wire.bend1_x;
        int local_best_bend_y = wire.bend1_y;
        
        #pragma omp for schedule(static) nowait
        for (int dx = 0; dx <= delta_x; dx++) {
          Wire test_wire = wire;
          test_wire.bend1_x = wire.start_x + (wire.end_x > wire.start_x ? dx : -dx);
          test_wire.bend1_y = wire.end_y;
          
          long long cost = calculate_wire_cost(test_wire, occupancy);
          if (cost < local_best_cost) {
            local_best_cost = cost;
            local_best_bend_x = test_wire.bend1_x;
            local_best_bend_y = test_wire.bend1_y;
          }
        }
        
        #pragma omp for schedule(static) nowait
        for (int dy = 0; dy <= delta_y; dy++) {
          Wire test_wire = wire;
          test_wire.bend1_x = wire.end_x;
          test_wire.bend1_y = wire.start_y + (wire.end_y > wire.start_y ? dy : -dy);
          
          long long cost = calculate_wire_cost(test_wire, occupancy);
          if (cost < local_best_cost) {
            local_best_cost = cost;
            local_best_bend_x = test_wire.bend1_x;
            local_best_bend_y = test_wire.bend1_y;
          }
        }
        
        #pragma omp critical
        {
          if (local_best_cost < best_cost) {
            best_cost = local_best_cost;
            best_bend_x = local_best_bend_x;
            best_bend_y = local_best_bend_y;
          }
        }
      }
      
      thread_local std::mt19937 gen(omp_get_thread_num());
      std::uniform_real_distribution<> prob_dist(0.0, 1.0);
      std::uniform_int_distribution<> route_dist(0, num_routes - 1);
      
      if (prob_dist(gen) < SA_prob) {
        int random_route = route_dist(gen);
        if (random_route < delta_x) {
          wire.bend1_x = wire.start_x + (wire.end_x > wire.start_x ? random_route : -random_route);
          wire.bend1_y = wire.end_y;
        } else {
          int dy = random_route - delta_x;
          wire.bend1_x = wire.end_x;
          wire.bend1_y = wire.start_y + (wire.end_y > wire.start_y ? dy : -dy);
        }
      } else {
        wire.bend1_x = best_bend_x;
        wire.bend1_y = best_bend_y;
      }
      
      update_occupancy(wire, occupancy, 1);
    }
  }
}

void route_across_wires(std::vector<Wire>& wires, std::vector<std::vector<int>>& occupancy,
                        int num_threads, double SA_prob, int SA_iters, int batch_size) {
  const int num_wires = wires.size();
  const int dim_y = occupancy.size();
  const int dim_x = occupancy[0].size();
  
  std::vector<omp_lock_t> row_locks(dim_y);
  for (int i = 0; i < dim_y; i++) {
    omp_init_lock(&row_locks[i]);
  }
  
  #pragma omp parallel for schedule(dynamic, batch_size) num_threads(num_threads)
  for (int i = 0; i < num_wires; i++) {
    Wire& wire = wires[i];
    wire.bend1_x = wire.start_x;
    wire.bend1_y = wire.start_y;
    
    int x1 = wire.start_x, y1 = wire.start_y;
    int x2 = wire.bend1_x, y2 = wire.bend1_y;
    int x3 = wire.end_x, y3 = wire.end_y;
    
    std::vector<int> rows_to_lock;
    if (x1 == x2) {
      int min_y = std::min(y1, y2);
      int max_y = std::max(y1, y2);
      for (int y = min_y; y <= max_y; y++) rows_to_lock.push_back(y);
    } else {
      rows_to_lock.push_back(y1);
    }
    
    if (x2 == x3) {
      int min_y = std::min(y2, y3);
      int max_y = std::max(y2, y3);
      for (int y = min_y; y <= max_y; y++) {
        if (y != y2 || x1 != x2) rows_to_lock.push_back(y);
      }
    } else {
      if (y2 != y1 || x1 == x2) rows_to_lock.push_back(y2);
    }
    
    std::sort(rows_to_lock.begin(), rows_to_lock.end());
    rows_to_lock.erase(std::unique(rows_to_lock.begin(), rows_to_lock.end()), rows_to_lock.end());
    
    for (int row : rows_to_lock) {
      omp_set_lock(&row_locks[row]);
    }
    
    if (x1 == x2) {
      int min_y = std::min(y1, y2);
      int max_y = std::max(y1, y2);
      for (int y = min_y; y <= max_y; y++) {
        occupancy[y][x1]++;
      }
    } else {
      int min_x = std::min(x1, x2);
      int max_x = std::max(x1, x2);
      for (int x = min_x; x <= max_x; x++) {
        occupancy[y1][x]++;
      }
    }
    
    if (x2 == x3) {
      int min_y = std::min(y2, y3);
      int max_y = std::max(y2, y3);
      for (int y = min_y; y <= max_y; y++) {
        if (y != y2) {
          occupancy[y][x2]++;
        }
      }
    } else {
      int min_x = std::min(x2, x3);
      int max_x = std::max(x2, x3);
      for (int x = min_x; x <= max_x; x++) {
        if (x != x2) {
          occupancy[y2][x]++;
        }
      }
    }
    
    for (int row : rows_to_lock) {
      omp_unset_lock(&row_locks[row]);
    }
  }
  
  for (int iter = 0; iter < SA_iters; iter++) {
    #pragma omp parallel num_threads(num_threads)
    {
      int tid = omp_get_thread_num();
      std::mt19937 gen(tid * 1000 + iter);
      std::uniform_real_distribution<> prob_dist(0.0, 1.0);
      
      #pragma omp for schedule(dynamic, batch_size)
      for (int batch_start = 0; batch_start < num_wires; batch_start += batch_size) {
        int batch_end = std::min(batch_start + batch_size, num_wires);
        std::vector<Wire> new_routes(batch_end - batch_start);
        
        for (int idx = 0; idx < batch_end - batch_start; idx++) {
          int i = batch_start + idx;
          Wire& wire = wires[i];
          
          if ((wire.start_x == wire.end_x) || (wire.start_y == wire.end_y)) {
            new_routes[idx] = wire;
            continue;
          }
          
          int delta_x = std::abs(wire.end_x - wire.start_x);
          int delta_y = std::abs(wire.end_y - wire.start_y);
          int num_routes = delta_x + delta_y;
          
          long long best_cost = calculate_wire_cost(wire, occupancy);
          Wire best_wire = wire;
          
          for (int dx = 0; dx <= delta_x; dx++) {
            Wire test_wire = wire;
            test_wire.bend1_x = wire.start_x + (wire.end_x > wire.start_x ? dx : -dx);
            test_wire.bend1_y = wire.end_y;
            
            long long cost = calculate_wire_cost(test_wire, occupancy);
            if (cost < best_cost) {
              best_cost = cost;
              best_wire = test_wire;
            }
          }
          
          for (int dy = 0; dy <= delta_y; dy++) {
            Wire test_wire = wire;
            test_wire.bend1_x = wire.end_x;
            test_wire.bend1_y = wire.start_y + (wire.end_y > wire.start_y ? dy : -dy);
            
            long long cost = calculate_wire_cost(test_wire, occupancy);
            if (cost < best_cost) {
              best_cost = cost;
              best_wire = test_wire;
            }
          }
          
          std::uniform_int_distribution<> route_dist(0, num_routes - 1);
          if (prob_dist(gen) < SA_prob) {
            int random_route = route_dist(gen);
            if (random_route < delta_x) {
              best_wire.bend1_x = wire.start_x + (wire.end_x > wire.start_x ? random_route : -random_route);
              best_wire.bend1_y = wire.end_y;
            } else {
              int dy = random_route - delta_x;
              best_wire.bend1_x = wire.end_x;
              best_wire.bend1_y = wire.start_y + (wire.end_y > wire.start_y ? dy : -dy);
            }
          }
          
          new_routes[idx] = best_wire;
        }
        
        for (int idx = 0; idx < batch_end - batch_start; idx++) {
          int i = batch_start + idx;
          Wire& old_wire = wires[i];
          Wire& new_wire = new_routes[idx];
          
          if (old_wire.bend1_x == new_wire.bend1_x && old_wire.bend1_y == new_wire.bend1_y) {
            continue;
          }
          
          std::vector<int> rows_to_lock;
          
          int x1 = old_wire.start_x, y1 = old_wire.start_y;
          int x2 = old_wire.bend1_x, y2 = old_wire.bend1_y;
          int x3 = old_wire.end_x, y3 = old_wire.end_y;
          
          if (x1 == x2) {
            int min_y = std::min(y1, y2);
            int max_y = std::max(y1, y2);
            for (int y = min_y; y <= max_y; y++) rows_to_lock.push_back(y);
          } else {
            rows_to_lock.push_back(y1);
          }
          if (x2 == x3) {
            int min_y = std::min(y2, y3);
            int max_y = std::max(y2, y3);
            for (int y = min_y; y <= max_y; y++) rows_to_lock.push_back(y);
          } else {
            rows_to_lock.push_back(y2);
          }
          
          x2 = new_wire.bend1_x;
          y2 = new_wire.bend1_y;
          
          if (x1 == x2) {
            int min_y = std::min(y1, y2);
            int max_y = std::max(y1, y2);
            for (int y = min_y; y <= max_y; y++) rows_to_lock.push_back(y);
          } else {
            rows_to_lock.push_back(y1);
          }
          if (x2 == x3) {
            int min_y = std::min(y2, y3);
            int max_y = std::max(y2, y3);
            for (int y = min_y; y <= max_y; y++) rows_to_lock.push_back(y);
          } else {
            rows_to_lock.push_back(y2);
          }
          
          std::sort(rows_to_lock.begin(), rows_to_lock.end());
          rows_to_lock.erase(std::unique(rows_to_lock.begin(), rows_to_lock.end()), rows_to_lock.end());
          
          for (int row : rows_to_lock) {
            omp_set_lock(&row_locks[row]);
          }
          
          update_occupancy(old_wire, occupancy, -1);
          update_occupancy(new_wire, occupancy, 1);
          
          for (int row : rows_to_lock) {
            omp_unset_lock(&row_locks[row]);
          }
          
          wires[i] = new_wire;
        }
      }
    }
  }
  
  for (int i = 0; i < dim_y; i++) {
    omp_destroy_lock(&row_locks[i]);
  }
}

int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();

  std::string input_filename;
  int num_threads = 0;
  double SA_prob = 0.1;
  int SA_iters = 5;
  char parallel_mode = '\0';
  int batch_size = 1;

  int opt;
  while ((opt = getopt(argc, argv, "f:n:p:i:m:b:")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;
      case 'n':
        num_threads = atoi(optarg);
        break;
      case 'p':
        SA_prob = atof(optarg);
        break;
      case 'i':
        SA_iters = atoi(optarg);
        break;
      case 'm':
        parallel_mode = *optarg;
        break;
      case 'b':
        batch_size = atoi(optarg);
        break;
      default:
        std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
        exit(EXIT_FAILURE);
    }
  }

  // Check if required options are provided
  if (empty(input_filename) || num_threads <= 0 || SA_iters <= 0 || (parallel_mode != 'A' && parallel_mode != 'W') || batch_size <= 0) {
    std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_threads [-p SA_prob] [-i SA_iters] -m parallel_mode -b batch_size\n";
    exit(EXIT_FAILURE);
  }

  std::cout << "Number of threads: " << num_threads << '\n';
  std::cout << "Simulated annealing probability parameter: " << SA_prob << '\n';
  std::cout << "Simulated annealing iterations: " << SA_iters << '\n';
  std::cout << "Input file: " << input_filename << '\n';
  std::cout << "Parallel mode: " << parallel_mode << '\n';
  std::cout << "Batch size: " << batch_size << '\n';

  std::ifstream fin(input_filename);

  if (!fin) {
    std::cerr << "Unable to open file: " << input_filename << ".\n";
    exit(EXIT_FAILURE);
  }

  int dim_x, dim_y;
  int num_wires;

  /* Read the grid dimension and wire information from file */
  fin >> dim_x >> dim_y >> num_wires;

  std::vector<Wire> wires(num_wires);
  std::vector occupancy(dim_y, std::vector<int>(dim_x));

  for (auto& wire : wires) {
    fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
    wire.bend1_x = wire.start_x;
    wire.bend1_y = wire.start_y;
  }

  /* Initialize any additional data structures needed in the algorithm */
  omp_set_num_threads(num_threads);

  const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
  std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  /** 
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Don't use global variables.
   * Use OpenMP to parallelize the algorithm. 
   */
  
  if (parallel_mode == 'W') {
    route_within_wires(wires, occupancy, num_threads, SA_prob, SA_iters);
  } else if (parallel_mode == 'A') {
    route_across_wires(wires, occupancy, num_threads, SA_prob, SA_iters, batch_size);
  }

  const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
  std::cout << "Computation time (sec): " << compute_time << '\n';

  /* Write wires and occupancy matrix to files */

  print_stats(occupancy);
  write_output(wires, num_wires, occupancy, dim_x, dim_y, num_threads, input_filename);
}

validate_wire_t Wire::to_validate_format(void) const {
  /* TODO(student): Implement this if you want to use the wr_checker. */
  /* See wireroute.h for details on validate_wire_t. */
  validate_wire_t vwire;
  
  vwire.p[0].x = start_x;
  vwire.p[0].y = start_y;
  
  int point_count = 1;
  
  if (bend1_x != start_x || bend1_y != start_y) {
    vwire.p[point_count].x = bend1_x;
    vwire.p[point_count].y = bend1_y;
    point_count++;
  }
  
  if (start_y == bend1_y) {
    if (end_x != bend1_x) {
      vwire.p[point_count].x = bend1_x;
      vwire.p[point_count].y = end_y;
      point_count++;
    }
  } else if (start_x == bend1_x) {
    if (end_y != bend1_y) {
      vwire.p[point_count].x = end_x;
      vwire.p[point_count].y = bend1_y;
      point_count++;
    }
  }
  
  if (vwire.p[point_count-1].x != end_x || vwire.p[point_count-1].y != end_y) {
    vwire.p[point_count].x = end_x;
    vwire.p[point_count].y = end_y;
    point_count++;
  }
  
  vwire.num_pts = point_count;
  
  return vwire;
}