/**
 * Parallel VLSI Wire Routing via OpenMP
 * Arjun Walia(arjunw)
 */

#include "wireroute.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <vector>
#include <climits>
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

inline void set_wire_route_by_index(Wire& w, int route_idx) {
  int dx = std::abs(w.end_x - w.start_x);
  int dy = std::abs(w.end_y - w.start_y);

  if (dx == 0 || dy == 0) {
    w.bend1_x = w.start_x;
    w.bend1_y = w.start_y;
    return;
  }

  if (route_idx < dx) {
    int offset = route_idx + 1;
    w.bend1_x = (w.end_x > w.start_x) ? w.start_x + offset : w.start_x - offset;
    w.bend1_y = w.start_y;
  } else {
    int offset = route_idx - dx + 1;
    w.bend1_x = w.start_x;
    w.bend1_y = (w.end_y > w.start_y) ? w.start_y + offset : w.start_y - offset;
  }
}

void randomize_wire_route(Wire &w) {
  int dx = std::abs(w.end_x - w.start_x);
  int dy = std::abs(w.end_y - w.start_y);

  if (dx == 0 || dy == 0) {
    w.bend1_x = w.start_x;
    w.bend1_y = w.start_y;
    return;
  }

  thread_local std::mt19937 device(std::random_device{}());
  std::uniform_int_distribution<int> route_distance(0, dx + dy - 1);
  set_wire_route_by_index(w, route_distance(device));
}

struct Point {
    int x, y;
};

std::vector<Point> get_wire_path(const Wire& w) {
  std::vector<Point> path;
  int x = w.start_x, y = w.start_y;

  bool vertical_first = (w.bend1_x == w.start_x);

  if (vertical_first) {
    int dir = (w.bend1_y > y) ? 1 : -1;

    while (y != w.bend1_y) {
      path.push_back({x, y});
      y += dir;
    }

    dir = (w.end_x > x) ? 1 : -1;

    while (x != w.end_x) {
      path.push_back({x, y});
      x += dir;
    }

    dir = (w.end_y > y) ? 1 : -1;

    while (y != w.end_y) {
      path.push_back({x, y});
      y += dir;
    }
  } else {
    int dir = (w.bend1_x > x) ? 1 : -1;
    
    while (x != w.bend1_x) {
      path.push_back({x, y});
      x += dir;
    }
    
    dir = (w.end_y > y) ? 1 : -1;
    
    while (y != w.end_y) {
      path.push_back({x, y});
      y += dir;
    }
    
    dir = (w.end_x > x) ? 1 : -1;
    
    while (x != w.end_x) {
      path.push_back({x, y});
      x += dir;
    }
  }
    
  path.push_back({x, y});
  return path;
}

long long calculate_path_cost(const std::vector<std::vector<int>>& grid, const Wire& w, int test_bend_x, int test_bend_y) {
  long long cost = 0;
  int x = w.start_x, y = w.start_y;

  bool vertical_first = (test_bend_x == w.start_x);

  auto add_cost = [&]() {
      int occ = grid[y][x];
      cost += (long long)(occ + 1) * (occ + 1);
  };
    
  if (vertical_first) {
    int dir = (test_bend_y > y) ? 1 : -1;

    while (y != test_bend_y) {
      add_cost();
      y += dir;
    }

    dir = (w.end_x > x) ? 1 : -1;

    while (x != w.end_x) {
      add_cost();
      x += dir;
    }

    dir = (w.end_y > y) ? 1 : -1;

    while (y != w.end_y) {
      add_cost();
      y += dir;
    }
  } else {
    int dir = (test_bend_x > x) ? 1 : -1;

    while (x != test_bend_x) {
      add_cost();
      x += dir;
    }

    dir = (w.end_y > y) ? 1 : -1;

    while (y != w.end_y) {
      add_cost();
      y += dir;
    }

    dir = (w.end_x > x) ? 1 : -1;

    while (x != w.end_x) {
      add_cost();
      x += dir;
    }
  }

  add_cost();
  return cost;
}

void optimize_single_wire_route(Wire& w, const std::vector<std::vector<int>>& grid) {
  int dx = std::abs(w.end_x - w.start_x);
  int dy = std::abs(w.end_y - w.start_y);

  if (dx == 0 || dy == 0) {
    w.bend1_x = w.start_x;
    w.bend1_y = w.start_y;
    return;
  }

  long long best_cost = LLONG_MAX;
  int best_route = 0;
  int total_routes = dx + dy;

  for (int route = 0; route < total_routes; ++route) {
    int bend_x, bend_y;

    if (route < dx) {
      int offset = route + 1;
      bend_x = (w.end_x > w.start_x) ? w.start_x + offset : w.start_x - offset;
      bend_y = w.start_y;
    } else {
      int offset = route - dx + 1;
      bend_x = w.start_x;
      bend_y = (w.end_y > w.start_y) ? w.start_y + offset : w.start_y - offset;
    }

    long long cost = calculate_path_cost(grid, w, bend_x, bend_y);

    if (cost < best_cost) {
      best_cost = cost;
      best_route = route;
    }
  }

  set_wire_route_by_index(w, best_route);
}

void route_within_wires(std::vector<std::vector<int>>& grid, std::vector<std::vector<int>>& grid_transposed, std::vector<Wire>& wires, double sa_prob, int num_iters) {
  std::random_device device;
  std::mt19937 generate_engine(device());
  std::bernoulli_distribution probability_check(1.0 - sa_prob);
    
  auto update_dual_grid = [](std::vector<std::vector<int>>& g_yx, std::vector<std::vector<int>>& g_xy, const Wire& w, int delta) {
    int x = w.start_x, y = w.start_y;
    bool vert_first = (w.bend1_x == w.start_x);
        
    if (vert_first) {
      int dir = (w.bend1_y > y) ? 1 : -1;

      while (y != w.bend1_y) {
        g_yx[y][x] += delta;
        g_xy[x][y] += delta;
        y += dir;
      }

      dir = (w.end_x > x) ? 1 : -1;

      while (x != w.end_x) {
        g_yx[y][x] += delta;
        g_xy[x][y] += delta;
        x += dir;
      }

      dir = (w.end_y > y) ? 1 : -1;

      while (y != w.end_y) {
        g_yx[y][x] += delta;
        g_xy[x][y] += delta;
        y += dir;
      }

      g_yx[y][x] += delta;
      g_xy[x][y] += delta;
    } else {
      int dir = (w.bend1_x > x) ? 1 : -1;

      while (x != w.bend1_x) {
        g_yx[y][x] += delta;
        g_xy[x][y] += delta;
        x += dir;
      }

      dir = (w.end_y > y) ? 1 : -1;

      while (y != w.end_y) {
        g_yx[y][x] += delta;
        g_xy[x][y] += delta;
        y += dir;
      }

      dir = (w.end_x > x) ? 1 : -1;
      while (x != w.end_x) {
        g_yx[y][x] += delta;
        g_xy[x][y] += delta;
        x += dir;
      }

      g_yx[y][x] += delta;
      g_xy[x][y] += delta;
    }
  };
    
  auto calc_cost_opt = [&grid, &grid_transposed](const Wire& w, int bx, int by) -> long long {
    long long cost = 0;
    int x = w.start_x, y = w.start_y;
    bool vert = (bx == w.start_x);

    auto sq = [](int v) { return (long long)v * v; };

    if (vert) {
      int dir = (by > y) ? 1 : -1;

      while (y != by) {
        cost += sq(grid_transposed[x][y] + 1);
        y += dir;
      }

      dir = (w.end_x > x) ? 1 : -1;

      while (x != w.end_x) {
        cost += sq(grid[y][x] + 1);
        x += dir;
      }

      dir = (w.end_y > y) ? 1 : -1;

      while (y != w.end_y) {
        cost += sq(grid_transposed[x][y] + 1);
        y += dir;
      }

      cost += sq(grid_transposed[x][y] + 1);
    } else {
      int dir = (bx > x) ? 1 : -1;

      while (x != bx) {
        cost += sq(grid[y][x] + 1);
        x += dir;
      }

      dir = (w.end_y > y) ? 1 : -1;

      while (y != w.end_y) {
        cost += sq(grid_transposed[x][y] + 1);
        y += dir;
      }

      dir = (w.end_x > x) ? 1 : -1;

      while (x != w.end_x) {
        cost += sq(grid[y][x] + 1);
        x += dir;
      }

      cost += sq(grid[y][x] + 1);
    }

    return cost;
    };
    
  for (auto& w : wires) {
    int dx = std::abs(w.end_x - w.start_x);
    int dy = std::abs(w.end_y - w.start_y);
        
    if (dx == 0 || dy == 0) {
      w.bend1_x = w.start_x;
      w.bend1_y = w.start_y;
    } else {
      long long best_cost = LLONG_MAX;
      int best_route = 0;
      int total_routes = dx + dy;

      for (int route = 0; route < total_routes; ++route) {
        int bend_x, bend_y;

        if (route < dx) {
          int offset = route + 1;
          bend_x = (w.end_x > w.start_x) ? w.start_x + offset : w.start_x - offset;
          bend_y = w.start_y;
        } else {
          int offset = route - dx + 1;
          bend_x = w.start_x;
          bend_y = (w.end_y > w.start_y) ? w.start_y + offset : w.start_y - offset;
        }

        long long cost = calc_cost_opt(w, bend_x, bend_y);
        if (cost < best_cost) {
          best_cost = cost;
          best_route = route;
        }
      }

      set_wire_route_by_index(w, best_route);
    }

    update_dual_grid(grid, grid_transposed, w, 1);
  }

  for (int iter = 1; iter < num_iters; ++iter) {
    for (auto& w : wires) {
      update_dual_grid(grid, grid_transposed, w, -1);

      if (probability_check(generate_engine)) {
        int dx = std::abs(w.end_x - w.start_x);
        int dy = std::abs(w.end_y - w.start_y);

        if (dx == 0 || dy == 0) {
          w.bend1_x = w.start_x;
          w.bend1_y = w.start_y;
        } else {
          int total_routes = dx + dy;
          long long global_best_cost = LLONG_MAX;
          int global_best_route = 0;

          #pragma omp parallel
          {
            long long local_best_cost = LLONG_MAX;
            int local_best_route = 0;
                        
            #pragma omp for schedule(static)
            for (int route = 0; route < total_routes; ++route) {
              int bend_x, bend_y;

              if (route < dx) {
                int offset = route + 1;
                bend_x = (w.end_x > w.start_x) ? w.start_x + offset : w.start_x - offset;
                bend_y = w.start_y;
              } else {
                int offset = route - dx + 1;
                bend_x = w.start_x;
                bend_y = (w.end_y > w.start_y) ? w.start_y + offset : w.start_y - offset;
              }

              long long cost = calc_cost_opt(w, bend_x, bend_y);
              
              if (cost < local_best_cost) {
                local_best_cost = cost;
                local_best_route = route;
              }
            }

            #pragma omp critical
            {
              if (local_best_cost < global_best_cost) {
                global_best_cost = local_best_cost;
                global_best_route = local_best_route;
              }
            }
          }

          set_wire_route_by_index(w, global_best_route);
        }
      } else {
        randomize_wire_route(w);
      }

      update_dual_grid(grid, grid_transposed, w, 1);
    }
  }
}

void route_across_wires(std::vector<std::vector<int>>& grid, std::vector<Wire>& wires, double sa_prob, int num_iters, int batch_size) {
  std::random_device rd;
  std::mt19937 generate_engine(rd());
  std::bernoulli_distribution probability_check(1.0 - sa_prob);

  size_t num_wires = wires.size();
  int grid_width = grid[0].size();

  std::vector<omp_lock_t> col_locks(grid_width);
  for (auto& lock : col_locks) {
    omp_init_lock(&lock);
  }
    
  for (int iter = 0; iter < num_iters; ++iter) {
    #pragma omp parallel
    {
      #pragma omp for schedule(dynamic, batch_size)

      for (size_t i = 0; i < num_wires; ++i) {
        Wire& w = wires[i];

        if (iter == 0) {
          randomize_wire_route(w);

          auto path = get_wire_path(w);
          std::vector<int> affected_cols;

          for (const auto& p : path) {
            if (affected_cols.empty() || affected_cols.back() != p.x) {
              affected_cols.push_back(p.x);
            }
          }

          std::sort(affected_cols.begin(), affected_cols.end());

          for (int col : affected_cols) {
            omp_set_lock(&col_locks[col]);
          }

          for (const auto& p : path) {
            grid[p.y][p.x] += 1;
          }

          for (int col : affected_cols) {
            omp_unset_lock(&col_locks[col]);
          }
        } else {
          auto old_path = get_wire_path(w);
          std::vector<int> old_cols;

          for (const auto& p : old_path) {
            if (old_cols.empty() || old_cols.back() != p.x) {
              old_cols.push_back(p.x);
            }
          }

          std::sort(old_cols.begin(), old_cols.end());

          for (int col : old_cols) {
            omp_set_lock(&col_locks[col]);
          }

          for (const auto& p : old_path) {
            grid[p.y][p.x] -= 1;
          }

          for (int col : old_cols) {
            omp_unset_lock(&col_locks[col]);
          }

          if (probability_check(generate_engine)) {
            optimize_single_wire_route(w, grid);
          } else {
            randomize_wire_route(w);
          }

          auto new_path = get_wire_path(w);
          std::vector<int> new_cols;

          for (const auto& p : new_path) {
            if (new_cols.empty() || new_cols.back() != p.x) {
              new_cols.push_back(p.x);
            }
          }

          std::sort(new_cols.begin(), new_cols.end());

          for (int col : new_cols) {
            omp_set_lock(&col_locks[col]);
          }

          for (const auto& p : new_path) {
            grid[p.y][p.x] += 1;
          }

          for (int col : new_cols) {
            omp_unset_lock(&col_locks[col]);
          }
        }
      }
    }
  }

  for (auto& lock : col_locks) {
    omp_destroy_lock(&lock);
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
  std::vector<std::vector<int>> occupancy(dim_y, std::vector<int>(dim_x));
  std::vector<std::vector<int>> occupancy_transposed(dim_x, std::vector<int>(dim_y));

  for (auto& wire : wires) {
    fin >> wire.start_x >> wire.start_y >> wire.end_x >> wire.end_y;
    wire.bend1_x = wire.start_x;
    wire.bend1_y = wire.start_y;
  }

  /* Initialize any additional data structures needed in the algorithm */

  const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
  std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

    omp_set_num_threads(num_threads);

    if (parallel_mode == 'W') {
        route_within_wires(occupancy, occupancy_transposed, wires, SA_prob, SA_iters);
    } else if (parallel_mode == 'A') {
        route_across_wires(occupancy, wires, SA_prob, SA_iters, batch_size);
    }
  /** 
   * Implement the wire routing algorithm here
   * Feel free to structure the algorithm into different functions
   * Don't use global variables.
   * Use OpenMP to parallelize the algorithm. 
   */

  const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
  std::cout << "Computation time (sec): " << compute_time << '\n';

  /* Write wires and occupancy matrix to files */

  print_stats(occupancy);
  write_output(wires, num_wires, occupancy, dim_x, dim_y, num_threads, input_filename);
}

validate_wire_t Wire::to_validate_format(void) const {
  /* TODO(student): Implement this if you want to use the wr_checker. */
  /* See wireroute.h for details on validate_wire_t. */
  throw std::logic_error("to_validate_format not implemented.");
}