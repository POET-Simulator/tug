#include <bits/stdint-uintn.h>
#include <diffusion/BTCSDiffusion.hpp>
#include <doctest/doctest.h>
#include <grid/BTCSBoundaryCondition.hpp>
#include <vector>

using namespace tug::boundary_condition;
using namespace tug::diffusion;

#define DIMENSION 2
#define N 51
#define M 51
#define MID 1300

static std::vector<double> alpha(N *M, 1e-3);

static BTCSDiffusion setupDiffu(uint32_t n, uint32_t m) {
  BTCSDiffusion diffu(DIMENSION);

  diffu.setXDimensions(n, n);
  diffu.setYDimensions(m, m);

  diffu.setTimestep(1.);

  return diffu;
}

TEST_CASE("closed boundaries - 1 concentration to 1 - rest 0") {
  std::vector<double> field(N * M, 0);

  field[MID] = 1;

  BTCSDiffusion diffu = setupDiffu(N, M);
  BoundaryCondition bc(N, M);

  uint32_t iterations = 1000;
  double sum = 0;

  for (int t = 0; t < iterations; t++) {
    diffu.simulate(field.data(), alpha.data(), bc);

    if (t == iterations - 1) {
      // iterate through rows
      for (int i = 0; i < M; i++) {
        // iterate through columns
        for (int j = 0; j < N; j++) {
          sum += field[i * N + j];
        }
      }
    }
  }
  CAPTURE(sum);
  // epsilon of 1e-8
  CHECK(sum == doctest::Approx(1).epsilon(1e-6));
}

TEST_CASE("constant boundaries (0) - 1 concentration to 1 - rest 0") {
  std::vector<double> field(N * M, 0);

  field[MID] = 1;

  BTCSDiffusion diffu = setupDiffu(N, M);
  BoundaryCondition bc(N, M);

  boundary_condition input = {BC_TYPE_CONSTANT, 0};

  bc.setSide(BC_SIDE_LEFT, input);
  bc.setSide(BC_SIDE_RIGHT, input);
  bc.setSide(BC_SIDE_TOP, input);
  bc.setSide(BC_SIDE_BOTTOM, input);

  uint32_t max_iterations = 20000;
  bool reached = false;

  int t = 0;

  for (t = 0; t < max_iterations; t++) {
    diffu.simulate(field.data(), alpha.data(), bc);

    if (field[N * M - 1] > 1e-15) {
      reached = true;
      break;
    }
  }

  if (!reached) {
    CAPTURE(field[N * M - 1]);
    FAIL_CHECK(
        "Concentration did not reach boundaries after count of iterations: ",
        t);
  }
}

TEST_CASE(
    "constant top and bottom (1 and 0) - left and right closed - 0 inlet") {
  std::vector<double> field(N * M, 0);

  BTCSDiffusion diffu = setupDiffu(N, M);
  BoundaryCondition bc(N, M);

  boundary_condition top =
      BoundaryCondition::returnBoundaryCondition(BC_TYPE_CONSTANT, 1);
  boundary_condition bottom =
      BoundaryCondition::returnBoundaryCondition(BC_TYPE_CONSTANT, 0);

  bc.setSide(BC_SIDE_TOP, top);
  bc.setSide(BC_SIDE_BOTTOM, bottom);

  uint32_t max_iterations = 100;

  for (int t = 0; t < max_iterations; t++) {
    diffu.simulate(field.data(), alpha.data(), bc);
  }

  for (int i = 0; i < N; i++) {
    double above = field[i];
    for (int j = 1; j < M; j++) {
      double curr = field[j * N + i];
      if (curr > above) {
        CAPTURE(curr);
        CAPTURE(above);
        FAIL("Concentration below is greater than above @ cell ", j * N + i);
      }
    }
  }
}

TEST_CASE("2D closed boundaries, 1 constant cell in the middle") {
  std::vector<double> field(N * M, 0);
  double val = 1e-2;

  BTCSDiffusion diffu = setupDiffu(N, M);
  BoundaryCondition bc(N, M);

  field[MID] = val;
  bc(BC_INNER, MID) = {BC_TYPE_CONSTANT, val};

  uint32_t max_iterations = 100;

  double sum = val;

  for (int t = 0; t < max_iterations; t++) {
    diffu.simulate(field.data(), alpha.data(), bc);

    CHECK_EQ(field[MID], val);

    double new_sum = .0;

    for (int i = 0; i < M; i++) {
      // iterate through columns
      for (int j = 0; j < N; j++) {
        new_sum += field[i * N + j];
      }
    }

    if (sum > new_sum) {
      CAPTURE(sum);
      CAPTURE(new_sum);
      FAIL("Sum of concentrations is equal or lower than to previous iteration "
           "after iteration ",
           t + 1);
    }

    sum = new_sum;
  }
}
