#include "affine_trans.h"
#include "fixture.h"
#include "tau/tau.h"
#include "utils.h"

TEST(AffineTransforms, Shift) {
  double input[10] = {[0 ... 9] = 1};

  double output_old[10];
  shiftfunc(input, output_old, 10, fixture_shift_data);

  cec_state_t fixture_state = cec_mk_state(CEC_2017, 10, NONE);
  double *output_new = malloc(10 * sizeof(double));
  output_new = shift_modern(input, 1, &fixture_state);

  for (int i = 0; i < 10; ++i) {
    REQUIRE(output_new[i] == output_old[i], "Shifted vectors are not equal.");
  }

  free(output_new);
}

TEST(AffineTransforms, Rotate) {
  double input[10] = {[0 ... 9] = 1};

  double output_old[10];
  rotatefunc(input, output_old, 10, fixture_rot_data);

  cec_state_t fixture_state = cec_mk_state(CEC_2017, 10, NONE);
  double *output_new = malloc(10 * sizeof(double));
  output_new = rotate_modern(input, 1, &fixture_state);

  for (int i = 0; i < 10; ++i) {
    REQUIRE(output_new[i] == output_old[i], "Rotated vectors are not equal.");
  }

  free(output_new);
}

void sr_func(double *x, double *sr_x, int nx, double *Os, double *Mr,
             double sh_rate, int s_flag, int r_flag, double *y);

TEST(AffineTransforms, ShiftRotate) {
  double input[10] = {[0 ... 9] = 1};

  double output_old[10], output_old_temp[10];
  sr_func(input, output_old, 10, fixture_shift_data, fixture_rot_data, 1.0, 1,
          1, output_old_temp);

  cec_state_t fixture_state = cec_mk_state(CEC_2017, 10, NONE);
  double *output_new = malloc(10 * sizeof(double));
  cec_affine_transforms_t transforms = {
      .bias_ = false,
      .rotate_ = true,
      .shift_ = true,
      .shift_rotate_ = false,
      .shuffle_ = false,
      .transform_rate_ = 1,
  };
  output_new = shift_rotate_modern(input, 1, &fixture_state, transforms);

  for (int i = 0; i < 10; ++i) {
    REQUIRE(output_new[i] == output_old[i],
            "Shifted and rotated vectors are not equal.");
  }

  free(output_new);
}
