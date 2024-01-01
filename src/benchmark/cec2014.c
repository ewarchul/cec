#include "cec2014.h"

#include "problem/basic.h"
#include "cec/allocator.h"
#include "geometry/trans.h"

PRIVATE f64 const SCORE[30] = {-1400};

PRIVATE optim_problem PROBLEM[30] = {sphere_fn};

f64* cec2014(cec_session* session, u8 fn_num, mdspan_f64 input) {
  f64* output = arena_alloc(&session->allocator, f64, input.ncol);
  for (u8 row = 0; row < input.nrow; ++row) {
    span_f64 input_s     = {input.ncol, &input.data[row]};
    f64*     transformed = scale_shift_rotate(
        input_s, &session->ctx[input_s.size][fn_num - 1], &session->allocator);
    output[row] = PROBLEM[fn_num - 1](SPAN_F64{input.ncol, transformed});
    output[row] += SCORE[fn_num - 1];
  }
  return output;
}

// case 1:
// ellips_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
// case 2:
// bent_cigar_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
// case 3:
// discus_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1);
// case 4:
// rosenbrock_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1); // scale_mul
// -> 2.048 / 100.0 case 5: ackley_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M,
// 1, 1); case 6: weierstrass_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1,
// 1); // scale_mul -> 0.5 / 100.0 case 7: griewank_func(&x[i * nx], &f[i], nx,
// cd.OShift, cd.M, 1, 1); // scale_mul -> 600.0 / 100.0 case 8:
// rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 0); // scale_mul
// -> 5.12 / 100.0 case 9: rastrigin_func(&x[i * nx], &f[i], nx, cd.OShift,
// cd.M, 1, 1); case 10: schwefel_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M,
// 1, 0); // scale_mul -> 1000.0 / 100.0 case 11: schwefel_func(&x[i * nx],
// &f[i], nx, cd.OShift, cd.M, 1, 1); // scale_mul -> 1000.0 / 100.0 case 12:
// katsuura_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1); // scale_mul
// -> 5.0 / 100.0, case 13: happycat_func(&x[i * nx], &f[i], nx, cd.OShift,
// cd.M, 1, 1); // scale mul -> 5.0 / 100.0 case 14: hgbat_func(&x[i * nx],
// &f[i], nx, cd.OShift, cd.M, 1, 1); // sclae_mul -> 5.0 / 100.0, case 15:
// grie_rosen_func(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1, 1); // scale_mul
// -> 5.0 / 100.0 f[i] += 1500.0; break; case 16: escaffer6_func(&x[i * nx],
// &f[i], nx, cd.OShift, cd.M, 1, 1); case 17: cec2014_hf01(&x[i * nx], &f[i],
// nx, cd.OShift, cd.M, cd.SS, 1, 1); case 18: cec2014_hf02(&x[i * nx], &f[i],
// nx, cd.OShift, cd.M, cd.SS, 1, 1); case 19: cec2014_hf03(&x[i * nx], &f[i],
// nx, cd.OShift, cd.M, cd.SS, 1, 1); case 20: cec2014_hf04(&x[i * nx], &f[i],
// nx, cd.OShift, cd.M, cd.SS, 1, 1); case 21: cec2014_hf05(&x[i * nx], &f[i],
// nx, cd.OShift, cd.M, cd.SS, 1, 1); case 22: cec2014_hf06(&x[i * nx], &f[i],
// nx, cd.OShift, cd.M, cd.SS, 1, 1); case 23: cec2014_cf01(&x[i * nx], &f[i],
// nx, cd.OShift, cd.M, 1); case 24: cec2014_cf02(&x[i * nx], &f[i], nx,
// cd.OShift, cd.M, 1); case 25: cec2014_cf03(&x[i * nx], &f[i], nx, cd.OShift,
// cd.M, 1); case 26: cec2014_cf04(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
// case 27:
// cec2014_cf05(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
// case 28:
// cec2014_cf06(&x[i * nx], &f[i], nx, cd.OShift, cd.M, 1);
// case 29:
// cec2014_cf07(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1);
// case 30:
// cec2014_cf08(&x[i * nx], &f[i], nx, cd.OShift, cd.M, cd.SS, 1);
