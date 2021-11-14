#ifndef io_h
#define io_h

#include "types.h"
#include "utils.h"
#include "log.h"
#include <stdio.h>
#include <stdlib.h>

cec_external_data_t cec_read_external_data(cec_version_t version, int dim,
                                           int problem, cec_affineT_type_t tag);
cec_problem_data_t cec_load_problems_data(cec_version_t version,
                                          int problem_nums, int dim,
                                          cec_affineT_type_t tag);
cec_problem_data_t *cec_load_affine_data(cec_version_t version,
                                         cec_benchmark_info_t info,
                                         cec_affineT_type_t tag);
cec_benchmark_data_t cec_load_benchmark_data(cec_version_t version);

#endif
