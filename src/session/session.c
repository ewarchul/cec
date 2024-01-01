#include "cec/session.h"

#include "benchmark_info.h"

cec_session* create_cec_session(cec_version version, byte* suite) {
  cec_session*   s    = generic_alloc(cec_session, 1000);
  benchmark_info info = mk_benchmark_info(version);
  s->version          = version;
  s->suite            = suite;
  s->allocator        = create_arena(120000);
  s->ctx              = generic_alloc(cec_ctx*, info.dim_num);
  for (u8 i = 0; i < info.dim_num; ++i) {
    s->ctx[info.dimensions[i]] =
        mk_ctx(info.storage_dir, info.dimensions[i], info.problem_num);
  }
  free(info.dimensions);
  return s;
}

void destroy_cec_session(cec_session* session) {
  benchmark_info info = mk_benchmark_info(session->version);
  for (u8 i = 0; i < info.dim_num; ++i) {
    free_cec_ctx(session->ctx[info.dimensions[i]]);
  }
  free_arena(&session->allocator);
  generic_free(info.dimensions);
  generic_free(session->ctx);
  generic_free(session);
}
