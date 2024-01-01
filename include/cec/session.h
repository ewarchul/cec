#pragma once

#include "types.h"
#include "context.h"
#include "allocator.h"

typedef struct cec_session cec_session;
struct cec_session {
  cec_version     version;
  byte*           suite;
  cec_ctx**       ctx;
  u8              dim_num;
  arena_allocator allocator;
};

cec_session* create_cec_session(cec_version version, byte* suite);
void         destroy_cec_session(cec_session* session);
