#pragma once

#include <stddef.h>
#include <stdlib.h>

#define generic_alloc(type, size) (type*)malloc(sizeof(type) * size)
#define generic_free(x)           ({ free(x); })

#define arena_alloc(allocator, type, size) \
  (type*)alloc_arena(allocator, sizeof(type), _Alignof(type), size)

typedef struct arena_allocator arena_allocator;
struct arena_allocator {
  char*     start;
  char*     end;
  ptrdiff_t capacity;
};

arena_allocator create_arena(ptrdiff_t size);
void*           alloc_arena(arena_allocator* a, ptrdiff_t size, ptrdiff_t align,
                            ptrdiff_t count);
ptrdiff_t       get_arena_free_mem(arena_allocator* a);
void            free_arena(arena_allocator* a);
