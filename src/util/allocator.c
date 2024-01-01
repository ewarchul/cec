#include "cec/allocator.h"

#include <stdint.h>
#include <string.h>

arena_allocator create_arena(ptrdiff_t capacity) {
  arena_allocator a = {0};
  a.start = generic_alloc(char, capacity);
  a.end   = a.start ? a.start + capacity : 0;
  a.capacity = capacity;
  return a;
}

void* alloc_arena(arena_allocator* a, ptrdiff_t size, ptrdiff_t align,
                      ptrdiff_t count) {
  ptrdiff_t avail   = a->end - a->start;
  ptrdiff_t padding = -(uintptr_t)a->start & (align - 1);
  if (count > (avail - padding) / size) {
    abort();
  }
  ptrdiff_t total      = size * count;
  char*     memory_ptr = a->start + padding;
  a->start += padding + total;
  return memset(memory_ptr, 0, total);
}

ptrdiff_t get_arena_free_mem(arena_allocator* a) {
  return a->end - a->start;
}

void free_arena(arena_allocator* a) {
  char* p = a->end - a->capacity;
  free(p);
}
