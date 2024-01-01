#pragma once

#include <json-c/json.h>
#include <stdio.h>

static json_object* get_key_js(json_object* js, const char* key) {
  return json_object_object_get(js, key);
}

#define same_as(lhs_t, rhs_t)                       \
  ({                                                \
    __builtin_types_compatible_p(lhs_t, rhs_t) &&   \
        __builtin_types_compatible_p(rhs_t, lhs_t); \
  })

#define get_key(type, js, key)                             \
  ({                                                       \
    json_object* obj_js = json_object_object_get(js, key); \
    type         value;                                    \
    if (same_as(type, bool)) {                             \
      value = json_object_get_boolean(obj_js);             \
    } else if (same_as(type, f64)) {                       \
      value = json_object_get_double(obj_js);              \
    } else {                                               \
      abort();                                             \
    }                                                      \
    value;                                                 \
  })

static json_object* read_json(const char* filepath) {
  json_object* obj = json_object_from_file(filepath);
  if (!obj) {
    fprintf(stderr, "Failed to read and JSON file [%s]. Exiting...\n",
            filepath);
    exit(EXIT_FAILURE);
  }
  return obj;
}
