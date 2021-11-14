#define BUFFER_SIZE 1024
char error_buf[BUFFER_SIZE];
char error_str[BUFFER_SIZE];

#define log_error(MSG, ...)                                                    \
  {                                                                            \
    snprintf(error_str, (BUFFER_SIZE - 1), "[ERROR] (%s:%s:%i) ", __FILE__,    \
             __func__, __LINE__);                                              \
    snprintf(error_buf, (BUFFER_SIZE - 1), MSG, ##__VA_ARGS__);                \
    strcat(error_str, error_buf);                                              \
    puts(error_str);                                                           \
  }

