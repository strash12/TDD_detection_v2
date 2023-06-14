#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>


bool dmac_read(size_t sample_count, size_t duration_ms, void *user_data, void (*data_reader)(const void *data, size_t data_size, const char *error, void *user_data));
