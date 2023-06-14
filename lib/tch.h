#ifndef TCH_H
#define TCH_H

#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>


typedef enum CpMode {
  Normal,
  Extended,
} CpMode;

typedef enum CorrelationKind {
  PSS,
  SSS,
} CorrelationType;

typedef struct TCH_FullTusurCommonConfig {
  float bandwidth;
  enum CpMode cp_mode;
  enum CorrelationKind correlation_kind;
  size_t correlation_level;
  bool automatic_threshold;
} FullTusurCommonConfig;


void init(const char* python_path, const char* log_level);
void write_common(const char* msg);
void run(const char* radius_control_software_path, void(*on_handle)(const char* name, const void* input, void* output, void* user_data), void* user_data);


#endif //TCH_H