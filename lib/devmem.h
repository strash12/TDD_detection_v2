#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>


typedef struct Devmem Devmem;

typedef enum DevmemCreateError {
  DevmemCreateError_NotFound = 0,
  DevmemCreateError_PermissionDenied = 1,
  DevmemCreateError_AlreadyExists = 2,
  DevmemCreateError_InvalidInput = 3,
  DevmemCreateError_Unknown = 4,
} DevmemCreateError;

typedef enum DevmemReadError {
  DevmemReadError_OutOfRange = 0,
} DevmemReadError;

typedef enum DevmemWriteError {
  DevmemWriteError_OutOfRange = 0,
} DevmemWriteError;



Devmem *devmem_create(const char *filename, size_t physical_addr, size_t len);
enum DevmemCreateError devmem_create_get_err(const Devmem *result);

void devmem_destroy(Devmem *devmem);

const void *devmem_write(Devmem *devmem, size_t offset, const void *data, size_t len);
enum DevmemWriteError devmem_write_get_err(const void *result);

const void *devmem_read(const Devmem *devmem, size_t offset, size_t len);
enum DevmemReadError devmem_read_get_err(const void *result);


