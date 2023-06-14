#ifndef DEVMEM_HPP
#define DEVMEM_HPP

#include <exception>
#include <string>
extern "C" {
	#include "devmem.h"
}


namespace devmem {
	class DevmemBaseException : public std::exception { };

	class DevmemCreateException : public DevmemBaseException { };
	class DevmemCreateNotFoundException : public DevmemCreateException { };
	class DevmemCreatePermissionDeniedException : public DevmemCreateException { };
	class DevmemCreateAlreadyExistsException : public DevmemCreateException { };
	class DevmemInvalidInputException : public DevmemCreateException { };
	class DevmemCreateUnknownException : public DevmemCreateException { };

	class DevmemWriteException : public DevmemBaseException { };
	class DevmemWriteOutOfRangeException : public DevmemWriteException { };

	class DevmemReadException : public DevmemBaseException { };
	class DevmemReadOutOfRangeException : public DevmemReadException { };


	class Devmem {
	private:
		::Devmem* _devmem;

	public:
		inline Devmem(const std::string& filename, size_t physical_addr, size_t len)
			: _devmem(devmem_create(filename.c_str(), physical_addr, len)) {
				if (reinterpret_cast<size_t>(_devmem) <= DevmemCreateError_Unknown) {
					DevmemCreateError err = devmem_create_get_err(_devmem);
					switch (err)
					{
					case DevmemCreateError_NotFound:
						throw DevmemCreateNotFoundException();
						break;
					case DevmemCreateError_PermissionDenied:
						throw DevmemCreatePermissionDeniedException();
						break;
					case DevmemCreateError_AlreadyExists:
						throw DevmemCreateAlreadyExistsException();
						break;
					case DevmemCreateError_InvalidInput:
						throw DevmemInvalidInputException();
						break;
					case DevmemCreateError_Unknown:
						throw DevmemCreateUnknownException();
						break;
					}
				}
			}

		~Devmem() { devmem_destroy(_devmem); }

		template<typename T>
		inline void write(size_t offset, T value) {
			const void* res = devmem_write(_devmem, offset, reinterpret_cast<const void*>(&value), sizeof(value));
			if (reinterpret_cast<size_t>(res) <= DevmemWriteError_OutOfRange) {
				throw DevmemWriteOutOfRangeException();
			}
		}

		template<typename T>
		inline T read(size_t offset = 0) {
			const void* res = devmem_read(_devmem, offset, sizeof(T));
			if (reinterpret_cast<size_t>(res) <= DevmemReadError_OutOfRange) {
				throw DevmemReadOutOfRangeException();
			}
			return *reinterpret_cast<const T*>(res);
		}

		
	};
}
#endif // DEVMEM_HPP