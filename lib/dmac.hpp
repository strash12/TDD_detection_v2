#include <vector>
#include <string>
#include <exception>
#include <iostream>
extern "C" {
	#include "dmac.h"
}

namespace dmac {
	class Exception: public std::exception {
	private:
		std::string _msg;

	public:
		Exception(): _msg() {}
		Exception(std::string msg): _msg(msg) {}

		const char* what() const throw() { return _msg.c_str(); }
	};

	namespace internal {
		template<typename T>
		struct UserData {
			std::vector<T> data;
			Exception exception;
		};

		template<typename T>
		void reader(const void *data, size_t data_size, const char *error, void *user_data) {
			UserData<T>* _user_data = reinterpret_cast<UserData<T>*>(user_data);

			if (error == nullptr) {
				const T* _data = reinterpret_cast<const T*>(data);
				size_t _data_size = data_size / sizeof(T);
				_user_data->data.reserve(_data_size);
				_user_data->data.insert(_user_data->data.begin(), &_data[0], &_data[_data_size]);
			} else {
				_user_data->exception = Exception(std::string(error));
			}
		}
	}
	
	template<typename T>
	std::vector<T> read(size_t sample_count, size_t duration_ms) {
		internal::UserData<T> user_data{std::vector<T>(), Exception()};
		
		bool res = dmac_read(sample_count, duration_ms, reinterpret_cast<void*>(&user_data), &internal::reader<T>);
		if (res) {
			return user_data.data;
		} else {
			throw user_data.exception;
		}
	}
}
