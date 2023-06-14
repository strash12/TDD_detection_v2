#ifndef TCH_HPP
#define TCH_HPP

#include <string>
#include <thread>
#include <functional>
#include <unordered_map>


namespace tch {
	namespace internal {
		extern "C" {
			#include "tch.h"
		}
	}

	class BaseHandle;
	void run(const std::string& radius_control_software_path, BaseHandle* handle);

	enum class CpMode {
		Normal,
		Extended
	};
	
	enum class CorrelationKind {
		PSS,
		SSS
	};

	struct FullTusurCommonConfig {
		float bandwidth = 1.92;
		CpMode cp_mode = CpMode::Normal;
		CorrelationKind correlation_kind = CorrelationKind::PSS;
		size_t correlation_level = 0x020100DC;
		bool automatic_threshold = false;
	};
	
	class BaseHandle {
		friend void run(const std::string& radius_control_software_path, BaseHandle* handle);

	public:
		virtual ~BaseHandle() {}

	protected:
		virtual void handle_set_bandwidth(float bandwidth) { printf("ffff\n"); }
		virtual float handle_get_bandwidth() {}

		virtual void handle_set_cp_mode(CpMode mode) {}
		virtual CpMode handle_get_cp_mode() {}

		virtual void handle_set_correlation_kind(CorrelationKind kind) {}
		virtual CorrelationKind handle_get_correlation_kind() {}
		
		virtual void handle_set_correlation_level(size_t level) {}
		virtual size_t handle_get_correlation_level() {}

		virtual void handle_set_automatic_threshold(bool enable) {}
		virtual bool handle_get_automatic_threshold() {}

		virtual FullTusurCommonConfig handle_configure() {}

	private:
		static void on_handle(const char* name, const void* input, void* output, BaseHandle* self) {
			static const std::unordered_map<std::string, std::function<void(BaseHandle* self, const void* input, void* output)>> string_switch = {
				{
					"set_bandwidth",
					[](BaseHandle* self, const void* input, void* output) {
						self->handle_set_bandwidth(*reinterpret_cast<const float*>(input));
					}
				},
				{
					"get_bandwidth",
					[](BaseHandle* self, const void* input, void* output) {
						*(reinterpret_cast<float*>(output)) = self->handle_get_bandwidth();
					}
				},
				{
					"set_cp_mode",
					[](BaseHandle* self, const void* input, void* output) {
						self->handle_set_cp_mode(*reinterpret_cast<const CpMode*>(input));
					}
				},
				{
					"get_cp_mode",
					[](BaseHandle* self, const void* input, void* output) {
						*(reinterpret_cast<CpMode*>(output)) = self->handle_get_cp_mode();
					}
				},
				{
					"set_correlation_kind",
					[](BaseHandle* self, const void* input, void* output) {
						self->handle_set_correlation_kind(*reinterpret_cast<const CorrelationKind*>(input));
					}
				},
				{
					"get_correlation_kind",
					[](BaseHandle* self, const void* input, void* output) {
						*(reinterpret_cast<CorrelationKind*>(output)) = self->handle_get_correlation_kind();
					}
				},
				{
					"set_correlation_level",
					[](BaseHandle* self, const void* input, void* output) {
						self->handle_set_correlation_level(*reinterpret_cast<const size_t*>(input));
					}
				},
				{
					"get_correlation_level",
					[](BaseHandle* self, const void* input, void* output) {
						*(reinterpret_cast<size_t*>(output)) = self->handle_get_correlation_level();
					}
				},
				{
					"set_automatic_threshold",
					[](BaseHandle* self, const void* input, void* output) {
						self->handle_set_automatic_threshold(*reinterpret_cast<const bool*>(input));
					}
				},
				{
					"get_automatic_threshold",
					[](BaseHandle* self, const void* input, void* output) {
						*(reinterpret_cast<bool*>(output)) = self->handle_get_automatic_threshold();
					}
				},
				{
					"configure",
					[](BaseHandle* self, const void* input, void* output) {
						*(reinterpret_cast<FullTusurCommonConfig*>(output)) = self->handle_configure();
					}
				}
			};

			std::string handle_name = name;
			auto search = string_switch.find(name);
			if (search != string_switch.end()) { search->second(self, input, output); }
		}
	};

	inline void init(const std::string& python_path, const std::string& log_level) {
		internal::init(python_path.c_str(), log_level.c_str());
	}

	inline void write_common(const std::string& msg) {
		internal::write_common(msg.c_str());
	}

	inline void run(const std::string& radius_control_software_path, BaseHandle* handle) {
		internal::run(radius_control_software_path.c_str(), reinterpret_cast<void(*)(const char*, const void*, void*, void*)>(BaseHandle::on_handle), handle);
	}
}


#endif //TCH_HPP