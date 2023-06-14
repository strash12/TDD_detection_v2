#ifndef header
#define header


#include "../TDD_detection/lib/tch.hpp"
#include "../TDD_detection/lib/devmem.hpp"
#include "../TDD_detection/lib/dmac.hpp"
//#include "../TDD_detection/lib/logger.hpp"
#include "../TDD_detection/lib/param_storage.h"
#include "../TDD_detection/lib/signal_create.h"
#include "../TDD_detection/lib/fpga_configure.h"
#include "../TDD_detection/lib/LTE_PROC.h"
#include "../TDD_detection/lib/SSS_PROC.h"
#include "../TDD_detection/lib/RS_PROC.h"
#include "../TDD_detection/lib/NRS_gen.h"
#include "../TDD_detection/lib/pss.h"

#include <vector>
#include <complex>
#include <string_view>
#include <string>
#include <exception>

typedef std::complex<double> xd;
typedef std::vector<xd> xvec;

namespace mark
{

    class TDDException:public std::exception
    {
        private:
            std::string m_error{};
        public:
            TDDException(std::string error)
                : m_error{error}
        {
        }
        const char* what() const noexcept override { return m_error.c_str(); }
    };

    class create_mark
    {

        public:
            create_mark(tch::FullTusurCommonConfig config);
            ~create_mark();
            param_store::_parametrs _param;
            tch::FullTusurCommonConfig _config;
            LTE::cor_data PSS_data;
            void calculate_mark(); 
        private:
            void set_params();
            void download_convert_signal(xvec& signal);
            void save_signal_to_file(std::vector<int32_t> signal);
            void corelate_signal(xvec signal);
            void search_maximum();
            void cell_id_calculate(xvec siganal);
    }; // class create_mark

    class TDD_create
    {
        public:
            TDD_create(param_store::_parametrs param);
            param_store::_parametrs _param;
            ~TDD_create();
            int Calculate_TDD();
        private:
            xvec set_mark(std::string path[3],std::string name);
            int second_corr(xvec signal);
            double calculate_angle(xvec signal,int first);
            xvec freq_shift_remove(xvec signal,double angle);
    };// class TDD_create

} // namespace mark


#endif //header