#ifndef defined FPGA_CONFIGURE
#define FPGA_CONFIGURE

#include <unistd.h>
#include "../../TDD_detection/header.h"

namespace fpga_configure
{
    class config_fpga
    {
        public:
            param_store::_parametrs _param;
            config_fpga(param_store::_parametrs param);
            ~config_fpga();

            void configure(uint32_t PSS);
            void set_PSS_fpga();
            bool read_validate_bit(int& Error_count);
        private:
            void MMCM();
            void config_corr();
            
    
    };//class config_fpga
    

    class SSS_upload
    {
        public:
            param_store::_parametrs _param;
            SSS_upload(param_store::_parametrs param);
            ~SSS_upload();
            void load_SSS(int CellId,std::string path);
            void shift_mark(int shiftt);
        private:


    };//class SSS_upload





}// namespace fpga_configure






#endif //FPGA_CONFIGURE