#include <fstream>
#include "header.h"
#include "bitset"

namespace fpga_configure
{
    devmem::Devmem MMCM_conf_1 = devmem::Devmem("/dev/mem", 0x83C30208, 4);
    devmem::Devmem MMCM_conf_2 = devmem::Devmem("/dev/mem", 0x83C30200, 4);
    devmem::Devmem MMCM_upd = devmem::Devmem("/dev/mem", 0x83C3025C, 2);
    devmem::Devmem COR_init = devmem::Devmem("/dev/mem", 0x83C20000, 2);
    devmem::Devmem COR_ctrl = devmem::Devmem("/dev/mem", 0x83C20004, 4);
    devmem::Devmem COR_tele = devmem::Devmem("/dev/mem", 0x83C20008, 4);

    int time = 10000000;
    
   config_fpga::config_fpga(param_store::_parametrs param)
    {
        _param = param;
    }
    config_fpga::~config_fpga() {}


    void config_fpga::MMCM()
    {
        fpga_configure::MMCM_conf_2.write<uint32_t>(0, _param.MMCM_2);
        fpga_configure::MMCM_conf_1.write<uint32_t>(0, _param.MMCM);
        fpga_configure::MMCM_upd.write<uint16_t>(0, 0x3u);
    }

    void config_fpga::config_corr()
    {
        uint32_t config = _param.adreses;
        fpga_configure::COR_ctrl.write(0, config);
        uint32_t automatic_on = (config & 0xffff0fff) | (1 << 12);
        fpga_configure::COR_ctrl.write<uint32_t>(0, automatic_on);
    }


    void config_fpga::configure(uint32_t PSS)
    {
        fpga_configure::COR_init.write<uint8_t>(0, 0x1u);
        fpga_configure::COR_ctrl.write<uint32_t>(0, PSS);
        fpga_configure::COR_init.write<uint8_t>(0, 0x3u);
        fpga_configure::COR_init.write<uint8_t>(0, 0x1u);
        usleep(time);

    }

    void config_fpga::set_PSS_fpga()
    {
        MMCM();
        config_corr();
        uint32_t PSS = _param.adreses;
        PSS = (PSS & 0xffff0fff) | (1 << 12);
        bool validate_bit;
        int ERROR_count = 0;
        
        
        for (int i = 1; i < 4; i++)
        {
            
            PSS = (PSS & 0xf0ffffff) | (i << 24);
            configure(PSS);
            validate_bit = read_validate_bit(ERROR_count);
            if (validate_bit == 1) 
            { 
                tch::write_common("PSS found");
                break; 
            }
            
            if (ERROR_count == 9) 
            { 
                tch::write_common("bad signal try to correlate");
                PSS = (PSS & 0x0fffffff) |(1<<29);
                fpga_configure::COR_ctrl.write(0,PSS);
                time = 20000000;
            }

            if (ERROR_count == 20 && validate_bit == 0 ) 
            { 
                tch::write_common("PSS not found, try again");
                throw mark::TDDException("PSS signal not found");
                break; 
            }

            if (i == 3) { i = 0; }
        }

    }

    bool config_fpga::read_validate_bit(int& Error_count)
    {
	
        float div_cor = 0;
        bool verefi_bit = 0, thr_hold = 0, asa = 0;
        uint32_t bits;
        for (int i = 0; i < 1000000; i++)
        {
            bits = 0x00000000;
            bits = fpga_configure::COR_tele.read<uint32_t>();
            verefi_bit = (bits >> 28) & 1;
            if (verefi_bit == 1)
            {
                return 1;
                break;
            }	
        }
	Error_count ++ ;    
        return 0;
    }

    
    SSS_upload::SSS_upload(param_store::_parametrs param)
    {
        _param = param;
    }

    SSS_upload::~SSS_upload() {};


    void SSS_upload::load_SSS(int CellId, std::string path)
    {
        uint32_t config = _param.adreses;
        int numOfLines = 0;
        std::string digits;
        COR_ctrl.write<uint32_t>(0, config);
        config = (config & 0xffff0fff) | (1 << 12);
        COR_ctrl.write<uint32_t>(0, config);
        config = (config & 0x0fffffff) | 1 << 28;
        COR_ctrl.write<uint32_t>(0, config);

        std::ifstream file;
        file.open(path);
        if (!file.is_open()) {
            std::cout << "Error on open SSS file" << std::endl;
        }

        while (numOfLines != CellId+1 && getline(file, digits)) {
            numOfLines++;
        }

        for (int i = 0; i < 1024; i++) {
            uint32_t SSS = 0x00000000;
            if (digits[i] > 0x39) {
                SSS = (SSS & 0xfffffff0) | (digits[i] - 0x37);
            }
            else {
                SSS = (SSS & 0xfffffff0) | (digits[i] - 0x30);
            }

            COR_tele.write<uint32_t>(0, SSS);
            COR_init.write<uint8_t>(0, 0x1u);
            COR_init.write<uint16_t>(0, 0x3u);
            COR_init.write<uint16_t>(0, 0x1u);

        }


    }


    void SSS_upload::shift_mark(int shiftt)
    {
        auto shift = devmem::Devmem("/dev/mem", 0x83C2000C, 4);
        uint32_t _shift = 0x00000000;
        shift.write<uint32_t>(0, _shift);
        _shift = 0x00000000 | (shiftt);
        shift.write<uint32_t>(0, _shift);
    }


}

