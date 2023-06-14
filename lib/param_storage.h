#ifndef defined PARAM_STORAGE
#define PARAM_STORAGE

#include "../../TDD_detection/header.h"
#include <string>


namespace param_store
{
    struct _parametrs
        {
            int fftsize,cp,cp1,Nrb,Lb,i1start,i1fin,i2start,i2fin,Nc,Nsc,Nrb_Max,Nrb_Min,Ns,OSC,fs,Cor_lvl,cyclic,windowing,band,shift,first,PSS_number,CellId,shift_mark,delay,refsamples,interp_coef,PSS_size;
	    int shift_second_cor_start, shift_second_cor_stop;
	    std::string path_first_part[3],path_second_part[3];
            uint32_t MMCM, MMCM_2,adreses;
	
        } ;
        
    class parametrs
    {
        private:
             _parametrs Norm_cp(int);
             _parametrs Ext_cp(int);
        public:
            parametrs();
            ~parametrs();
            _parametrs set_param(tch::FullTusurCommonConfig _config);
    };
} //namespace param


#endif //PARAM_STORAGE