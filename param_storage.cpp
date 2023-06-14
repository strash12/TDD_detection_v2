#include "header.h"

namespace param_store
{
    param_store::parametrs::parametrs(){}
    param_store::parametrs::~parametrs(){}

    
    _parametrs parametrs::Ext_cp(int band)
    {
    _parametrs param;    
    switch (band)
    {
    case 1:
        param.fftsize       = 1024;
        param.band          = 4;
        param.cp1           = 256;
        param.cp            = 256;
        param.fs            = 15.36e6;
        param.Nrb           = 50;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param.fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;               
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 600;
        param.PSS_size      = 1280;
        param.adreses      = 0x000401D2;
        param .MMCM         = 0x00000010;
        param .MMCM_2       = 0x00002004;
        param.shift_mark         = param.fs/100-param.fs/1000-6*param.fftsize;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B4.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B4_n2.txt";
	    param.shift_second_cor_start = 4700;
	    param.shift_second_cor_stop = 5000;
        param.windowing     = 5092;
        param.cyclic        = 2;
        param.shift         = 50;
        param.delay         = -14;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;
    break;

    
    case 3:
        param.fftsize        = 256;
        param.band           = 2;
        param.cp1            = 64;
        param.cp             = 64;
        param.fs             = 3.84e6;
        param.Nrb            = 15;
        param.i1start        = param.fftsize/2-31;
        param.i1fin          = param.fftsize/2-1;
        param.i2start        = param.fftsize/2+1;
        param.i2fin          = param.fftsize/2+31;
        param.Nc             = 1600;             
        param.Nsc            = 12;               
        param.Nrb_Max        = 110;            
        param.Nrb_Min        = 6;             
        param.Ns             = 20;
        param.OSC            = 150;
        param.PSS_size      = 320;
        param.adreses       = 0x000200D2;
        param .MMCM          = 0x00000040;
        param .MMCM_2        = 0x00002004;
        param.shift_mark         = param.fs/100-param.fs/1000-4 * param.fftsize;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B2.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B2_n2.txt";
	    param.shift_second_cor_start = -200;
	    param.shift_second_cor_stop = 20;
        param.windowing     = 1242;
        param.cyclic        = 2;
        param.shift         = 50;
        param.delay         = -14;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;
    break;

    case 7:
        param.fftsize       = 512;
        param.band          = 3;
        param.cp1           = 128;
        param.cp            = 128;
        param.fs            = 7.68e6;
        param.Nrb           =  25;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param.fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;               
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 300;
        param.PSS_size      = 640;
        param.adreses     = 0x00030136;
        param .MMCM         = 0x00000020;
        param .MMCM_2       = 0x00002004;
        param.shift_mark         = param.fs/100-param.fs/1000- 4 * param.fftsize;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B3.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B3_n2.txt";
	    param.shift_second_cor_start = 800;
	    param.shift_second_cor_stop = 1000;
        param.windowing     = 2528;
        param.cyclic        = 2;
        param.shift         = 50;
        param.delay         = -14;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;
    break;

    case 15:
        param.fftsize       = 1024;
        param.band          = 4;
        param.cp1           = 256;
        param.cp            = 256;
        param.fs            = 15.36e6;
        param.Nrb           = 50;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param.fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;               
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 600;
        param.PSS_size      = 1280;
        param.adreses      = 0x000401D2;
        param .MMCM         = 0x00000010;
        param .MMCM_2       = 0x00002004;
        param.shift_mark         = param.fs/100-param.fs/1000-param.fftsize;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B4.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B4_n2.txt";
	    param.shift_second_cor_start = -500;
	    param.shift_second_cor_stop = -200;
        param.windowing     = 5099;
        param.cyclic        = 2;
        param.shift         = 50;
        param.delay         = -14;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;
    break;

    case 23:
        param.fftsize       = 1536;
        param.band          = 5;
        param.cp1           = 384;
        param.cp            = 384;
        param.fs            = 23.04e6;
        param.Nrb           = 75;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param. fftsize/2+1;
        param.i2fin         = param. fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;               
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 900;
        param.PSS_size      = 2560;
        param.adreses       = 0x000501A4;
        param .MMCM         = 0x0002FE0A;
        param .MMCM_2       = 0x01702004;
        param.shift_mark         = param.fs/100-param.fs/1000-param.fftsize;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B5_n1_1.txt";
        param.path_first_part[1]   = "./signal_storage/SSS_array_B5_n1_2.txt";
	    param.path_first_part[2]   = "./signal_storage/SSS_array_B5_n1_3.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B5_n2_1.txt";
        param.path_second_part[1]   = "./signal_storage/SSS_array_B5_n2_2.txt";
	    param.path_second_part[2]   = "./signal_storage/SSS_array_B5_n2_3.txt";
	    param.shift_second_cor_start = -200;
	    param.shift_second_cor_stop = 200;
        param.windowing     = 7645;
        param.cyclic        = 2;
        param.shift         = 50;
        param.delay         = -14;
        param.refsamples    = param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;
    break;

    case 30:
        param.fftsize       = 2048;
        param.band          = 6;
        param.cp1           = 512;
        param.cp            = 512;
        param.fs            = 30.72e6;
        param.Nrb           = 100;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param. fftsize/2+1;
        param.i2fin         = param. fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;               
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 1200;
        param.PSS_size      = 2560;
        param.adreses       = 0x00060800;
        param .MMCM         = 0x00000008;
        param .MMCM_2       = 0x00002004;
        param.shift_mark         = param.fs/100-param.fs/1000-param.fftsize;
        param.path_first_part [0]   = "./signal_storage/SSS_array_B6_n1_1.txt";
        param.path_first_part [1]   = "./signal_storage/SSS_array_B6_n1_2.txt";
        param.path_first_part [2]   = "./signal_storage/SSS_array_B6_n1_3.txt";
        param.path_second_part [0]   = "./signal_storage/SSS_array_B6_n2_1.txt";
        param.path_second_part [1]   = "./signal_storage/SSS_array_B6_n2_2.txt";
        param.path_second_part [2]   = "./signal_storage/SSS_array_B6_n2_3.txt";
	    param.shift_second_cor_start = -200;
	    param.shift_second_cor_stop = 200;
        param.windowing     = 10000;
        param.cyclic        = 2;
        param.shift         = 1000;
        param.delay         = -14;
        param.refsamples    = param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;
        break;
    
    default:
        mark::TDDException("Error set parametrs to param_storage");
    break;
    }
    return param;
    }

    _parametrs parametrs::Norm_cp(int band)
    {
        _parametrs param;
        switch (band)
        {
    case 1:
        param.fftsize       = 1024;
        param.band          = 4;
        param.cp1           = 80;
        param.cp            = 72;
        param.fs            = 15.36e6;
        param.Nrb           = 50;
        param.Lb            = 7;
        param.i1start       = param. fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param. fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;             
        param.Nrb_Max       = 110;           
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 600;
        param.PSS_size      = 1096;
        param.adreses     = 0x000401D2;
        param .MMCM         = 0x00000010;
        param .MMCM_2       = 0x00002004;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B4.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B4_n2.txt";
        param.shift_mark    = param.fs/100-param.fs/1000-6*param.fftsize;
	    param.shift_second_cor_start = 2500;
	    param.shift_second_cor_stop = 2800;
        param.windowing     = 5265;
        param.cyclic        = 1;
        param.shift         = 50;
        param.delay         = 2402;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;

    break;
    
    case 3:
        param.fftsize       = 256;
        param.band          = 2;
        param.cp1           = 20;
        param.cp            = 18;
        param.fs            = 3.84e6;
        param.Nrb           = 15;
        param.Lb            = 7;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param.fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;              
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 150;
        param.PSS_size      = 274;
        param.adreses      = 0x000200D2;
        param .MMCM         = 0x00000040;
        param .MMCM_2       = 0x00002004;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B2.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B2_n2.txt";
        param.shift_mark    = param.fs/100-param.fs/1000- 4 *param.fftsize;
	    param.shift_second_cor_start = -3000;
	    param.shift_second_cor_stop = 3000;
        param.windowing     = 1288;
        param.cyclic        = 1;
        param.shift         = 50;
        param.delay         = 2402;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;


    break;
    
    case 7:
        param.fftsize       = 512;
        param.band          = 3;
        param.cp1           = 40;
        param.cp            = 36;
        param.fs            = 7.68e6;
        param.Nrb           = 25;
        param.Lb            = 7;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param.fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;             
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;            
        param.Ns            = 20;
        param.OSC           = 300;
        param.PSS_size      = 548;
        param.adreses     = 0x00030136;
        param.MMCM         = 0x00000020;
	    param.MMCM_2       = 0x00002004;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B3.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B3_n2.txt";
        param.shift_mark    = param.fs/100-param.fs/1000- 3 * param.fftsize;
	    param.shift_second_cor_start = -2000;
	    param.shift_second_cor_stop = -1800;
        param.windowing     = 2619;
        param.cyclic        = 1;
        param.shift         = 50;
        param.delay         = 2402;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;

    break;
    
    case 15:
        param.fftsize       = 1024;
        param.band          = 4;
        param.cp1           = 80;
        param.cp            = 72;
        param.fs            = 15.36e6;
        param.Nrb           = 50;
        param.Lb            = 7;
        param.i1start       = param. fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param. fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;             
        param.Nrb_Max       = 110;           
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 600;
        param.PSS_size      = 1096;
        param.adreses     = 0x000401D2;
        param .MMCM         = 0x00000010;
        param .MMCM_2       = 0x00002004;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B4.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B4_n2.txt";
	    param.shift_second_cor_start = -2700;
	    param.shift_second_cor_stop = -2400;
        param.shift_mark    = param.fs/100-param.fs/1000-param.fftsize;
        param.windowing     = 5265;
        param.cyclic        = 1;
        param.shift         = 50;
        param.delay         = 2402;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;

    break;

    case 23:
        param.fftsize       = 1536;
        param.band          = 5;
        param.cp1           = 120;
        param.cp            = 108;
        param.fs            = 23.04e6;
        param.Nrb           = 100;
        param.Lb            = 7;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param.fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;               
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 1200;
        param.PSS_size      = 2192;
        param.adreses       = 0x000501A4;
        param .MMCM         = 0x0002FE0A;
        param .MMCM_2       = 0x01702004;
        param.shift_mark         = param.fs/100-param.fs/1000-param.fftsize;
        param.path_first_part[0]   = "./signal_storage/SSS_array_B5.txt";
        param.path_first_part[0]   = "./signal_storage/SSS_array_B5_n1_1.txt";
        param.path_first_part[1]   = "./signal_storage/SSS_array_B5_n1_2.txt";
	    param.path_first_part[2]   = "./signal_storage/SSS_array_B5_n1_3.txt";
        param.path_second_part[0]   = "./signal_storage/SSS_array_B5_n2_1.txt";
        param.path_second_part[1]   = "./signal_storage/SSS_array_B5_n2_2.txt";
	    param.path_second_part[2]   = "./signal_storage/SSS_array_B5_n2_3.txt";
	    param.shift_second_cor_start = -3000;
	    param.shift_second_cor_stop = -1500;
        param.windowing     = 7938;
        param.cyclic        = 1;
        param.shift         = 50;
        param.delay         = 2042;
	    param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;
    break;

    case 30:
        param.fftsize       = 2048;
        param.band          = 6;
        param.cp1           = 160;
        param.cp            = 144;
        param.fs            = 30.72e6;
        param.Nrb           = 100;
        param.Lb            = 7;
        param.i1start       = param.fftsize/2-31;
        param.i1fin         = param.fftsize/2-1;
        param.i2start       = param.fftsize/2+1;
        param.i2fin         = param.fftsize/2+31;
        param.Nc            = 1600;             
        param.Nsc           = 12;               
        param.Nrb_Max       = 110;            
        param.Nrb_Min       = 6;             
        param.Ns            = 20;
        param.OSC           = 1200;
        param.PSS_size      = 2192;
        param.adreses     = 0x000601AE;
        param .MMCM         = 0x00000008;
        param .MMCM_2       = 0x00002004;
        param.path_first_part [0]   = "./signal_storage/SSS_array_B6_n1_1.txt";
        param.path_first_part [1]   = "./signal_storage/SSS_array_B6_n1_2.txt";
        param.path_first_part [2]   = "./signal_storage/SSS_array_B6_n1_3.txt";
        param.path_second_part [0]   = "./signal_storage/SSS_array_B6_n2_1.txt";
        param.path_second_part [1]   = "./signal_storage/SSS_array_B6_n2_2.txt";
        param.path_second_part [2]   = "./signal_storage/SSS_array_B6_n2_3.txt";
        param.shift_mark    = param.fs/100-param.fs/1000-param.fftsize;
	    param.shift_second_cor_start = -2000;
	    param.shift_second_cor_stop = -1000;
        param.windowing     = 10000;
        param.cyclic        = 1;
        param.shift         = 1000;
        param.delay         = 2402;
        param.refsamples    =   param.fs/1000+3*param.fftsize+param.cp1+2*param.cp+param.delay;

    break;
    
    default:
        mark::TDDException("Error set parametrs to param_storage");
    break;
    }
    return param;
    }

         
    _parametrs  parametrs::set_param(tch::FullTusurCommonConfig _config)
    {
        int band = _config.bandwidth;
        if(_config.cp_mode == tch::CpMode::Normal)
            {   
                return Norm_cp(band);
            }
        else
            {
                return Ext_cp(band);
            }
            
    }

    

}//namespace param




