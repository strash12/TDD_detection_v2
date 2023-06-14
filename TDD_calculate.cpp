#include "header.h"


namespace mark
{
    TDD_create::TDD_create(param_store::_parametrs param)
    {
        _param = param;
    }

    TDD_create::~TDD_create(){};

    int TDD_create::Calculate_TDD()
    {

            xvec sig_part_2 = set_mark(_param.path_second_part,"fr2.txt");
            xvec sig_part_1 = set_mark(_param.path_first_part,"fr1.txt");
            int first_1 = second_corr(sig_part_1);
            int first_2 = second_corr(sig_part_2);
		std::cout<<"first 1 = "<<first_1<<std::endl;
std::cout<<"first 2 = "<<first_2<<std::endl;
            double angle_1 = calculate_angle(sig_part_1,first_1);
            double angle_2 = calculate_angle(sig_part_2,first_2);
std::cout<<"angle_1 = "<<angle_1<<std::endl;
std::cout<<"angle_2 = "<<angle_2<<std::endl;
            xvec sig_shift_1 = freq_shift_remove(sig_part_1,angle_1);
            xvec sig_shift_2 = freq_shift_remove(sig_part_2,angle_2);
            RS::RS_cut RS_proc(_param);
            return RS_proc.set_RS(sig_shift_1,sig_shift_2,first_1,first_2);
        

    }

    xvec TDD_create::set_mark(std::string path[3],std::string name)
    {
        int ERROR_count;
        bool validate_bit = 0;
        int i = 0;

        fpga_configure::SSS_upload SSS(_param);
        fpga_configure::config_fpga check(_param);

        while (!validate_bit)
        {
            SSS.load_SSS(_param.CellId,path[i++]);
            SSS.shift_mark(_param.shift_mark);
	        usleep(10000000);
            validate_bit = check.read_validate_bit(ERROR_count);
            if(_param.band > 4 && i == 3 ){i = 0;}
	        else if(_param.band < 5 && i == 1) {i = 0;}

            if(ERROR_count == 10 )
            {
                tch::write_common("SSS not found");
                throw mark::TDDException("SSS signal not found, try again");
                break;
            }
        }

        signal::signal_create download_sig(_param.PSS_size);
        xvec signal = download_sig.download_convert_signal(name,_param.interp_coef);
        return signal;
    }

    int TDD_create::second_corr(xvec signal)
    {
        LTE::xcorr corr;
        xvec PSS = PSS_generate::pss::get_pss(_param.PSS_number,_param.fftsize,_param.cp);
        return corr.corr_coef(signal,PSS,_param.refsamples,_param.fftsize,_param.cp,_param.shift_second_cor_start,_param.shift_second_cor_stop);
    }

    
    double TDD_create::calculate_angle(xvec signal, int first)
    {
        LTE::freq_shift_estimate freq;
        return freq.get_angle(signal,_param.cp,_param.fftsize,first);
    }

    
    xvec TDD_create::freq_shift_remove(xvec signal,double angle)
    {
        xvec sig(signal.size());
        for(double i = 0; i <signal.size();i++)
        {
            sig[i] = signal[i] * exp(J*(i)*(angle/_param.fftsize));
        }
        return sig;
    }

}