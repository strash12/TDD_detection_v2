#include "header.h"

namespace mark
{
        

        create_mark::create_mark(tch::FullTusurCommonConfig config)
        {
            _config = config;
        };
        create_mark::~create_mark(){};

        void create_mark::calculate_mark()
        {

                set_params();
                signal::signal_create download_sig(_param.PSS_size);
                xvec signal = download_sig.download_convert_signal("real1.txt",_param.interp_coef);
                corelate_signal(signal);
                cell_id_calculate(signal);    

        }

        void create_mark::set_params()
        {
            param_store::parametrs params;//setting the selected parameters
            _param = params.set_param(_config);
            fpga_configure::config_fpga fpga(_param);//search PSS in the fpga correlator
            fpga.set_PSS_fpga();

        } //set_params

        void create_mark::corelate_signal(xvec signal)
        {
            //xcorr signal widh ref_signal
            LTE::xcorr corr;
            double lvl = 0;
            xvec PSS(_param.fftsize+_param.cp);
                for(int i = 0; i < 3; i++)
                {   
                    PSS = PSS_generate::pss::get_pss(i,_param.fftsize,_param.cp);
                    PSS_data = corr.correlate(signal,PSS,_param.windowing,_param.windowing,_param.shift);
                    if(PSS_data.maximum > lvl)
                    {
                        lvl = PSS_data.maximum;
                        _param.first = PSS_data.first;
                        _param.PSS_number = i;
                    }
                }

        } //corelate_signal

        
        void create_mark::cell_id_calculate(xvec signal)
        {
                SSS::SSS_PROC SSS;
                _param.CellId = SSS.CellId_calculate(signal,_param,_param.first);
        } // cell_id_calculate


} // namespace mark