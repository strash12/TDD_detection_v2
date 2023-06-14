#include "header.h"

namespace RS
{
    RS_cut::RS_cut(param_store::_parametrs param)
    {
        _param = param;
        d = 0;
    }
    
    RS_cut::~RS_cut(){};

    
    
    int RS_cut::set_RS(xvec sig_shift_1, xvec sig_shift_2, int first_1, int first_2)
    {   
        int sub = _param.fs/1e3;
        xvec Frame(sig_shift_1.begin()+first_1-(2*_param.fftsize)-_param.cp-_param.cp1-sub-1,sig_shift_1.begin()+first_1-(2*_param.fftsize)-_param.cp-_param.cp1-1+(4*sub));
        Frame.insert(Frame.end(),sig_shift_2.begin()+first_2-(2*_param.fftsize)-_param.cp-_param.cp1-sub-1,sig_shift_2.begin()+first_2-(2*_param.fftsize)-_param.cp-_param.cp1-1+(4*sub));
        gen_NRS();
        cut_slot(Frame);
        
	return search_maximum();
        
    }

    void RS_cut::gen_NRS()
    {
        try
        {
            rs_signal_fixed_real = new float*[2*_param.Nrb];
                for(int i = 0; i <2*_param.Nrb;i++)
                    {
                        rs_signal_fixed_real[i] = new float [_param.Ns*2];
                    }
     
            rs_signal_fixed_imag = new float*[2*_param.Nrb];
                for(int i = 0; i <2*_param.Nrb;i++)
                    {
                        rs_signal_fixed_imag[i] = new float [_param.Ns*2];
                    }
            
            rs_space_index = new int*[2*_param.Nrb];
                for(int i = 0; i <2*_param.Nrb;i++)
                    {
                        rs_space_index[i] = new int [_param.Ns*2];
                    }
        }
    catch (std::bad_alloc ba)
        {
            std::cout<<"no memory allocated"<<std::endl;
            std::cout<<ba.what()<<std::endl;
        }
        
        nrs::NRS_GEN RS(_param,_param.CellId);

       RS.NRC(_param,rs_signal_fixed_real,rs_signal_fixed_imag,rs_space_index);

    }

    void  RS_cut::cut_slot(xvec Frame)
    {
        int Slot_size = _param.fs/2000;
        LTE::FFT_transform fft(_param.fftsize);

        for (int snum = 0; snum<20;snum++)
        {
	        xvec SLot(Frame.begin()+snum*Slot_size,Frame.begin()+(snum+1)*Slot_size);
            xvec Slot_sig_1(SLot.begin()+_param.cp1,SLot.begin()+_param.cp1+_param.fftsize);
            xvec Slot_sig_2(SLot.begin()+_param.cp1 + 4*(_param.fftsize + _param.cp),SLot.begin()+_param.cp1 + 5*(_param.fftsize + _param.cp) - _param.cp);
            
            xvec SLOT_sig1_fft = fft.fft(Slot_sig_1);
            xvec SLOT_sig2_fft = fft.fft(Slot_sig_2);

            xvec slot_sig_1_cut(SLOT_sig1_fft.begin()+_param.fftsize/2-(2*_param.Nrb)*3,SLOT_sig1_fft.begin()+_param.fftsize/2);
                slot_sig_1_cut.insert(slot_sig_1_cut.end(),SLOT_sig1_fft.begin()+_param.fftsize/2+1,SLOT_sig1_fft.begin()+_param.fftsize/2+(2*_param.Nrb)*3+1);
                
                xvec slot_sig_2_cut(SLOT_sig2_fft.begin()+_param.fftsize/2-(2*_param.Nrb)*3,SLOT_sig2_fft.begin()+_param.fftsize/2);
                slot_sig_2_cut.insert(slot_sig_2_cut.end(),SLOT_sig2_fft.begin()+_param.fftsize/2+1,SLOT_sig2_fft.begin()+_param.fftsize/2+(2*_param.Nrb)*3+1);
            
	            
      gen_array(slot_sig_1_cut,slot_sig_2_cut);
        }
		
    }


    void RS_cut::gen_array(xvec Slot1,xvec Slot2)
    {
        xvec framesig_RS;
            framesig_RS.reserve(_param.Nrb*2);
                for (int i = 0; i < _param.Nrb*2;i++)
                {
                    int a = rs_space_index[i][0]-1;
                    framesig_RS.push_back(Slot1[a]);
                }

            xvec framesig_RS_2;
            framesig_RS_2.reserve(_param.Nrb*2);
                for (int i = 0; i < _param.Nrb*2;i++)
                {
                    int a = rs_space_index[i][1]-1;
                    framesig_RS_2.push_back(Slot2[a]);
                }
        
        fin_cor(framesig_RS);
        fin_cor(framesig_RS_2);
    }
    

    inline void RS_cut::fin_cor(xvec signal)
    {

        double CorRe = 0;
        double CorIm = 0;
        double Corlvln1 = 0 ;
        double maximum  = 0 ;    
            for (int i = 0; i<_param.Nrb*2 ; i++)
            {
                CorRe += signal[i].real()*rs_signal_fixed_real[i][d]-signal[i].imag()*(rs_signal_fixed_imag[i][d]);
                CorIm += signal[i].imag()*rs_signal_fixed_real[i][d]+signal[i].real()*(rs_signal_fixed_imag[i][d]);
            }    
        Corlvln1 = sqrt(pow(CorRe,2)+pow(CorIm,2));  
        cor_arr[d] = Corlvln1;
        CorRe = 0;
        CorIm = 0;
        Corlvln1 = 0;
    d++;
    
    }

int RS_cut::search_maximum()
    {
        
    std::cout<<"size1 = "<< 40<<std::endl;

            for(int i = 0; i<cor_arr.size(); i ++)
            {
                std::cout<<cor_arr[i]<<std::endl;
            } 
        
        
        
        double fincor[10]{{0}};
        
        for(int i = 0;i<10;i++)
        {
            fincor [i] = cor_arr[i*4]+cor_arr[i*4+1]+cor_arr[i*4+2]+cor_arr[i*4+3];
            std::cout<<" --"<<fincor [i]<<std::endl;
        }
 
        double maximume=0;

        for(int i=0;i<10;i++)
        {
        if (fincor[i]>maximume)
            {
                maximume = fincor[i];
            }
        }


        for (int i = 0;i<10;i++)
        {

            if (fincor[i]<(maximume/5))
            {
                TDD_conf[i] = 0;
            }

            else
            {
                TDD_conf[i] = 1;
            }

        }
        return TDD_config();
    }

    int RS_cut::TDD_config()
    {
    const std::array<std::array<int, 10>, 7> TDD_patterns {{
        {1,1,0,0,0,1,1,0,0,0},
        {1,1,0,0,1,1,1,0,0,1},
        {1,1,0,1,1,1,1,0,1,1},
        {1,1,0,0,0,1,1,1,1,1},
        {1,1,0,0,1,1,1,1,1,1},
        {1,1,0,1,1,1,1,1,1,1},
        {1,1,0,0,0,1,1,0,0,1}
    }};
    
    int TDD_config = -1;
    
    for (int i = 0; i < TDD_patterns.size(); ++i)
    {
        if (TDD_conf == TDD_patterns[i])
        {
            TDD_config = i;
            break;
        }
    }

    return TDD_config;
}

} //namespace RS