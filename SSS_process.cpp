#include "header.h"
#include "../TDD_detection/signal_storage/SSS.h"
namespace SSS
{
    SSS::SSS_PROC::SSS_PROC(){}
    SSS::SSS_PROC::~SSS_PROC(){}

    int SSS::SSS_PROC::CellId_calculate(xvec signal, param_store::_parametrs& param,int first)
    {
        xvec SSS_cut(param.fftsize);
        cut_sss_signal(signal,SSS_cut,param,first);
        xvec SSS_fft(param.fftsize);
        fft_sss_signal(param.fftsize,SSS_cut,SSS_fft);
        xvec SSS_62(62);
        cut_62_elemets_sss(param.fftsize,SSS_62,SSS_fft);
        double maxSSS1 = 0,maxSSS2 = 0;
        int Cellid1 = sss_corelatate(SSS_first_part,maxSSS1,SSS_62);
        int Cellid2 = sss_corelatate(SSS_second_part,maxSSS2,SSS_62);
        int CellId = maxSSS1>maxSSS2?Cellid1:Cellid2;
        std::cout<<"max 1 = "<<maxSSS1<< "cel 1  = "<< Cellid1<<std::endl;
        std::cout<<"max 2 = "<<maxSSS2<< "cel 2  = "<< Cellid2<<std::endl;
        std::cout<<"cel   = "<< CellId<<std::endl;
        return CellId;

    }

    void SSS::SSS_PROC::cut_sss_signal(xvec signal,xvec& SSS_cut,param_store::_parametrs param,int first)
    {
        int start = first - (param.fftsize + param.cp)- (param.fftsize + param.cp1) - param.fftsize-1;
        int stop  = first - 2 * param.fftsize -  param.cp-  param.cp1 -1;
        std::copy(signal.begin()+start,signal.begin()+stop,SSS_cut.begin());
    }

    void SSS::SSS_PROC::fft_sss_signal(int fft_size,xvec SSS_cut,xvec& SSS_fft)
    {
        LTE::FFT_transform fft(fft_size);
        SSS_fft = fft.fft(SSS_cut);
    }

    void SSS::SSS_PROC::cut_62_elemets_sss(int fft_size,xvec& SSS_62,xvec SSS_fft)
    {
        int begin_1 = fft_size/2 - 31;
        int end_1   = fft_size/2;
        int begin_2 = fft_size/2 + 1;
        int end_2   = fft_size/2 + 32; 

            std::copy(SSS_fft.begin()+begin_1,SSS_fft.begin()+end_1,SSS_62.begin());
            std::copy(SSS_fft.begin()+begin_2,SSS_fft.begin()+end_2,SSS_62.begin()+31);

    }

    int SSS::SSS_PROC::sss_corelatate(double (&SSS_signal)[504][62],double& maxSSS,xvec SSS_62)
    {   
    //double CorRe = 0,CorIm = 0;
    double Corlvln = 0;
    int Cellid = 0;       
    xd CorRe;
        for ( int i = 0; i < 504;i++)
            {
                for (int j = 0; j < 62; j++)
                    {
                        CorRe  +=  SSS_62[j]*SSS_signal[i][j];      
                    }

           Corlvln =  sqrt(pow(CorRe.real(),2)+pow(CorRe.imag(),2));
                    if(Corlvln > maxSSS)
                        {
                            maxSSS = Corlvln;
                            Cellid = i;
                        }
            CorRe  = 0;
            }
    return Cellid;
    }

}