#ifndef defined LTE_PROC
#define LTE_PROC

#include "../../TDD_detection/header.h"
#include <string>



const xd J(0, 1);

namespace LTE
{
    struct cor_data
    {
        double maximum ;
        double level   ;
        int     first   ;
        int     number  ;
        xd cor_lvl;
        
    };//cor_data

    class xcorr
    {
        public:
        xcorr();
        ~xcorr();
        cor_data correlate(xvec signal,xvec ref_signal,int start,int stop,int shift);
        int corr_coef(xvec signal,xvec ref_signal,int ref_samples,int fftsize,int cp,int shift_second_cor_start, int shift_second_cor_stop);
        private:
        
    };//class LTE
        
    class FFT_transform
    {
        private:
            int _FFT_SIZE;
            //xd omega(const double &p, const double &q);
            //xvec _fft(xvec &f);
        public:
            FFT_transform(int fftsize);
            ~FFT_transform();
            xvec fft(xvec x);
           
    };//class FFT

    class freq_shift_estimate
    {
        public:
            freq_shift_estimate();
            ~freq_shift_estimate();
            double get_angle(xvec signal,int cp,int fftsize,int first);
        private:
            double calculate_angle(xvec cp1,xvec cp2,int cp);
    };//class freq_shift_estimate

} //namespace LTE


#endif //LTE_PROC