#ifndef defined SSS_PROCESS
#define SSS_PROCESS
#include "../../TDD_detection/header.h"

namespace SSS
{
    class SSS_PROC
    {
        public:
	     SSS_PROC();
            ~SSS_PROC();
            int CellId_calculate(xvec signal, param_store::_parametrs& param,int first);
        private:
            void cut_sss_signal(xvec signal,xvec& SSS_cut,param_store::_parametrs param,int first);
            void fft_sss_signal(int fft_size,xvec SSS_cut,xvec& SSS_fft);
            void cut_62_elemets_sss(int fft_size,xvec& SSS_62,xvec SSS_fft);
            int  sss_corelatate(double(&SSS)[504][62],double& maxSSS,xvec SSS_62);
    
    };//class SSS_PROC

}//namespace SSS





#endif //SSS_PROCESS