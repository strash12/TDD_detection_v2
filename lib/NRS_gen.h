#ifndef defined NRS_gen
#define NRS_gen

#include "../../TDD_detection/header.h"
#include <bitset>


namespace nrs
{
class NRS_GEN
{
    protected:
        int M;
        int Ncp,Ns_fix,fix_1,C_init,C_init_fix,r,m_,Lb;
        int Ns = 20;
        int _CellId;
        int rs_space_1_fix_normal[40];
        int rs_space_1_fix_extended[40];
        int NS_1_normal[2];
        int NS_1_extended[2];

    public:
        NRS_GEN(param_store::_parametrs param,int);
        ~NRS_GEN();
        void NRC(param_store::_parametrs param,float** rs_signal_fixed_real, float** rs_signal_fixed_imag, int** rs_space_index );
};


} //namespace RS


#endif //NRS_gen