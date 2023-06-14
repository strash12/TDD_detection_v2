#ifndef INCLUDE_LTE_PSS_H
#define INCLUDE_LTE_PSS_H

#include <cmath>
#include <fftw3.h>
#include <cstring>
#include <complex>
#include "../../TDD_detection/header.h"
namespace PSS_generate
{
class pss
{
public:
    pss();
    ~pss();

    static xvec get_pss(int Cell_id,int len,int cp);

private:    
    static const xd d_C_I;
    static const double d_PI;
    static void normalize(xd *zc_t,xd *pss_seq, int len);
    static void gen_pss_t(xd *zc_t,int Cell_id,int len);
    static void gen_conj_pss(xd *zc_t,int Cell_id,int len);
    static void zc(xd *zc, int Cell_id);

};
}

#endif /*INCLUDE_LTE_PSS_H*/