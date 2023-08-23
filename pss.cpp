#include "header.h"

namespace PSS_generate
{
pss::pss(){}
pss::~pss(){}


int huuui = 0;
const xd
pss::d_C_I = xd(0,1);

const double
pss::d_PI = double(M_PI);

xvec pss::get_pss(int Cell_id,int len,int cp)
{
    if (Cell_id > 2)
        throw "invalid PSS number";
    xvec zc_t(len+cp);
    xd pss_seq[len];
    xd pss_seq_normalize[len];
    gen_pss_t(pss_seq,Cell_id,len);
    normalize(pss_seq_normalize,pss_seq,len);

    std::copy(pss_seq_normalize+len-cp+1,pss_seq_normalize+len,zc_t.begin());
    std::copy(pss_seq_normalize,pss_seq_normalize+len,zc_t.begin()+cp-1);
    return (zc_t);



    //memcpy(zc_t,pss_seq_normalize+len-cp+1,sizeof(xd)*cp);
    //memcpy(zc_t+cp-1,pss_seq_normalize,sizeof(xd)*len);

}

void pss::gen_pss_t(xd *zc_t,int Cell_id,int len)
{
    xd zc_f[62];
    zc(zc_f,Cell_id);

    xd* d_in  = (xd*) fftw_malloc(sizeof(xd)*len);
    xd* d_out  = (xd*) fftw_malloc(sizeof(xd)*len);

    memset(d_in,0,sizeof(xd)*len);
    memcpy(d_in+len-31, zc_f, sizeof(xd)*31);
    memcpy(d_in+1, zc_f+31, sizeof(xd)*31);
                  

    fftw_plan p = fftw_plan_dft_1d(len,reinterpret_cast<fftw_complex*>(d_in),reinterpret_cast<fftw_complex*>(d_out),FFTW_BACKWARD,FFTW_ESTIMATE);
    
    fftw_execute(p);

    memcpy(zc_t, d_out, sizeof(xd)*len);

    fftw_destroy_plan(p);
    fftw_free(d_in);
    fftw_free(d_out);
    
}

void pss::zc(xd *zc, int cell_id)
{
    
    float u=0;
    
    switch (cell_id)
    {
    case 0:
        u=25.0;
        break;
    case 1:
        u=29.0;
        break;
    case 2:
        u=34.0;
        break;
    }

   
    for(int n = 0; n < 31; n++)
    {
        zc[n]=exp(d_C_I* xd(d_PI*u* double(-1*n*(n+1))/63.0 ) );
    }
    for(int n = 31; n < 62; n++)
    {
        zc[n]=exp(d_C_I* xd(d_PI*u* double(-1*(n+1)*(n+2))/63.0 ) );
    }

}

void pss::normalize(xd *zc_t,xd *pss_seq, int len)
{
    xd pss_normalize[len];
    for(int i = 0; i<len;i++){
       pss_normalize[i] = xd(pss_seq[i].real()/len,pss_seq[i].imag()/len);
   }

   memcpy(zc_t, pss_normalize, sizeof(xd)*len);
}

}
