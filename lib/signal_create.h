#ifndef defined SIGNAL_CREATE
#define SIGNAL_CREATE

#include<complex>
#include<vector>

typedef std::complex<double> xd;
typedef std::vector<xd> xvec;

namespace signal
{
    class signal_create
    {
        private:
            xvec convert_signal_normal(int,int);
            xvec convert_signal_extended(int,int);
            void save_signal_to_file(std::vector<int32_t> samples,std::string name);
            int size;
        public:
            signal_create(int PSS_signal_size);
            ~signal_create();
            xvec PSS(int,int,int);
            xvec download_convert_signal(std::string name,int coef);
            xvec PSS_signal;
            
    };
} //namespace signal


#endif //SIGNAL_CREATE