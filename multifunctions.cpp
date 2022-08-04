#include"multicomplex.cpp"
namespace MComplex {
    static constexpr int EXPANSION_TERMS = 10;
    static constexpr int HALLEY_ITER = 5;


    template<class T>
    T sin(const T& val) { return std::sin(val); }
    template<class T>
    T cos(const T& val) { return std::cos(val); }
    template<class T>
    T tan(const T& val) { return std::tan(val); }
    template<class T>
    T sinh(const T& val) { return std::sinh(val); }
    template<class T>
    T cosh(const T& val) { return std::cosh(val); }
    template<class T>
    T tanh(const T& val) { return std::tanh(val); }
    template<class T>
    T asin(const T& val) { return std::asin(val); }
    template<class T>
    T acos(const T& val) { return std::acos(val); }
    template<class T>
    T atan(const T& val) { return std::atan(val); }
    template<class T>
    T asinh(const T& val) { return std::asinh(val); }
    template<class T>
    T acosh(const T& val) { return std::acosh(val); }
    template<class T>
    T atanh(const T& val) { return std::atanh(val); }

    
    template<class T>
    T log(const T& val) { return std::log(val); }
    template<class T>
    T exp(const T& val) { return std::exp(val); }




    template<unsigned N, class T>
    MultiComplex<N, T> tan(const MultiComplex<N, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> sin(const MultiComplex<N, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> cos(const MultiComplex<N, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> sinh(const MultiComplex<N, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> cosh(const MultiComplex<N, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> exp(const MultiComplex<N, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> log(const MultiComplex<N, T>& z);
    template<unsigned N, class T>
    MultiComplex<N-1, T> arg(const MultiComplex<N, T>& z);
    template<class T>
    T arg(const MultiComplex<0, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const MultiComplex<N, T>& w);
    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const T& w);
    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const int& n);
    template<unsigned N, class T>
    MultiComplex<N-1, T> abs(const MultiComplex<N, T>& z);
    template<class T>
    T abs(const MultiComplex<0, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> atan(const MultiComplex<N, T>& z);
    




    template<unsigned N, class T>
    MultiComplex<N, T> tan(const MultiComplex<N, T>& z) {
        MultiComplex<N, T> ret = sin(z)/cos(z);      
        return ret;
    }

    template<unsigned N, class T>
    MultiComplex<N, T> sin(const MultiComplex<N, T>& z) {
        MultiComplex<N, T> ret{sin(z.real()) * cosh(z.imag()), cos(z.real()) * sinh(z.imag())};        
        return ret;
    }
    template<unsigned N, class T>
    MultiComplex<N, T> cos(const MultiComplex<N, T>& z) {
        MultiComplex<N, T> ret{cos(z.real()) * cosh(z.imag()), sin(z.real()) * sinh(z.imag())};        
        return ret;
    }
    template<unsigned N, class T>
    MultiComplex<N, T> sinh(const MultiComplex<N, T>& z) {
        MultiComplex<N, T> ret{sinh(z.real()) * cos(z.imag()), cosh(z.real()) * sin(z.imag())};        
        return ret;
    }
    template<unsigned N, class T>
    MultiComplex<N, T> cosh(const MultiComplex<N, T>& z) {
        MultiComplex<N, T> ret{cosh(z.real()) * cos(z.imag()), sinh(z.real()) * sin(z.imag())};        
        return ret;
    }
    template<unsigned N, class T>
    MultiComplex<N, T> exp(const MultiComplex<N, T>& z) {
        auto r  = exp(z.real());
        MultiComplex<N, T> ret{r * cos(z.imag()), r * sin(z.imag())};  
        return ret;
    }
    template<unsigned N, class T>
    MultiComplex<N, T> log(const MultiComplex<N, T>& z) {  
        MultiComplex<N, T> y = promote<N, T>(log(z[0]));
        MultiComplex<N, T> gamma;
        for(unsigned i = 1; i < HALLEY_ITER; i++) {
            gamma = z * exp(-y);
            y -= 2 * (1 - gamma) / (1 + gamma);
        }
        return y;
    }
    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const int& n) {
        MultiComplex<N, T> ret;
        if (n == 0) return ret;
        ret = z;
        for (int i = 1; i < n; i++) {
            ret *= z;
        }
        return ret;
    }
    template<unsigned N, class T>
    MultiComplex<N, T> atan(const MultiComplex<N, T>& z) {
        const int terms = 5;
        MultiComplex<N, T> ret;
        for (int n = 0; n < EXPANSION_TERMS; n++) {
            ret += ((-1)^n) * pow(z, 2*n+1) / (2*n+1);
        }
        return ret;
    }
};

