#include"multicomplex.cpp"
namespace MComplex {
    static constexpr int EXPANSION_TERMS = 10;
    static constexpr int HALLEY_ITER = 5;


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
        std::cout << "//// exp ("<<N<< ")  args{"<< z <<"}" << std::endl;  
        MultiComplex<N-1, T> r  = exp(z.real());
        MultiComplex<N, T> ret{r * cos(z.imag()), r * sin(z.imag())};  
        std::cout << "exp(" << N << "): " << ret << std::endl;  
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
};

