#include"multicomplex.cpp"
# define M_PI 3.14159265358979323846
namespace MComplex {



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
    T atan2(const T& val) { return std::atan2(val); } 
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
    template<class T>
    T pow(const T& val, const T& exp) { return std::pow(val, exp); }
    template<class T>
    T abs(const T& val) { return std::abs(val); }
  




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
    T arg(const MultiComplex<1, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const MultiComplex<N, T>& w);
    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const T& w);
    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const int& n);
    template<unsigned N, class T>
    MultiComplex<N-1, T> abs(const MultiComplex<N, T>& z);
    template<class T>
    T abs(const MultiComplex<1, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> atan(const MultiComplex<N, T>& z);
    template<unsigned N, class T>
    MultiComplex<N, T> atan2(const MultiComplex<N, T>& z, const MultiComplex<N, T>& w);
    

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
        MultiComplex<N, T> ret{cos(z.real()) * cosh(z.imag()), -sin(z.real()) * sinh(z.imag())};        
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
        return MultiComplex<N, T>{log(abs(z)), arg(z)};
    }

    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const T& val) {
        return exp(val * log(z));
    }
    template<unsigned N, class T>
    MultiComplex<N, T> pow(const MultiComplex<N, T>& z, const int& n) {
        if (n == 0) return MultiComplex<N, T>{1};
        MultiComplex<N, T> ret = z;
        for (int i = 1; i < n; i++) ret *= z;
        return ret;
    }
    template<unsigned N, class T>
    MultiComplex<N, T> atan(const MultiComplex<N, T>& z) {
        return I<N, T>() * 0.5 * log((I<N, T>() + z) / (I<N, T>() - z));
    }
    template<unsigned N, class T>
    MultiComplex<N, T> atan2(const MultiComplex<N, T>& z,const MultiComplex<N, T>& w) { // z/w
        T y = z[0], x = w[0]; // only make decision on the normal complex branch cuts, so log(exp(x)) only works for complex, etc.
        if (x > 0) {
            return 2 * atan( z / (pow(z * z + w * w, 0.5) + w) );
        } else if (x <= 0 && y != 0) {
            return 2 * atan(     (pow(z * z + w * w, 0.5) - w) / z );
        } else if (x <  0 && y == 0) {
            return promote<N, T>(M_PI);
        } else  { // (x == 0 && y == 0)
            return MultiComplex<N, T>{};
        }
    }
    template<unsigned N, class T>
    MultiComplex<N-1, T> abs(const MultiComplex<N, T>& z) {
        return pow(z.real() * z.real() + z.imag() * z.imag(), 0.5);
    }
    template<class T>
    T abs(const MultiComplex<1, T>& z) {
        return pow(z.real() * z.real() + z.imag() * z.imag(), 0.5);
    }
    template<unsigned N, class T>
    MultiComplex<N-1, T> arg(const MultiComplex<N, T>& z) {
        return atan2(z.imag(), z.real());
    }
    template<class T>
    T arg(const MultiComplex<1, T>& z) {
        return std::atan2(z.imag(), z.real());
    }

};

