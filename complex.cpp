#include<iostream>
#include<string>
#include<complex.h>
#include<random>
#include<chrono>

namespace mComplex {
#define MAX(N1, N2, T1, T2) mComplex::max<(N1 > N2), N1, N2, T1, T2>
#define ENFORCE(x) typename = typename std::enable_if<(x)>::type
#define RIGHT_TYPE(N1, N2, T1, T2) MultiComplex<MAX(N1, N2, T1, T2)::value, typename MAX(N1, N2, T1, T2)::type>

template<class T>
using Complex = std::complex<T>;
template<unsigned N, class T>
class MultiComplex;

template<bool, unsigned N1, unsigned N2, class T1, class T2>
struct max {
    typedef T1 type;
    static constexpr unsigned value = N1;
};
template<unsigned N1, unsigned N2, class T1, class T2> 
struct max <true, N1, N2, T1, T2>{
    typedef T2 type;
    static constexpr unsigned value = N2;
};

//////////////////////////////////////////////








template<unsigned N, class T> // bicomplex as N = 2
struct MultiComplex { // only allow operations between same order MultiComplex

    MultiComplex<N-1, T> z1, z2;  
    MultiComplex<N-1, T>&real() { return z1; } 
    MultiComplex<N-1, T>&imag() { return z2; }
    MultiComplex<N-1, T> real() const { return z1; }
    MultiComplex<N-1, T> imag() const { return z2; }
    void printN() const { std::cout << N << std::endl; }

    MultiComplex<N, T> operator-() const { return MultiComplex<N, T>(-z1, -z2); }
    // overload operator unary +
    MultiComplex<N, T> operator+() const { return MultiComplex<N, T>(z1, z2); }
    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<N, T>& operator+=(const MultiComplex<N, T2>& w) { z1 += w.z1; z2 += w.z2;return *this; }
    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator+=(T2 w) { z1 += w; return *this; }

    template<class T2>
    MultiComplex<N, T>& operator+=(const Complex<T2>& w) { z1 += w; return *this; }

    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<N, T>& operator-=(const MultiComplex<N, T2>& w) { z1 -= w.z1; z2 -= w.z2; return *this; }

    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator-=(T2 w) { z1 -= w; return *this; }

    template<class T2>
    MultiComplex<N, T>& operator-=(const std::complex<T2>& w) { z1 -= w; return *this; }

    template<class T2>
    MultiComplex<N, T>& operator*=(const MultiComplex<N, T2>& w) {
        MultiComplex<N-1, T> tmp = z1 * w.z1 - z2 * w.z2;
        z2  = z1 * w.z2 + z2 * w.z1; z1 = tmp; return *this; }

    template<class T2>
    MultiComplex<N, T>& operator*=(const std::complex<T2>& w) { z1 *= w; z2 *= w; return *this; }

    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator*=(T2 w) { z1 *= w; z2 *= w; return *this; }

  

    // overload cout MultiComplex
    friend std::ostream& operator<<(std::ostream& os, const MultiComplex<N, T>& w) {
        os << "(" << w.z1 << " + " << w.z2 << "i" << N << ")";
        return os;
    }


};
// all additions
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator+(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z += w; return z; }
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator+(MultiComplex<N, T1> z, const std::complex<T2>& w) { z += w; return z; }
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator+(const std::complex<T2>& w, MultiComplex<N, T1> z) { z += w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator+(MultiComplex<N, T1> z, const T2& w) { z += w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator+(const T2& w, MultiComplex<N, T1> z) { z += w; return z; }
// all substractions
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator-(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z -= w; return z; }
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator-(MultiComplex<N, T1> z, const std::complex<T2>& w) { z -= w; return z; }
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator-(const std::complex<T2>& w, MultiComplex<N, T1> z) { z -= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator-(MultiComplex<N, T1> z, const T2& w) { z -= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator-(const T2& w, MultiComplex<N, T1> z) { z -= w; return z; }
// all multiplications
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator*(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z *= w; return z; }
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator*(MultiComplex<N, T1> z, const std::complex<T2>& w) { z *= w; return z; }
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator*(const std::complex<T2>& w, MultiComplex<N, T1> z) { z *= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator*(MultiComplex<N, T1> z, const T2& w) { z *= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator*(const T2& w, MultiComplex<N, T1> z) { z *= w; return z; }



template<class T>
struct MultiComplex<0, T> : public std::complex<T> { // special shit for bicomplex N = 0 (2^(1+N) complex numbers)
    using std::complex<T>::complex; 
};

template<unsigned N, class T>
void superImag(MultiComplex<N, T> &z, T *val) { 
   superImag(z.imag(), val);
}
template <class T>
void superImag(MultiComplex<0, T> &z, T *val) { 
    *val = z.imag();
}



template<unsigned N, class T>
void randomize(MultiComplex<N, T> &z) {
    
    randomize(z.real());
    randomize(z.imag());
}
template<class T>
void randomize(MultiComplex<0, T> &z) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dist(-1, 1);
    z = MultiComplex<0, T>(dist(gen), dist(gen));

    
}

};




int main() {
    using namespace mComplex;
    MultiComplex<1, long double> z1;
    MultiComplex<1, long double> z2;



    randomize(z1);
    randomize(z2);
    std::cout << z1 << std::endl;
    std::cout << z2 << std::endl;
    z1 *= z2;
    std::cout << z1 << std::endl;
 
    return 0;
}