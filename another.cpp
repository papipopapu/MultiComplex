#include<iostream>
#include<string>
#include<complex.h>
#include<array>
#include<algorithm>

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


template<unsigned N1, unsigned N2, class T1, class T2, ENFORCE(N1 != N2)>
void plus_equals(MultiComplex<N1, T1> &a, const MultiComplex<N1, T1> &b) {
    a.real() += b;
}
template<unsigned N1, unsigned N2, class T1, class T2, ENFORCE(N1 == N2)>
void plus_equals(MultiComplex<N1, T1> &a, const MultiComplex<N2, T2> &b) {
    a.imag() += b.imag();
    a.real() += b.real();
}
template<unsigned N1, unsigned N2, class T1, class T2, ENFORCE(N1 >= N2)>
RIGHT_TYPE(N1, N2, T1, T2) plus(MultiComplex<N1, T1> a, const MultiComplex<N2, T2> &b) {
    a += b;
    return a;
}
template<unsigned N1, unsigned N2, class T1, class T2, ENFORCE(N1 < N2)>
RIGHT_TYPE(N1, N2, T1, T2) plus(const MultiComplex<N1, T1>& a, MultiComplex<N2, T2> b) {
    b += a;
    return b;
}





template<unsigned N, class T> // bicomplex as N = 1
class MultiComplex { // only allow operations between same order MultiComplex
private:
    MultiComplex<N-1, T> z1, z2;  
public:
    MultiComplex<N, T>() {} 
    MultiComplex<N, T>(MultiComplex<N-1, T> z1_, MultiComplex<N-1, T> z2_) : z1(z1_), z2(z2_) {}   

    MultiComplex<N-1, T>&real() { return z1; } 
    MultiComplex<N-1, T>&imag() { return z2; }
    MultiComplex<N-1, T> real() const { return z1; }
    MultiComplex<N-1, T> imag() const { return z2; }
    void printN() const { std::cout << N << std::endl; }
    
    // interact with lower order MultiComplex
    template<unsigned N2, class T2, ENFORCE(N>N2)>
    MultiComplex<N, T>& operator+=(const MultiComplex<N2, T2>& w) {
        plus_equals<N, N2, T, T2>(*this, w);
        return *this;
    }
    // we only interact with other types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator+=(T2 w) {
        z1 += w;
        return *this;
    }
    template<class T2>
    MultiComplex<N, T>& operator+=(const Complex<T2>& w) {
        z1 += w;
        return *this;
    }

};
template<unsigned N1, unsigned N2, class T1, class T2>
RIGHT_TYPE(N1, N2, T1, T2) operator+(MultiComplex<N1, T1> z, MultiComplex<N2, T2> w) {
    return plus<N1, N2, T1, T2>(z, w);
}

template<class T>
class MultiComplex<0U, T> : public std::complex<T> { 
public:
    using Complex<T>::Complex;
};




};
int main() {
    using namespace mComplex;


    MultiComplex<0U, float> z0(1.0, 0.0);
    MultiComplex<1U, float> z1(z0, z0);
    auto z3 = z1 + z1;

    z1 += 1;
    
    std::cout << z1.real() << std::endl;
  
 
    return 0;
}