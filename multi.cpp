#include<iostream>
#include<string>
#include<complex.h>
#include<array>
#include<algorithm>

using Complex = std::complex<double>;
#define MAX(N1, N2) compare<(N1>N2), N1, N2>::ans


template <bool, unsigned N1, unsigned N2>
struct compare {
    static const unsigned ans = N1;
};
template <unsigned N1, unsigned N2>
struct compare<false, N1, N2> {
    static const unsigned ans = N2;
};






template<unsigned N>
class MultiComplex {
    public:
    uint n;
    MultiComplex<N-1> z1, z2;    
    MultiComplex<N>() {n = N;}
    void set_z1(Complex z) {
        z1.set_z1(z);
    }
    void set_z2(Complex z) {
        z2.set_z1(z);
    }

    
};
template <>
class MultiComplex<2> {
    public:
    uint n;
    Complex z1, z2;
    MultiComplex<2>() {n = 2;}
    void set_z1(Complex z) {
        z1 = z;
    }
    void set_z2(Complex z) {
        z2 = z;
    }
   
};

template <unsigned N1, unsigned N2>
MultiComplex<MAX(N1, N2)> test(MultiComplex<N1> z1, MultiComplex<N2> z2) {
    MultiComplex<MAX(N1, N2)> z;
    return z1;
}

int main() {
    

    MultiComplex<3> x;
    MultiComplex<2> y;
    auto z = test(x, y);
    std::cout << z.n << std::endl; 

    MultiComplex<4> x2;
    MultiComplex<3> y2;
    auto w = test(x2, y2);
    std::cout << w.n << std::endl; 

    
    return 0;
}