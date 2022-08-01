#include<iostream>
#include<string>
#include<complex.h>
#include<array>
#include<algorithm>

namespace mComplex {
template<class T>
using Complex = std::complex<T>;


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

    
    
};
template<class T>
class MultiComplex<0U, T> : public std::complex<T> { 
public:
    using Complex<T>::Complex;
};


};
int main() {
    using namespace mComplex;


    MultiComplex<0U, float> z0(1.0, 0.0);
    std::cout << std::exp(z0) << std::endl;
  
 
    return 0;
}