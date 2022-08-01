#include<iostream>
#include<string>
#include<complex.h>
#include<array>
#include<algorithm>

using Complex = std::complex<double>;

template<unsigned N>
class MultiComplex;

template<unsigned N> 
struct compare {
    typedef MultiComplex<N> type;
};

template<>
struct compare<0U> {
    typedef Complex type;
};

#define TYPE(N) typename compare<N>::type


template<unsigned N> // bicomplex as N = 1
class MultiComplex { // only allow operations between same order MultiComplex
private:
    TYPE(N-1) z1, z2;  
    uint n;
public:
    MultiComplex<N>() : n(N) {} 
    MultiComplex<N>(TYPE(N-1) z1_, TYPE(N-1) z2_) : n(N), z1(z1_), z2(z2_) {}   
    MultiComplex<N>& operator=(MultiComplex<N> w) {
        std::swap(w);
        return *this;
    }
    void setRe(Complex z) { 
        auto next = this;
        for (uint i = 1; i < N; i++) {
            auto tmp = next;
            auto next = tmp->Re();
        }
        next->Re() = z;
    }
    void setIm(Complex z) { 
        auto next = this;
        for (uint i = 1; i < N; i++) {
            auto tmp = next;
            if (i == 1) {
                auto next = tmp->Re();
            } else {
            auto next = tmp->Re();
            }
        }
        if (N == 1) {
            next->Im() = z;
        } else {
            next->Re() = z;
        }
    }
    
    TYPE(N-1)&Re() { return z1; } 
    TYPE(N-1)&Im() { return z2; }
    TYPE(N-1) Re() const { return z1; }
    TYPE(N-1) Im() const { return z2; }
};


int main() {
    

    
    MultiComplex<1> x;
    x.setRe(Complex(1, 2));
    x.setIm(Complex(0, 1));

    MultiComplex<1> y;
    y.setRe(Complex(2, 3));
    y.setIm(Complex(3, 2));


    std::cout << x.Re() << std::endl;



    
    return 0;
}