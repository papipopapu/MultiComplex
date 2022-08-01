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
public:
    MultiComplex<N>() {} 
    MultiComplex<N>(TYPE(N-1) z1_, TYPE(N-1) z2_) : z1(z1_), z2(z2_) {}   

    TYPE(N-1)&real() { return z1; } 
    TYPE(N-1)&imag() { return z2; }
    TYPE(N-1) real() const { return z1; }
    TYPE(N-1) imag() const { return z2; }
    void printN() const { std::cout << N << std::endl; }

    
    // operators
    // overload operator unary -
    MultiComplex<N> operator-() const {
        return MultiComplex<N>(-z1, -z2);
    }
    // overload operator unary +
    MultiComplex<N> operator+() const {
        return MultiComplex<N>(z1, z2);
    }

    // overload operator +
        // MultiComplex += MultiComplex
    MultiComplex<N>& operator+=(const MultiComplex<N>& w) {
        z1 += w.z1;
        z2 += w.z2;
        return *this;
    }
        // MultiComplex += Complex
    MultiComplex<N>& operator+=(const Complex& w) {
        z1 += w;
        return *this;
    }
    // overload operator -
        // MultiComplex += MultiComplex
    MultiComplex<N>& operator-=(const MultiComplex<N>& w) {
        z1 -= w.z1;
        z2 -= w.z2;
        return *this;
    }
        // MultiComplex += Complex
    MultiComplex<N>& operator-=(const Complex& w) {
        z1 -= w;
        return *this;
    }
    // overload operator *
        // MultiComplex *= MultiComplex
    MultiComplex<N>& operator*=(const MultiComplex<N>& w) {
        TYPE(N-1) tmp = z1 * w.z1 - z2 * w.z2;
        z2  = z1 * w.z2 + z2 * w.z1;
        z1 = tmp;
        return *this;
    }
        // MultiComplex * Complex
    MultiComplex<N>& operator*=(const Complex& w) {
        z1 *= w;
        z2 *= w;
        return *this;
    }

    std::string to_string() const {
        return "(" + std::to_string(z1.real()) + " + " + std::to_string(z2.real()) + "i" + std::to_string(N) + ")";
    }
};
template<unsigned N>
void setReal(MultiComplex<N> &z, Complex w) { 
        setImag(z.real(), w);
}
  
template<unsigned N>
void setImag(MultiComplex<N> &z, Complex w) { 
        setReal(z.imag(), w);
}
template<>
void setReal(MultiComplex<1> &z, Complex w) { 
        z.real() = w;
}
  
template<>
void setImag(MultiComplex<1> &z, Complex w) { 
        z.imag() = w;
}

template<unsigned N>
void superImag(MultiComplex<N> &z, double *val) { 
        superImag(z.imag(), val);
}
template <>
void superImag(MultiComplex<1> &z, double *val) { 
        *val = z.imag().imag();
}
  
int main() {
    

    
    MultiComplex<2> x;
    setReal(x, Complex(1, 2));
    setImag(x, Complex(0, 1));

    MultiComplex<1> y;
    setReal(y, Complex(2, 3));
    setImag(y, Complex(3, 2));

    x.imag() = y;

    double ans;
    superImag(x, &ans);
    std::cout << x.imag().real() << " " << ans << std::endl;



    
    return 0;
}