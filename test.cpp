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

    std::string to_string() const {
        return "(" + z1.to_string() + ")" + " + " + "(" + z2.to_string() + ")" + "i" + std::to_string(N+2) ;
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
struct MultiComplex<0, T>  { // special shit for bicomplex N = 0 (2^(1+N) complex numbers)
private:
    std::complex<T> z1, z2;
public:
    MultiComplex<0, T>() {} 
    MultiComplex<0, T>(std::complex<T> z1_, std::complex<T> z2_) : z1(z1_), z2(z2_) {}

    std::complex<T>&real() { return z1; } 
    std::complex<T>&imag() { return z2; }
    std::complex<T> real() const { return z1; }
    std::complex<T> imag() const { return z2; }
    void printN() const { std::cout << 0 << std::endl; }
   
    MultiComplex<0, T> operator-() const { return MultiComplex<0, T>(-z1, -z2); }
    // overload operator unary +
    MultiComplex<0, T> operator+() const { return MultiComplex<0, T>(z1, z2); }
    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<0, T>& operator+=(const MultiComplex<0, T2>& w) { z1 += w.z1; z2 += w.z2;return *this; }
    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<0, T>& operator+=(T2 w) { z1 += w; return *this; }

    template<class T2>
    MultiComplex<0, T>& operator+=(const Complex<T2>& w) { z1 += w; return *this; }

    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<0, T>& operator-=(const MultiComplex<0, T2>& w) { z1 -= w.z1; z2 -= w.z2; return *this; }

    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<0, T>& operator-=(T2 w) { z1 -= w; return *this; }

    template<class T2>
    MultiComplex<0, T>& operator-=(const std::complex<T2>& w) { z1 -= w; return *this; }

    template<class T2>
    MultiComplex<0, T>& operator*=(const MultiComplex<0, T2>& w) {
        std::complex<T> tmp = z1 * w.z1 - z2 * w.z2;
        z2  = z1 * w.z2 + z2 * w.z1; z1 = tmp; return *this; }

    template<class T2>
    MultiComplex<0, T>& operator*=(const std::complex<T2>& w) { z1 *= w; z2 *= w; return *this; }

    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<0, T>& operator*=(T2 w) { z1 *= w; z2 *= w; return *this; }

    std::string to_string() const {
        return "("+std::to_string(z1.real())+" + "+std::to_string(z1.imag()) +"i1)"+" + "+"("+std::to_string(z2.real())+" + "+std::to_string(z2.imag()) +"i1)i2";
    }
   
};
// all additions
template<class T1, class T2>
MultiComplex<0, T1> operator+(MultiComplex<0, T1> z, const MultiComplex<0, T2>& w) { z += w; return z; }
template<class T1, class T2>
MultiComplex<0, T1> operator+(MultiComplex<0, T1> z, const std::complex<T2>& w) { z += w; return z; }
template<class T1, class T2>
MultiComplex<0, T1> operator+(const std::complex<T2>& w, MultiComplex<0, T1> z) { z += w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator+(MultiComplex<0, T1> z, const T2& w) { z += w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator+(const T2& w, MultiComplex<0, T1> z) { z += w; return z; }
// all substractions
template<class T1, class T2>
MultiComplex<0, T1> operator-(MultiComplex<0, T1> z, const MultiComplex<0, T2>& w) { z -= w; return z; }
template<class T1, class T2>
MultiComplex<0, T1> operator-(MultiComplex<0, T1> z, const std::complex<T2>& w) { z -= w; return z; }
template<class T1, class T2>
MultiComplex<0, T1> operator-(const std::complex<T2>& w, MultiComplex<0, T1> z) { z -= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator-(MultiComplex<0, T1> z, const T2& w) { z -= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator-(const T2& w, MultiComplex<0, T1> z) { z -= w; return z; }
// all multiplications
template<class T1, class T2>
MultiComplex<0, T1> operator*(MultiComplex<0, T1> z, const MultiComplex<0, T2>& w) { z *= w; return z; }
template<class T1, class T2>
MultiComplex<0, T1> operator*(MultiComplex<0, T1> z, const std::complex<T2>& w) { z *= w; return z; }
template<class T1, class T2>
MultiComplex<0, T1> operator*(const std::complex<T2>& w, MultiComplex<0, T1> z) { z *= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator*(MultiComplex<0, T1> z, const T2& w) { z *= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator*(const T2& w, MultiComplex<0, T1> z) { z *= w; return z; }

template<unsigned N, class T>
void superImag(MultiComplex<N, T> &z, T *val) { 
   superImag(z.imag(), val);
}
template <class T>
void superImag(MultiComplex<0, T> &z, T *val) { 
    *val = z.imag().imag();
}

template<unsigned N, class T>
void setReal(MultiComplex<N, T> &z) { 
    setReal(z.real());
}
  
template<unsigned N, class T>
void setImag(MultiComplex<N, T> &z) { 
    setReal(z.imag());
}
template<class T>
void setReal(MultiComplex<0, T> &z) { 
    z.real() = std::complex<T>(1, 0);
}
  
template<class T>
void setImag(MultiComplex<0, T> &z) { 
    z.imag() = std::complex<T>(1, 0);
}

template<unsigned N, class T> 
void setUnit(MultiComplex<N, T> &z) {
    setImag(z);
    setUnit(z.real());
}
template<class T> 
void setUnit(MultiComplex<0, T> &z) {
    z.imag() = std::complex<T>(1, 0);
    z.real() = std::complex<T>(1, 1);
}

template<unsigned N, class T>
void randomize(MultiComplex<N, T> &z) {
    
    randomize(z.real());
    randomize(z.imag());
}
template<class T>
void randomize(MultiComplex<0, T> &z) {
    const double rand_max = 10, rand_min = -10;
    std::uniform_real_distribution<double> dist(rand_min, rand_max);
    std::default_random_engine gen;
    gen.seed(std::chrono::system_clock::now().time_since_epoch().count());
    z.real() = std::complex<T>(dist(gen), dist(gen));
    z.imag() = std::complex<T>(dist(gen), dist(gen));
    
}

};




int main() {
    using namespace mComplex;
    MultiComplex<6, float> z1;
    MultiComplex<6, float> z2;
    srand(0);
    randomize(z1); 
    randomize(z2);
    //std::cout << z1.to_string() << std::endl;
    // meassure loop execution time
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 100; i++) {
        z1 *= z2;
        z1 = z2;
    }
    auto end = std::chrono::high_resolution_clock::now();
    //std::cout << z1.to_string() << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << " micro s" << std::endl;
  
 
    return 0;
}