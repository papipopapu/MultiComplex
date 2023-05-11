#ifndef MULTICOMPLEX_HPP
#define MULTICOMPLEX_HPP


#include<iostream>
#include<random>



namespace MComplex {
#define ENFORCE(x) typename = typename std::enable_if<(x)>::type


template<unsigned N, class T> // C^N
struct MultiComplex { // only allow operations between same order MultiComplex

    MultiComplex<N-1, T> z1, z2;  
    MultiComplex<N-1, T>& real()       { return z1; }
    MultiComplex<N-1, T>& imag()       { return z2; }
    MultiComplex<N-1, T>  real() const { return z1; }
    MultiComplex<N-1, T>  imag() const { return z2; }
    void printN() const { std::cout << N << std::endl; }
    MultiComplex<N, T> operator-() const { return MultiComplex<N, T>{-z1, -z2}; }
    // overload operator unary +
    MultiComplex<N, T> operator+() const { return MultiComplex<N, T>{z1, z2}; }
    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<N, T>& operator+=(const MultiComplex<N, T2>& w) { z1 += w.z1; z2 += w.z2;return *this; }
    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator+=(T2 w) { z1 += w; return *this; }
    template<class T2>
    MultiComplex<N, T>& operator-=(const MultiComplex<N, T2>& w) { z1 -= w.z1; z2 -= w.z2; return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator-=(T2 w) { z1 -= w; return *this; }
    template<class T2>
    MultiComplex<N, T>& operator*=(const MultiComplex<N, T2>& w) {
        MultiComplex<N-1, T> tmp = z1 * w.z1 - z2 * w.z2;
        z2  = z1 * w.z2 + z2 * w.z1; z1 = tmp; return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator*=(T2 w) { z1 *= w; z2 *= w; return *this; }
    template<class T2>
    MultiComplex<N, T>& operator/=(const MultiComplex<N, T2>& w) {
        MultiComplex<N-1, T> den = w.z1 * w.z1 + w.z2 * w.z2;
        MultiComplex<N-1, T> tmp = z2 * w.z1 - z1 * w.z2;
        z1 = (z1 * w.z1 + z2 * w.z2) / den; 
        z2 = tmp                     / den;
        return *this; }

    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator/=(T2 w) { z1 /= w; z2 /= w; return *this; }
    // overload cout MultiComplex
    friend std::ostream& operator<<(std::ostream& os, const MultiComplex<N, T>& w) {
        os << "(" << w.z1 << " + " << w.z2 << "i" << N << ")";
        return os;
    }
    T& operator[](int i) {
        int K = N;
        int bruh = 2<<(K-1);
        if (i >= bruh) {
            return z2[i - bruh];
        } else {
            return z1[i];
        }
    } 
    T operator[](int i) const {
        int K = N;
        int bruh = 1<<(K-1);
        if (i >= bruh) {
            return z2[i - bruh];
        } else {
            return z1[i];
        }
    } 

};
// all additions
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator+(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z += w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator+(MultiComplex<N, T1> z, const T2& w) { z += w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator+(const T2& w, MultiComplex<N, T1> z) { z += w; return z; }
// all substractions
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator-(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z -= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator-(MultiComplex<N, T1> z, const T2& w) { z -= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator-(const T2& w, MultiComplex<N, T1> z) { z -= w; return z; }
// all multiplications
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator*(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z *= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator*(MultiComplex<N, T1> z, const T2& w) { z *= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator*(const T2& w, MultiComplex<N, T1> z) { z *= w; return z; }
// all divisions
template<unsigned N, class T1, class T2>
MultiComplex<N, T1> operator/(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z /= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator/(MultiComplex<N, T1> z, const T2& w) { z /= w; return z; }
template<unsigned N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator/(const T2& w, MultiComplex<N, T1> z) { z /= w; return z; }

template<unsigned N, class T1, class T2>
bool operator==(const MultiComplex<N, T1>& z, const MultiComplex<N, T2>& w) { return z.z1 == w.z1 && z.z2 == w.z2; }

template<unsigned N, class T1, class T2>
bool operator!=(const MultiComplex<N, T1>& z, const MultiComplex<N, T2>& w) { return !(z==w); }

template<class T>
struct MultiComplex<1, T> { // special shit for bicomplex N = 0 (2^(1+N) complex numbers)
    T z1, z2;
    T& real()       { return z1; }
    T& imag()       { return z2; }
    T  real() const { return z1; }
    T  imag() const { return z2; }
    void printN() const { std::cout << 0 << std::endl; }
    MultiComplex<1, T> operator-() const { return MultiComplex<1, T>{-z1, -z2}; }
    // overload operator unary +
    MultiComplex<1, T> operator+() const { return MultiComplex<1, T>{z1, z2}; }
    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<1, T>& operator+=(const MultiComplex<1, T2>& w) { z1 += w.z1; z2 += w.z2;return *this; }
    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<1, T>& operator+=(T2 w) { z1 += w; return *this; }
    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<1, T>& operator-=(const MultiComplex<1, T2>& w) { z1 -= w.z1; z2 -= w.z2; return *this; }
    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<1, T>& operator-=(T2 w) { z1 -= w; return *this; }
    template<class T2>
    MultiComplex<1, T>& operator*=(const MultiComplex<1, T2>& w) {
        T tmp = z1 * w.z1 - z2 * w.z2;
        z2  = z1 * w.z2 + z2 * w.z1; z1 = tmp; return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<1, T>& operator*=(T2 w) { z1 *= w; z2 *= w; return *this; }
    template<class T2>
    MultiComplex<1, T>& operator/=(const MultiComplex<1, T2>& w) {
        T den = w.z1 * w.z1 + w.z2 * w.z2;
        T tmp = z2 * w.z1 - z1 * w.z2;
        z1 = (z1 * w.z1 + z2 * w.z2) / den; 
        z2 = tmp                    / den;
        return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<1, T>& operator/=(T2 w) { z1 /= w; z2 /= w; return *this; }
    // overload cout MultiComplex
    friend std::ostream& operator<<(std::ostream& os, const MultiComplex<1, T>& w) {
        os << "(" << w.z1 << " + " << w.z2 << "i" << 1 << ")";
        return os;
    }
    T& operator[](int i) {
        if (i < 1) {
            return z1;
        } else {
            return z2;
        }
    } 
    T operator[](int i) const {
        if (i < 1) {
            return z1;
        } else {
            return z2;
        }
    } 
};
// all additions
template<class T1, class T2>
MultiComplex<1, T1> operator+(MultiComplex<1, T1> z, const MultiComplex<1, T2>& w) { z += w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<1, T1> operator+(MultiComplex<1, T1> z, const T2& w) { z += w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<1, T1> operator+(const T2& w, MultiComplex<1, T1> z) { z += w; return z; }
// all substractions
template<class T1, class T2>
MultiComplex<1, T1> operator-(MultiComplex<1, T1> z, const MultiComplex<1, T2>& w) { z -= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<1, T1> operator-(MultiComplex<1, T1> z, const T2& w) { z -= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<1, T1> operator-(const T2& w, MultiComplex<1, T1> z) { z -= w; return z; }
// all multiplications
template<class T1, class T2>
MultiComplex<1, T1> operator*(MultiComplex<1, T1> z, const MultiComplex<1, T2>& w) { z *= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<1, T1> operator*(MultiComplex<1, T1> z, const T2& w) { z *= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<1, T1> operator*(const T2& w, MultiComplex<1, T1> z) { z *= w; return z; }
// all divisions
template<class T1, class T2>
MultiComplex<1, T1> operator/(MultiComplex<1, T1> z, const MultiComplex<1, T2>& w) { z /= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<1, T1> operator/(MultiComplex<1, T1> z, const T2& w) { z /= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<1, T1> operator/(const T2& w, MultiComplex<1, T1> z) { z /= w; return z; }

template<class T1, class T2>
bool operator==(const MultiComplex<1, T1>& z, const MultiComplex<1, T2>& w) { return z.z1 == w.z1 && z.z2 == w.z2; }

template<class T1, class T2>
bool operator!=(const MultiComplex<1, T1>& z, const MultiComplex<1, T2>& w) { return !(z==w); }


template<unsigned N, class T>
void superImag(MultiComplex<N, T> &z, T *val) { 
   superImag(z.imag(), val);
}
template <class T>
void superImag(MultiComplex<1, T> &z, T *val) { 
    *val = z.imag();
}

template<unsigned N, class T>
void randomize(MultiComplex<N, T> &z) {
    
    randomize(z.real());
    randomize(z.imag());
}
template<class T>
void randomize(MultiComplex<1, T> &z) {
    // create random MultiComplex with deterministic seed
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dist(-1, 1);
    z = MultiComplex<1, T>{dist(gen), dist(gen)};
}

template<unsigned N, class T>
void toOne(MultiComplex<N, T> &z) {
    toOne(z.real());
    toOne(z.imag());
}
template<class T>
void toOne(MultiComplex<1, T> &z) {
    z = MultiComplex<1, T>{1, 0};
}

template<unsigned N, class T, ENFORCE(N==1)>
MultiComplex<1, T> I() {
    return MultiComplex<1, T>{0, 1};
}
template<unsigned N, class T, ENFORCE(N==1)>
MultiComplex<1, T> R() {
    return MultiComplex<1, T>{1, 0};
}

template<unsigned N, class T, ENFORCE(N>1)>
MultiComplex<N, T> R() {
    return MultiComplex<N, T>{R<N-1, T>(), MultiComplex<N-1, T>()};
}
template<unsigned N, class T, ENFORCE(N>1)>
MultiComplex<N, T> I() {
    return MultiComplex<N, T>{MultiComplex<N-1, T>(), R<N-1, T>()};
}
template<unsigned N1, unsigned N2, class T, ENFORCE(N1==N2)>
MultiComplex<N2, T> promote(const MultiComplex<N2, T>& z) {
    return z;
}
template<unsigned N1, unsigned N2, class T, ENFORCE(N1>N2)>
MultiComplex<N1, T> promote(const MultiComplex<N2, T>& z) {
    return MultiComplex<N1, T>{promote<N1-1, N2, T>(z), MultiComplex<N1-1, T>()};
}
template<unsigned N, class T>
MultiComplex<N, T> promote(const T& val) {
    MultiComplex<N, T> ret;
    ret[0] = val;
    return ret;
}

template<unsigned N, class T>
MultiComplex<N, T> ImagUnit() {
    MultiComplex<N, T> ret{};
    for (int i = 0; i < N; ++i) {
        ret[1<<i] = 1;
    }
    return ret;
}
// Actual functions

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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

#endif


