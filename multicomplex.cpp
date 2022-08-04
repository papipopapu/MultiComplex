#include<iostream>
#include<string>
#include<random>
#include<chrono>

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
        z2 = tmp                    / den;
        return *this; }

    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator/=(T2 w) { z1 /= w; z2 /= w; return *this; }
    // overload cout MultiComplex
    friend std::ostream& operator<<(std::ostream& os, const MultiComplex<N, T>& w) {
        os << "(" << w.z1 << " + " << w.z2 << "i" << N << ")";
        return os;
    }
    T& operator[](const unsigned& i) {
        if (i > N) {
            return z2[i - N];
        } else {
            return z1[i];
        }
    } 
    T operator[](const unsigned& i) const {
        if (i > N) {
            return z2[i - N];
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
struct MultiComplex<0, T> { // special shit for bicomplex N = 0 (2^(1+N) complex numbers)
    T z1, z2;
    T& real()       { return z1; }
    T& imag()       { return z2; }
    T  real() const { return z1; }
    T  imag() const { return z2; }
    void printN() const { std::cout << 0 << std::endl; }
    MultiComplex<0, T> operator-() const { return MultiComplex<0, T>{-z1, -z2}; }
    // overload operator unary +
    MultiComplex<0, T> operator+() const { return MultiComplex<0, T>{z1, z2}; }
    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<0, T>& operator+=(const MultiComplex<0, T2>& w) { z1 += w.z1; z2 += w.z2;return *this; }
    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<0, T>& operator+=(T2 w) { z1 += w; return *this; }
    // interact with same order MultiComplex
    template<class T2>
    MultiComplex<0, T>& operator-=(const MultiComplex<0, T2>& w) { z1 -= w.z1; z2 -= w.z2; return *this; }
    // we only interact with other (lower) types if they are Complex or Real
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<0, T>& operator-=(T2 w) { z1 -= w; return *this; }
    template<class T2>
    MultiComplex<0, T>& operator*=(const MultiComplex<0, T2>& w) {
        T tmp = z1 * w.z1 - z2 * w.z2;
        z2  = z1 * w.z2 + z2 * w.z1; z1 = tmp; return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<0, T>& operator*=(T2 w) { z1 *= w; z2 *= w; return *this; }
    template<class T2>
    MultiComplex<0, T>& operator/=(const MultiComplex<0, T2>& w) {
        T den = w.z1 * w.z1 + w.z2 * w.z2;
        T tmp = z2 * w.z1 - z1 * w.z2;
        z1 = (z1 * w.z1 + z2 * w.z2) / den; 
        z2 = tmp                    / den;
        return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<0, T>& operator/=(T2 w) { z1 /= w; z2 /= w; return *this; }
    // overload cout MultiComplex
    friend std::ostream& operator<<(std::ostream& os, const MultiComplex<0, T>& w) {
        os << "(" << w.z1 << " + " << w.z2 << "i" << 0 << ")";
        return os;
    }
    T& operator[](const unsigned& i) {
        if (i < 1) {
            return z1;
        } else {
            return z2;
        }
    } 
    T operator[](const unsigned& i) const {
        if (i < 1) {
            return z1;
        } else {
            return z2;
        }
    } 
};
// all additions
template<class T1, class T2>
MultiComplex<0, T1> operator+(MultiComplex<0, T1> z, const MultiComplex<0, T2>& w) { z += w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator+(MultiComplex<0, T1> z, const T2& w) { z += w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator+(const T2& w, MultiComplex<0, T1> z) { z += w; return z; }
// all substractions
template<class T1, class T2>
MultiComplex<0, T1> operator-(MultiComplex<0, T1> z, const MultiComplex<0, T2>& w) { z -= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator-(MultiComplex<0, T1> z, const T2& w) { z -= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator-(const T2& w, MultiComplex<0, T1> z) { z -= w; return z; }
// all multiplications
template<class T1, class T2>
MultiComplex<0, T1> operator*(MultiComplex<0, T1> z, const MultiComplex<0, T2>& w) { z *= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator*(MultiComplex<0, T1> z, const T2& w) { z *= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator*(const T2& w, MultiComplex<0, T1> z) { z *= w; return z; }
// all divisions
template<class T1, class T2>
MultiComplex<0, T1> operator/(MultiComplex<0, T1> z, const MultiComplex<0, T2>& w) { z /= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator/(MultiComplex<0, T1> z, const T2& w) { z /= w; return z; }
template<class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<0, T1> operator/(const T2& w, MultiComplex<0, T1> z) { z /= w; return z; }

template<class T1, class T2>
bool operator==(const MultiComplex<0, T1>& z, const MultiComplex<0, T2>& w) { return z.z1 == w.z1 && z.z2 == w.z2; }

template<class T1, class T2>
bool operator!=(const MultiComplex<0, T1>& z, const MultiComplex<0, T2>& w) { return !(z==w); }


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
    // create random MultiComplex with deterministic seed
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dist(-1, 1);
    z = MultiComplex<0, T>{dist(gen), dist(gen)};
}

template<unsigned N, class T>
void toOne(MultiComplex<N, T> &z) {
    toOne(z.real());
    toOne(z.imag());
}
template<class T>
void toOne(MultiComplex<0, T> &z) {
    z = MultiComplex<0, T>{1, 0};
}

template<unsigned N, class T, ENFORCE(N==0)>
MultiComplex<0, T> I() {
    return MultiComplex<0, T>{0, 1};
}
template<unsigned N, class T, ENFORCE(N==0)>
MultiComplex<0, T> R() {
    return MultiComplex<0, T>{1, 0};
}

template<unsigned N, class T, ENFORCE(N>0)>
MultiComplex<N, T> R() {
    return MultiComplex<N, T>{R<N-1, T>(), MultiComplex<N-1, T>()};
}
template<unsigned N, class T, ENFORCE(N>0)>
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
    MultiComplex<N, T> ret = val * R<N, T>();
    return ret;
}



};


