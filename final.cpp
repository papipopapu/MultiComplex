#include<iostream>
#include<string>
#include<array>
#include<complex>
#include <algorithm>
#include <iterator>
#include <random>
#include <chrono>

template<size_t N, class T>
using inner = std::array<std::complex<T>, N>;




#define ENFORCE(x) typename = typename std::enable_if<(x)>::type

template<size_t N, class T>
struct MultiComplex;

////////////////////////////////////////////////////////////////////////////////
template<size_t N, class T1, class T2>
typename std::enable_if<(N==1), void>::type
inner_timesequals(std::complex<T1> *z, const std::complex<T2> *w) { *z *= *w; }

template<size_t N, class T1, class T2>
typename std::enable_if<(N>1), void>::type
inner_timesequals(std::complex<T1> *z, const std::complex<T2> *w) {

    std::complex<T1> cache_z1[N/2];
    std::complex<T1> cache_z2[N/2];

    std::copy(z, z + N/2, cache_z1);
    std::copy(z + N/2, z + N, cache_z2);

    inner_timesequals<N/2, T1, T2>(z, w);
    inner_timesequals<N/2, T1, T2>(cache_z2, w + N/2);
    for (size_t i = 0; i < N/2; i++) { z[i] -= cache_z2[i]; }

    inner_timesequals<N/2, T1, T2>(cache_z1, w + N/2);
    inner_timesequals<N/2, T1, T2>(z + N/2, w);
    for (size_t i = 0; i < N/2; i++) { z[i + N/2] += cache_z1[i]; }
}
template<size_t N, class T1, class T2>
typename std::enable_if<(N==1), void>::type
inner_divequals(std::complex<T1> *z, const std::complex<T2> *w) { *z /= *w; }

template<size_t N, class T1, class T2>
typename std::enable_if<(N>1), void>::type
inner_divequals(std::complex<T1> *z, const std::complex<T2> *w) {

    std::complex<T1> cache_z1[N/2];
    std::complex<T1> cache_z2[N/2];

    std::copy(z, z + N/2, cache_z1);
    std::copy(z + N/2, z + N, cache_z2);

    inner_timesequals<N/2, T1, T2>(z, w);
    inner_timesequals<N/2, T1, T2>(cache_z2, w + N/2);
    for (size_t i = 0; i < N/2; i++) { z[i] += cache_z2[i]; }

    inner_timesequals<N/2, T1, T2>(cache_z1, w + N/2);
    inner_timesequals<N/2, T1, T2>(z + N/2, w);
    for (size_t i = 0; i < N/2; i++) { z[i + N/2] -= cache_z1[i]; }

    // now we have numerator, we need to build denominator
    std::copy(w, w + N/2, cache_z1);
    std::copy(w + N/2, w + N, cache_z2);
    inner_timesequals<N/2, T1, T2>(cache_z1, w);
    inner_timesequals<N/2, T1, T2>(cache_z2, w + N/2);
    for (size_t i = 0; i < N/2; i++) { cache_z1[i] += cache_z2[i]; }

    // divide!
    inner_divequals<N/2, T1, T2>(z, cache_z1);
    inner_divequals<N/2, T1, T2>(z + N/2, cache_z1);
}

template<size_t N, class T>
struct MultiComplex {
public:
    inner<(1<<N), T> values;
    
    MultiComplex<N, T> operator-() const { MultiComplex<N, T> ret;
        for (size_t i = 0; i < (1<<N); i++) { ret.values[i] = -values[i]; } return  ret; }

    MultiComplex<N, T> operator+() const { MultiComplex<N, T> ret = *this; return  ret; }


    template<class T2>
    MultiComplex<N, T>& operator+=(const MultiComplex<N, T2>& other) {
        for(size_t i = 0; i < (1<<N); i++) {
            values[i] += other.values[i]; } return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator+=(const T2& other) {
        values[0] += other; return *this; }
    template<class T2>
    MultiComplex<N, T>& operator+=(const std::complex<T2>& other) {
        values[0] += other; return *this; }

    template<class T2>
    MultiComplex<N, T>& operator-=(const MultiComplex<N, T2>& other) {
        for(size_t i = 0; i < (1<<N); i--) {
            values[i] -= other.values[i]; } return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator-=(const T2& other) {
        values[0] -= other; return *this; }
    template<class T2>
    MultiComplex<N, T>& operator-=(const std::complex<T2>& other) {
        values[0] -= other; return *this; }

    template<class T2>
    MultiComplex<N, T>& operator*=(const MultiComplex<N, T2>& other) {
        inner_timesequals<(1<<N), T, T2>(values.data(), other.values.data()); return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator*=(const T2& other) {
        values[0] *= other; return *this; }
    template<class T2>
    MultiComplex<N, T>& operator*=(const std::complex<T2>& other) {
        values[0] *= other; return *this; }

    template<class T2>
    MultiComplex<N, T>& operator/=(const MultiComplex<N, T2>& other) {
        inner_divequals<(1<<N), T, T2>(values.data(), other.values.data()); return *this; }
    template<class T2, ENFORCE(std::is_arithmetic<T2>::value)>
    MultiComplex<N, T>& operator/=(const T2& other) {
        values[0] /= other; return *this; }
    template<class T2>
    MultiComplex<N, T>& operator/=(const std::complex<T2>& other) {
        values[0] /= other; return *this; }

    friend std::ostream& operator<<(std::ostream& os, const MultiComplex<N, T>& z) {
        os << "{ ";
        for(size_t i = 0; i < (1<<N); i++) {
            os << z.values[i] << " "; } os << "}"; return os; }
};
// all additions
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator+(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z += w; return z; }
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator+(MultiComplex<N, T1> z, const std::complex<T2>& w) { z += w; return z; }
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator+(const std::complex<T2>& w, MultiComplex<N, T1> z) { z += w; return z; }
template<size_t N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator+(MultiComplex<N, T1> z, const T2& w) { z += w; return z; }
template<size_t N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator+(const T2& w, MultiComplex<N, T1> z) { z += w; return z; }
// all substractions
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator-(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z -= w; return z; }
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator-(MultiComplex<N, T1> z, const std::complex<T2>& w) { z -= w; return z; }
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator-(const std::complex<T2>& w, MultiComplex<N, T1> z) { z -= w; return z; }
template<size_t N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator-(MultiComplex<N, T1> z, const T2& w) { z -= w; return z; }
template<size_t N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator-(const T2& w, MultiComplex<N, T1> z) { z -= w; return z; }
// all multiplications
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator*(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z *= w; return z; }
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator*(MultiComplex<N, T1> z, const std::complex<T2>& w) { z *= w; return z; }
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator*(const std::complex<T2>& w, MultiComplex<N, T1> z) { z *= w; return z; }
template<size_t N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator*(MultiComplex<N, T1> z, const T2& w) { z *= w; return z; }
template<size_t N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator*(const T2& w, MultiComplex<N, T1> z) { z *= w; return z; }
// all divisions
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator/(MultiComplex<N, T1> z, const MultiComplex<N, T2>& w) { z /= w; return z; }
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator/(MultiComplex<N, T1> z, const std::complex<T2>& w) { z /= w; return z; }
template<size_t N, class T1, class T2>
MultiComplex<N, T1> operator/(const std::complex<T2>& w, MultiComplex<N, T1> z) { z /= w; return z; }
template<size_t N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator/(MultiComplex<N, T1> z, const T2& w) { z /= w; return z; }
template<size_t N, class T1, class T2, ENFORCE(std::is_arithmetic<T2>::value)>
MultiComplex<N, T1> operator/(const T2& w, MultiComplex<N, T1> z) { z /= w; return z; }

// randomize MultiComplex<N, T>
template<size_t N, class T>
void randomize(MultiComplex<N, T>& z) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<T> dis(-1, 1);
    for(size_t i = 0; i < (1<<N); i++) {
        z.values[i] = std::complex<T>(dis(gen), dis(gen)); }
}


int main() {
    using cx = std::complex<double>;

    MultiComplex<10, long double> z3 = {cx(6, 3)}; // from 10 up needs long double precision
    MultiComplex<10, long double> z4 = {cx(7, -5)}; 

    z3 /= z4;

    std::cout << z3 << std::endl;

  

  
    
    return 0;
}







