#include <iostream>
#include <complex.h>

using Complex = std::complex<double>;

template<class T> class MultiComplex;

template<class V>
class MultiComplex<MultiComplex<V>>{
public:
    MultiComplex<V> z1, z2;
    MultiComplex(MultiComplex<V> z1_, MultiComplex<V> z2_) : z1(z1_), z2(z2_) {}  

    // overload operator +
        // MultiComplex += MultiComplex
    MultiComplex<MultiComplex<V>>& operator+=(const MultiComplex<MultiComplex<V>>& w) {
        z1 += w.z1;
        z2 += w.z2;
        return *this;
    }
        // MultiComplex += Complex
    MultiComplex<MultiComplex<V>>& operator+=(const Complex& w) {
        z1 += w;
        return *this;
    }
    // overload operator -
        // MultiComplex += MultiComplex
    MultiComplex<MultiComplex<V>>& operator-=(const MultiComplex<MultiComplex<V>>& w) {
        z1 -= w.z1;
        z2 -= w.z2;
        return *this;
    }
        // MultiComplex += Complex
    MultiComplex<MultiComplex<V>>& operator-=(const Complex& w) {
        z1 -= w;
        return *this;
    }
    // overload operator *
        // MultiComplex *= MultiComplex
    MultiComplex<MultiComplex<V>>& operator*=(const MultiComplex<MultiComplex<V>>& w) {
        MultiComplex<V> tmp = z1 * w.z1 - z2 * w.z2;
        z2  = z1 * w.z2 + z2 * w.z1;
        z1 = tmp;
        return *this;
    }
        // MultiComplex * Complex
    MultiComplex<MultiComplex<V>>& operator*=(const Complex& w) {
        z1 *= w;
        z2 *= w;
        return *this;
    }
};
template<class V>
MultiComplex<MultiComplex<V>> operator+(MultiComplex<MultiComplex<V>> z, const MultiComplex<MultiComplex<V>>& w) {
    z += w;
    return z;
}
template<class V>
MultiComplex<MultiComplex<V>> operator+(MultiComplex<MultiComplex<V>> z, const Complex& w) {
    z += w;
    return z;
}
template<class V>
MultiComplex<MultiComplex<V>> operator-(MultiComplex<MultiComplex<V>> z, const MultiComplex<MultiComplex<V>>& w) {
    z -= w;
    return z;
}
template<class V>
MultiComplex<MultiComplex<V>> operator-(MultiComplex<MultiComplex<V>> z, const Complex& w) {
    z -= w;
    return z;
}
template<class V>
MultiComplex<MultiComplex<V>> operator*(MultiComplex<MultiComplex<V>> z, const MultiComplex<MultiComplex<V>>& w) {
    z *= w;
    return z;
}
template<class V>
MultiComplex<MultiComplex<V>> operator*(MultiComplex<MultiComplex<V>> z, const Complex& w) {
    z *= w;
    return z;
}
template<>
class MultiComplex<Complex>{
public:
    Complex z1, z2;
    MultiComplex(Complex z1_, Complex z2_) : z1(z1_), z2(z2_) {}  
    // overload operator +
        // MultiComplex += MultiComplex
    MultiComplex<Complex>& operator+=(const MultiComplex<Complex>& w) {
        z1 += w.z1;
        z2 += w.z2;
        return *this;
    }
        // MultiComplex += Complex
    MultiComplex<Complex>& operator+=(const Complex& w) {
        z1 += w;
        return *this;
    }
    // overload operator -
        // MultiComplex += MultiComplex
    MultiComplex<Complex>& operator-=(const MultiComplex<Complex>& w) {
        z1 -= w.z1;
        z2 -= w.z2;
        return *this;
    }
        // MultiComplex += Complex
    MultiComplex<Complex>& operator-=(const Complex& w) {
        z1 -= w;
        return *this;
    }
    // overload operator *
        // MultiComplex *= MultiComplex
    MultiComplex<Complex>& operator*=(const MultiComplex<Complex>& w) {
        Complex tmp = z1 * w.z1 - z2 * w.z2;
        z2  = z1 * w.z2 + z2 * w.z1;
        z1 = tmp;
        return *this;
    }
        // MultiComplex * Complex
    MultiComplex<Complex>& operator*=(const Complex& w) {
        z1 *= w;
        z2 *= w;
        return *this;
    }
};
template<class V>
MultiComplex<MultiComplex<V>> operator+(MultiComplex<MultiComplex<V>> z, const MultiComplex<MultiComplex<V>>& w) {
    z += w;
    return z;
}
template<class V>
MultiComplex<MultiComplex<V>> operator+(MultiComplex<MultiComplex<V>> z, const Complex& w) {
    z += w;
    return z;
}
template<class V>
MultiComplex<MultiComplex<V>> operator-(MultiComplex<MultiComplex<V>> z, const MultiComplex<MultiComplex<V>>& w) {
    z -= w;
    return z;
}
template<class V>
MultiComplex<MultiComplex<V>> operator-(MultiComplex<MultiComplex<V>> z, const Complex& w) {
    z -= w;
    return z;
}
template<class V>
MultiComplex<Complex> operator*(MultiComplex<Complex> z, const MultiComplex<Complex>& w) {
    z *= w;
    return z;
}
template<class V>
MultiComplex<Complex> operator*(MultiComplex<Complex> z, const Complex& w) {
    z *= w;
    return z;
}





template<class T> // rape the compiler and the c++ standard lmao
MultiComplex<T> multiComplex(T z1, T z2) {
    return MultiComplex<T>(z1, z2);
}

int main()
{
    auto z1 = multiComplex(Complex(1, 2), Complex(3, 4));
    auto z2 = multiComplex(Complex(5, 6), Complex(7, 8));
    auto z3 = multiComplex(z1, z2);
    auto z4 = z1 * z2;

    std::cout << z4.z1 << std::endl;

    return 0;
}
