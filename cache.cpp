#include <iostream>
#include <string>
#include <complex.h>

using Complex = std::complex<double>;

template<class V>
class MultiComplex {
private:
    V z1_, z2_;
public:
    template<class T>
    MultiComplex(MultiComplex<T> z1_, MultiComplex<T> z2_) : z1_(z1_), z2_(z2_) {}
    MultiComplex(Complex z1_, Complex z2_) : z1_(z1_), z2_(z2_) {}
    MultiComplex() {}

    V& z1() { return z1_; } 
    V& z2() { return z2_; }
    V z1() const { return z1_; }
    V z2() const { return z2_; }

    template<unsigned N>
    friend void swap(MultiComplex<N>& lhs, MultiComplex<N>& rhs) {
        using std::swap;
        swap(lhs.z1, rhs.z1);
        swap(lhs.z2, rhs.z2);
    }
    // overload operator unary -
    MultiComplex<N> operator-() const {
        return MultiComplex<N>(-z1_, -z2_);
    }
    // overload operator unary +
    MultiComplex<N> operator+() const {
        return MultiComplex<N>(z1_, z2_);
    }

    // overload operator = 
    MultiComplex<N>& operator=(MultiComplex<N> w) {
        swap(*this, w);
        return *this;
    }
    // overload operator +
        // MultiComplex += MultiComplex
    MultiComplex<V>& operator+=(const MultiComplex<V>& w) {
        z1_ += w.z1_;
        z2_ += w.z2_;
        return *this;
    }
        // MultiComplex += Complex
    MultiComplex<V>& operator+=(const Complex& w) {
        z1_ += w;
        return *this;
    }
    // overload operator -
        // MultiComplex += MultiComplex
    MultiComplex<V>& operator-=(const MultiComplex<V>& w) {
        z1_ -= w.z1_;
        z2_ -= w.z2_;
        return *this;
    }
        // MultiComplex += Complex
    MultiComplex<V>& operator-=(const Complex& w) {
        z1_ -= w;
        return *this;
    }
    // overload operator *
        // MultiComplex *= MultiComplex
    MultiComplex<V>& operator*=(const MultiComplex<V>& w) {
        V tmp = z1_ * w.z1_ - z2_ * w.z2_;
        z2_  = z1_ * w.z2_ + z2_ * w.z1_;
        z1_ = tmp;
        return *this;
    }
        // MultiComplex * Complex
    MultiComplex<V>& operator*=(const Complex& w) {
        z1_ *= w;
        z2_ *= w;
        return *this;
    }
};


template<class V>
MultiComplex<V> operator+(MultiComplex<V> z, const MultiComplex<V>& w) {
    z += w;
    return z;
}

template<class V> 
MultiComplex<MultiComplex<V>> operator+(const MultiComplex<V>& w, MultiComplex<MultiComplex<V>> z) {
    z.z1() += w;
    return z;
}
template<class V>
MultiComplex<V> operator+(MultiComplex<V> z, const Complex& w) {
    z += w;
    return z;
}
template<class V>
MultiComplex<V> operator+(const Complex& w, MultiComplex<V> z) {
    z += w;
    return z;
}
template<class V>
MultiComplex<V> operator-(MultiComplex<V> z, const MultiComplex<V>& w) {
    z -= w;
    return z;
}
template<class V>
MultiComplex<V> operator-(MultiComplex<V> z, const Complex& w) {
    z -= w;
    return z;
}
template<class V>
MultiComplex<V> operator-(const Complex& w, MultiComplex<V> z) {
    z = -z;
    z += w;
    return z;
}
template<class V>
MultiComplex<V> operator*(MultiComplex<V> z, const MultiComplex<V>& w) {
    z *= w;
    return z;
}
template<class V>
MultiComplex<V> operator*(MultiComplex<V> z, const Complex& w) {
    z *= w;
    return z;
}
template<class V>
MultiComplex<V> operator*(const Complex& w, MultiComplex<V> z) {
    z *= w;
    return z;
}
template<class V>
bool operator==(const MultiComplex<V>& z, const MultiComplex<V>& w) {
    return z.z1 == w.z1 && z.z2 == w.z2;
}
template<class V>
bool operator!=(const MultiComplex<V>& z, const MultiComplex<V>& w) {
    return !(z == w);
}
template<class V>
bool operator< (const MultiComplex<V>& z, const MultiComplex<V>& w) {
    return z.z1 < w.z1;
}
template<class V>
bool operator> (const MultiComplex<V>& z, const MultiComplex<V>& w) {
    return w.z1 < z.z1;
}
template<class V>
bool operator<=(const MultiComplex<V>& z, const MultiComplex<V>& w) {
    return !(z > w);
}
template<class V>
bool operator>=(const MultiComplex<V>& z, const MultiComplex<V>& w) {
    return !(z < w);
}

template<class T> // rape the compiler and the c++11 standard lmao
MultiComplex<T> multiComplex(T z1, T z2) {
    return MultiComplex<T>(z1, z2);
}
typedef typename
        conditional<N==0, Complex, MultiComplex<N-1>>::type
        Bruh;
template<class T> 
MultiComplex<T> getNextZero() {
    MultiComplex<T> compare;
    MultiComplex<Complex> ret(0, 0);
    while (typeid(ret) != typeid(compare)) {
        auto temp = ret;
        auto ret = multiComplex(temp, temp);
    }
    return ret;
}

template <bool, typename T, typename F>
struct conditional {
    typedef T type;
};
template <typename T, typename F>
struct conditional<false, T, F> {
    typedef F type;
};

int main()
{
    auto z1 = multiComplex(Complex(1, 2), Complex(3, 4));
    auto z2 = multiComplex(Complex(5, 6), Complex(7, 8));
    auto z3 = multiComplex(z1, z2);
    auto z4 = multiComplex(z2, z1);
    auto z5 = multiComplex(z3, z4);

    //int N = getN(z5);
    std::cout << z5.z1().z1().z1() << std::endl;
    //std::cout << "N: " << N << std::endl;



   
    return 0;
}
