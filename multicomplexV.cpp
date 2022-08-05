
#include <iostream>
#include <initializer_list>
#include <vector>
#include <complex>
#include <chrono>


template<class T>
struct MultiComplex  {
    size_t N;
    std::vector<std::complex<T>> values;

    MultiComplex() {}
    template<size_t provided>
    MultiComplex(size_t N_, const std::complex<T>(&values_)[provided]) : N(N_) {
        values.resize(1<<N_);
        std::copy(values_, values_ + provided, values.data());
    }
    MultiComplex(size_t N_, const std::complex<T>* values_) : N(N_) {
        values.resize(1<<N_);
        std::copy(values_, values_ + (1<<N_), values.data());
    }

    std::complex<T>& operator[](size_t i) { return values[i]; }
    const std::complex<T> operator[](size_t i) const { return values[i]; }

    MultiComplex operator+= (const MultiComplex& other) {    
        if (N < other.N) throw "Invalid operation";
        for (size_t i = 0; i < values.size(); i++) {
            values[i] += other.values[i];
        }
        return *this;
    }
    MultiComplex operator-= (const MultiComplex& other) {
        if (N < other.N) throw "Invalid operation";
        for (size_t i = 0; i < values.size(); i++) {
            values[i] -= other.values[i];
        }
        return *this;
    }
    MultiComplex operator*=(const MultiComplex& other) {
        if (N < other.N) {
            values.resize(other.size()); N = 1<<size();
            MultiComplex<T> tmp = real();
            real(tmp * other.real() - imag() * other.imag());       
            imag(tmp * other.imag() + imag() * other.real());
            
            }
        else if (other.N == N) {
            if (N == 0) values[0] *= other.values[0];
            else {
                
                MultiComplex<T> tmp = real();
                real(tmp * other.real() - imag() * other.imag());       
                imag(tmp * other.imag() + imag() * other.real());
 }
        } else {
            MultiComplex<T> proxy;
            for (size_t i=0; i<size(); i+=other.size()) {
                proxy = MultiComplex<T>(other.N, values.data() + i);
                proxy *= other; std::copy(proxy.values.begin(), proxy.values.end(), values.begin() + i);
            }   
        }
        return *this;
    }
    
    const size_t size() const { return values.size();}
    const size_t order()    const { return N;}

    MultiComplex imag() const {
        if (N == 0) throw "Invalid operation";
        return MultiComplex<T>(N-1, values.data() + size()/2);
    }
    void imag(const MultiComplex& other) {
        size_t n_copy = std::min(other.size(), size()/2);
        std::copy(other.values.data(), other.values.data() + n_copy, values.data() + size()/2);
    }
    MultiComplex real() const {
        return MultiComplex<T>(N-1, values.data());
    }
    void real(const MultiComplex& other) {
        size_t n_copy = std::min(other.size(), size()/2);
        std::copy(other.values.data(), other.values.data() + n_copy, values.data());
    }


    friend std::ostream& operator<<(std::ostream& os, const MultiComplex<T>& m) {
        for (size_t i = 0; i < m.size(); i++) {
            os << m[i] << " ";
        }
        return os;
    }

};

template<class T>
MultiComplex<T> operator+(const MultiComplex<T>& a, const MultiComplex<T>& b) {
    if (a.order() > b.order()) {
        MultiComplex<T> result(a);
        result += b;
        return result;
    } else {
        MultiComplex<T> result(b);
        result += a;
        return result;
    }
}
template<class T>
MultiComplex<T> operator-(const MultiComplex<T>& a, const MultiComplex<T>& b) {
    if (a.order() > b.order()) {
        MultiComplex<T> result(a);
        result -= b;
        return result;
    } else {
        MultiComplex<T> result(b);
        result *= MultiComplex<T>(0, {-1});
        result += a;
        return result;
    }
}
template<class T>
MultiComplex<T> operator*(const MultiComplex<T>& a, const MultiComplex<T>& b) {
    if (a.order() > b.order()) {
        MultiComplex<T> result(a);
        result *= b;
        return result;
    } else {
        MultiComplex<T> result(b);
        result *= a;
        return result;
    }
}


int main() {
    MultiComplex<double> b(2, {1,2,3,4});
    MultiComplex<double> a(1, {0, 1});
    // time and print multiplication a = a * a
    auto start = std::chrono::high_resolution_clock::now();
    a = b * a;
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms" << std::endl;

    std::cout << a << std::endl;
    return 0;
}




/*
template<class std::complex<T>>
struct MultiComplex {
    
    bool is_parent;
    size_t N; // true array size = 2 ^ N real numbers
    std::complex<T>* values;
   ~MultiComplex() {if (is_parent) delete[] values;}
    MultiComplex() : is_parent(false) {}
    MultiComplex(size_t N_) :  is_parent(true), N(N_), values(new std::complex<T>[N]) {}
    MultiComplex(size_t N_, std::complex<T>* values_, bool is_parent_ = false) : is_parent(is_parent_), N(N_) {
        if (is_parent_) {
            values = new std::complex<T>[N];
            std::copy(values_, values_ + (1<<N), values);
        } else {
            values = values_;
        }
    } // dangerous!
    MultiComplex(size_t N_, std::initializer_list<std::complex<T>> values_) : is_parent(true), N(N_) {
        values = new std::complex<T>[1 << N];
        std::copy(values_.begin(), values_.end(), values);
    }
    MultiComplex(const MultiComplex& other) {
        is_parent = true;
        N = other.N;
        values = new std::complex<T>[1 << N];
        std::copy(other.values, other.values + (1 << N), values);
        std::cout << "copying()" << std::endl;
    }
    MultiComplex(MultiComplex&& other) {
        if (other.is_parent) {
            is_parent = true;
            N = other.N;
            values = new std::complex<T>[1 << N];
            std::copy(other.values, other.values + (1 << N), values);
            std::cout << "from moving copying()" << std::endl;
        } else {
            is_parent = false;
            N = other.N;
            values = other.values;
            std::cout << "moving()" << std::endl;
        }
    }
    MultiComplex& operator=(const MultiComplex& other) {
        if (this == &other) return *this;
        if (is_parent) delete[] values;
        is_parent = true;
        N = other.N;
        values = new std::complex<T>[1 << N];
        std::copy(other.values, other.values + (1 << N), values);
        std::cout << "copying=" << std::endl;
        return *this;
    }
    MultiComplex& operator=(MultiComplex&& other) {
        if (is_parent) delete[] values;
        if (other.is_parent) {
            is_parent = true;
            N = other.N;
            values = new std::complex<T>[1 << N];
            std::copy(other.values, other.values + (1 << N), values);
            std::cout << "copying from moving=" << std::endl;
        } else {
            is_parent = false;
            N = other.N;
            values = other.values;
            std::cout << "moving=" << std::endl;
            
        }
        return *this;
    }
    MultiComplex& operator=(std::initializer_list<std::complex<T>> values_) {
        if (is_parent) delete[] values;
        is_parent = true;
        N = values_.size();
        values = new std::complex<T>[1 << N];
        std::copy(values_.begin(), values_.end(), values);
        return *this;
    }
    MultiComplex IMAG() {
        if (N == 0) return *this;
        MultiComplex child(N/2, values + (1<<N)); return child;
    }
    MultiComplex REAL() {
        MultiComplex child(N/2, values         ); return child;
    }
    MultiComplex imag() {
        if (N == 0) return *this;
        MultiComplex child(N/2, values + (1<<N), true); return child;
    }
    MultiComplex real() {
        MultiComplex child(N/2, values         , true); return child;
    }


};
*/
