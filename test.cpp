
#include <iostream>
#include <initializer_list>




template<class T>
struct MultiComplex {
    
    bool is_parent;
    size_t N; // true array size = 2 ^ N real numbers
    T* values;
   ~MultiComplex() {if (is_parent) delete[] values;}
    MultiComplex() : is_parent(false) {}
    MultiComplex(size_t N_) :  is_parent(true), N(N_), values(new T[N]) {}
    MultiComplex(size_t N_, T* values_, bool is_parent_ = false) : is_parent(is_parent_), N(N_), values(values_) {} // dangerous!
    MultiComplex(size_t N_, std::initializer_list<T> values_) : is_parent(true), N(N_) {
        values = new T[1 << N];
        std::copy(values_.begin(), values_.end(), values);
    }
    MultiComplex(const MultiComplex& other) {
        is_parent = true;
        N = other.N;
        values = new T[1 << N];
        std::copy(other.values, other.values + (1 << N), values);
    }
    MultiComplex(MultiComplex&& other) {
        if (other.is_parent) {
            is_parent = true;
            N = other.N;
            values = new T[1 << N];
            std::copy(other.values, other.values + (1 << N), values);
        } else {
            is_parent = false;
            N = other.N;
            values = other.values;
        }
    }
    MultiComplex& operator=(const MultiComplex& other) {
        if (is_parent) delete[] values;
        is_parent = true;
        N = other.N;
        values = new T[1 << N];
        std::copy(other.values, other.values + (1 << N), values);
        return *this;
    }
    MultiComplex& operator=(MultiComplex&& other) {
        if (is_parent) delete[] values;
        if (other.is_parent) {
            is_parent = true;
            N = other.N;
            values = new T[1 << N];
            std::copy(other.values, other.values + (1 << N), values);
        } else {
            is_parent = false;
            N = other.N;
            values = other.values;
            std::cout << "moving" << std::endl;
        }
        return *this;
    }
    MultiComplex& operator=(std::initializer_list<T> values_) {
        if (is_parent) delete[] values;
        is_parent = true;
        N = values_.size();
        values = new T[1 << N];
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


};



MultiComplex<double> f(const MultiComplex<double>& x) {
    return x;
}


int main() {
    
    MultiComplex<double> y(2, {1,2,3,4});
    MultiComplex<double> x;
    x = y.IMAG();
    x.values[0] = 69;
    std::cout << y.values[0] << std::endl;


}