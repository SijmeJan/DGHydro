#ifndef DG_ARRAY_HPP
#define DG_ARRAY_HPP

#include "../rhs/temp.hpp"

namespace DGHydro {

  template<class T, int N>
  class Array {
  public:
    // Default constructor
    Array() {
      data = new T[N];
    };
    // Copy constructor
    Array(const Array& s) {
      data = new T[N];
      for (int i = 0; i < N; i++) data[i] = s.data[i];
    }
    // Move constructor
    Array(Array&& s) {
      data = s.data;
      s.data = nullptr;
    }
    ~Array() {
      delete[] data;
    };

    // Copy assignment
    Array& operator=(const Array& s) {
      //data = new T[N];
      for (int i = 0; i < N; i++) data[i] = s.data[i];
      return *this;
    }
    // Move assignment
    Array& operator=(Array&& s) {
      if (this != &s) {
        data = s.data;
        s.data = nullptr;
      }
      return *this;
    }
    // Assign to constant
    Array& operator=(const double& s) {
      //data = new T[N];
      for (int i = 0; i < N; i++) data[i] = s;
      return *this;
    }

    // Element-wise addition
    Array& operator+=(const Array& rhs) {
      for (int i = 0; i < N; i++) data[i] += rhs.data[i];
      return *this;
    }
    friend Array operator+(Array lhs, const Array& rhs) {
      lhs += rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Add scalar to every element
    Array& operator+=(const double& rhs) {
      for (int i = 0; i < N; i++) data[i] += rhs;
      return *this;
    }
    friend Array operator+(Array lhs, const double& rhs) {
      lhs += rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise subtraction
    Array& operator-=(const Array& rhs) {
      for (int i = 0; i < N; i++) data[i] -= rhs.data[i];
      return *this;
    }
    friend Array operator-(Array lhs, const Array& rhs) {
      lhs -= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Subtract scalar from every element
    Array& operator-=(const double& rhs) {
      for (int i = 0; i < N; i++) data[i] -= rhs;
      return *this;
    }
    friend Array operator-(Array lhs, const double& rhs) {
      lhs -= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise multiplication
    Array& operator*=(const Array& rhs) {
      for (int i = 0; i < N; i++) data[i] *= rhs.data[i];
      return *this;
    }
    friend Array operator*(Array lhs, const Array& rhs) {
      lhs *= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Multiply every element by double
    Array& operator*=(const double& rhs) {
      for (int i = 0; i < N; i++) data[i] *= rhs;
      return *this;
    }
    friend Array operator*(Array lhs, const double& rhs) {
      lhs *= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise division
    Array& operator/=(const Array& rhs) {
      for (int i = 0; i < N; i++) data[i] /= rhs.data[i];
      return *this;
    }
    friend Array operator/(Array lhs, const Array& rhs) {
      lhs /= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Divide all elements by double
    Array& operator/=(const double& rhs) {
      for (int i = 0; i < N; i++) data[i] /= rhs;
      return *this;
    }
    friend Array operator/(Array lhs, const double& rhs) {
      lhs /= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Array subscript operator
    T& operator[](std::size_t idx) { return data[idx]; }
    const T& operator[](std::size_t idx) const { return data[idx]; }

    // Negation operator
    Array operator-() {
      Array s;
      for (int i = 0; i < N; i++) s.data[i] = -data[i];

      return s;
    }

  private:
    T *data;
  };

  // Scalar addition, subtraction and multiplication from the left
  template<class T, int N>
  Array<T, N> operator+(double const& scalar, Array<T, N> rhs) {
    return rhs += scalar; // calls rhs.operator+=(scalar);
  }
  template<class T, int N>
  Array<T, N> operator-(double const& scalar, Array<T, N> rhs) {
    return (-rhs) += scalar;
  }
  template<class T, int N>
  Array<T, N> operator*(double const& scalar, Array<T, N> rhs) {
    return rhs *= scalar; // calls rhs.operator*=(scalar);
  }
}

#endif // DG_ARRAY_HPP
