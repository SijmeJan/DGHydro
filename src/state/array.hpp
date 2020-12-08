#ifndef DG_ARRAY_HPP
#define DG_ARRAY_HPP

#include "../rhs/temp.hpp"

namespace DGHydro {

  template<class T, int N>
  class Array {
  public:
    // Default constructor
    Array() {
      data = nullptr;
    };

    // Copy constructor
    Array(const Array& s) {
      //std::cout << "Copy: allocating " << N << " elements of type " << type_name<decltype(data[0])>() << "\n";
      data = new T[N];
      std::copy(s.data, &s.data[N-1], data);
    }

    // Move constructor
    Array(Array&& s) {
      //std::cout << "Move constructor for type " << type_name<decltype(data[0])>() << "\n";
      data = s.data;
      s.data = nullptr;
    }

    // Constructor from double
    Array(const double& s) {
      //std::cout << "Copy from double: allocating " << N
      //          << " elements of type " << type_name<decltype(data[0])>() << "\n";
      data = new T[N];

      for (int i = 0; i < N; i++) data[i] = s;
    }

    ~Array() {
      if (data != nullptr) {
        //std::cout << "Deallocating " << N << " elements of type " << type_name<decltype(data[0])>() << "\n";

        delete[] data;
      }
    };

    // Copy assignment
    Array& operator=(const Array<T,N>& s) {
      //std::cout << "Copy assignment from " << type_name<decltype(s)>()
      //          << " to " << type_name<decltype(this)>() << "\n";

      if (data == nullptr) {
        //std::cout << "Copy assigment: allocating " << N
        //          << " elements of type " << type_name<decltype(data[0])>() << "\n";
        data = new T[N];
      }

      std::copy(s.data, &s.data[N-1], data);
      //for (int i = 0; i < N; i++) data[i] = s.data[i];

      return *this;
    }
    // Move assignment
    Array& operator=(Array&& s) {
      //std::cout << "Move assignment\n";

      if (this != &s) {
        if (data != nullptr) {
          //std::cout << "Deallocating " << N << " elements of type " << type_name<decltype(data[0])>() << " in move assignment\n";
          delete[] data;
        }

        data = s.data;
        s.data = nullptr;
      }
      return *this;
    }
    // Assign to constant
    Array& operator=(const double& s) {
      //std::cout << "Double assignment\n";

      if (data == nullptr) {
        //std::cout << "Allocating " << N << " elements of type " << type_name<decltype(data[0])>() << " in double assignment\n";
        data = new T[N];
      }
      for (int i = 0; i < N; i++) data[i] = s;
      return *this;
    }

    // Element-wise addition
    Array& operator+=(const Array& rhs) {
      std::cout << "Adding array\n";
      for (int i = 0; i < N; i++) data[i] += rhs.data[i];
      return *this;
    }
    friend Array operator+(Array lhs, const Array& rhs) {
      std::cout << "Adding array 2\n";
      lhs += rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Add scalar to every element
    Array& operator+=(const double& rhs) {
      std::cout << "Adding scalar\n";
      for (int i = 0; i < N; i++) data[i] += rhs;
      return *this;
    }
    friend Array operator+(Array lhs, const double& rhs) {
      std::cout << "Adding scalar 2\n";
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
      std::cout << "Multiplying array by double\n";
      for (int i = 0; i < N; i++) data[i] *= rhs;
      return *this;
    }
    friend Array operator*(Array lhs, const double& rhs) {
      std::cout << "Multiplying array by double 2\n";
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
    T *data = nullptr;
    //std::unique_ptr<T> data;
  };

  // Scalar addition, subtraction and multiplication from the left
  template<class T, int N>
  Array<T, N> operator+(double const& scalar, Array<T, N> rhs) {
    std::cout << "Adding scalar from left\n";
    return rhs += scalar; // calls rhs.operator+=(scalar);
  }
  template<class T, int N>
  Array<T, N> operator-(double const& scalar, Array<T, N> rhs) {
    return (-rhs) += scalar;
  }
  template<class T, int N>
  Array<T, N> operator*(double const& scalar, Array<T, N> rhs) {
    std::cout << "Multiplying array by double from left " << rhs[0] << "\n";
    return rhs *= scalar; // calls rhs.operator*=(scalar);
  }
}

#endif // DG_ARRAY_HPP
