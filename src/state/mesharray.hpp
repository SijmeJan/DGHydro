#ifndef DG_MESHARRAY_HPP
#define DG_MESHARRAY_HPP

#include "./dynarray.hpp"

namespace DGHydro {

  template<class T>
  class MeshArray: public DynArray<T> {
  public:
    // Default constructor
    MeshArray(int Nx, int Ny, int Nz) : DynArray<T>(Nx*Ny*Nz), Nx(Nx), Ny(Ny), Nz(Nz)  {
      //std::cout << "MeshArray default constructor\n";
      N = Nx*Ny*Nz;
      //data = new T[N];
    };
    // Copy constructor
    MeshArray(const MeshArray& s) : DynArray<T>(s) {
      //std::cout << "MeshArray copy constructor\n";
      Nx = s.Nx;
      Ny = s.Ny;
      Nz = s.Nz;
      N = Nx*Ny*Nz;

      data = new T[N];
      for (int i = 0; i < N; i++) data[i] = s.data[i];
    }

    // Move constructor
    MeshArray(MeshArray&& s) : DynArray<T>(s) {
      //std::cout << "MeshArray move constructor\n";
      Nx = s.Nx;
      Ny = s.Ny;
      Nz = s.Nz;
      N = Nx*Ny*Nz;

      data = s.data;
      s.data = nullptr;
    }

    ~MeshArray() {
      //std::cout << "MeshArray destructor" << std::endl;
      //delete[] data;
    };


    // Copy assignment
    MeshArray& operator=(const MeshArray& s) {
      //std::cout << "MeshArray copy assignment\n";
      Nx = s.Nx;
      Ny = s.Ny;
      Nz = s.Nz;
      N = Nx*Ny*Nz;

      data = new T[N];
      for (int i = 0; i < N; i++) data[i] = s.data[i];

      return *this;
    }
    // Move assignment
    MeshArray& operator=(MeshArray&& s) {
      //std::cout << "MeshArray move assigment\n";
      if (this != &s) {
        Nx = s.Nx;
        Ny = s.Ny;
        Nz = s.Nz;
        N = Nx*Ny*Nz;

        data = s.data;
        s.data = nullptr;
      }
      return *this;
    }
    // Assign to constant
    MeshArray& operator=(const double& s) {
      //std::cout << "MeshArray constant assignment\n";

      data = new T[N];
      for (int i = 0; i < N; i++) data[i] = s;

      return *this;
    }

    /*
    // Element-wise addition
    MeshArray& operator+=(const MeshArray& rhs) {
      data += rhs.data;
      return *this;
    }
    friend MeshArray operator+(MeshArray lhs, const MeshArray& rhs) {
      lhs += rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Add scalar to every element
    MeshArray& operator+=(const double& rhs) {
      data += rhs;
      return *this;
    }
    friend MeshArray operator+(MeshArray lhs, const double& rhs) {
      lhs += rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise subtraction
    MeshArray& operator-=(const MeshArray& rhs) {
      data -= rhs.data;
      return *this;
    }
    friend MeshArray operator-(MeshArray lhs, const MeshArray& rhs) {
      lhs -= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Subtract scalar from every element
    MeshArray& operator-=(const double& rhs) {
      data -= rhs;
      return *this;
    }
    friend MeshArray operator-(MeshArray lhs, const double& rhs) {
      lhs -= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise multiplication
    MeshArray& operator*=(const MeshArray& rhs) {
      data[i] *= rhs.data;
      return *this;
    }
    friend MeshArray operator*(MeshArray lhs, const MeshArray& rhs) {
      lhs *= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Multiply every element by double
    MeshArray& operator*=(const double& rhs) {
      data *= rhs;
      return *this;
    }
    friend MeshArray operator*(MeshArray lhs, const double& rhs) {
      lhs *= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise division
    MeshArray& operator/=(const MeshArray& rhs) {
      data /= rhs.data;
      return *this;
    }
    friend MeshArray operator/(MeshArray lhs, const MeshArray& rhs) {
      lhs /= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Divide all elements by double
    MeshArray& operator/=(const double& rhs) {
      data /= rhs;
      return *this;
    }
    friend MeshArray operator/(MeshArray lhs, const double& rhs) {
      lhs /= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // MeshArray subscript operator
    T& operator[](std::size_t idx) { return data[idx]; }
    const T& operator[](std::size_t idx) const { return data[idx]; }

    // Negation operator
    MeshArray operator-() {
      MeshArray s;
      s.data = -data;

      return s;
    }
    */
  protected:
    using DynArray<T>::N;
    using DynArray<T>::data;

    int Nx, Ny, Nz;
    //T *data;
  };

  /*
  // Scalar addition, subtraction and multiplication from the left
  template<class T>
  MeshArray<T> operator+(double const& scalar, MeshArray<T> rhs) {
    return rhs += scalar; // calls rhs.operator+=(scalar);
  }
  template<class T>
  MeshArray<T> operator-(double const& scalar, MeshArray<T> rhs) {
    return (-rhs) += scalar;
  }
  template<class T>
  MeshArray<T> operator*(double const& scalar, MeshArray<T> rhs) {
    return rhs *= scalar; // calls rhs.operator*=(scalar);
  }
  */
}

#endif // DG_MESHARRAY_HPP
