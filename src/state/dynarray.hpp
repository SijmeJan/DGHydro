#ifndef DG_DYNARRAY_HPP
#define DG_DYNARRAY_HPP


namespace DGHydro {

  template<class T>
  class DynArray {
  public:
    // Default constructor
    DynArray(int N) : N(N) {
      //std::cout << "Default constructor" << std::endl;
      data = nullptr;
    };
    // Copy constructor
    DynArray(const DynArray& s) {
      N = s.N;

      data = new T[N];
      for (int i = 0; i < N; i++) data[i] = s.data[i];

      //std::cout << "Copy constructor for " << s.N << " elements of type "
      //          << type_name<decltype(data[0])>() << "\n";

    }
    // Move constructor
    DynArray(DynArray&& s) {
      //std::cout << "Move constructor" << std::endl;
      N = s.N;

      data = s.data;
      s.data = nullptr;
    }
    ~DynArray() {
      //std::cout << "Destructor" << std::endl;
      delete[] data;
    };

    // Copy assignment
    DynArray& operator=(const DynArray& s) {
      //std::cout << "Copy assignment" << std::endl;
      N = s.N;

      if (data == nullptr)
        data = new T[N];
      for (int i = 0; i < N; i++) data[i] = s.data[i];
      return *this;
    }
    // Move assignment
    DynArray& operator=(DynArray&& s) {
      //std::cout << "Move assignment" << std::endl;
      if (this != &s) {
        if (data != nullptr)
          delete[] data;

        N = s.N;
        data = s.data;
        s.data = nullptr;
      }
      return *this;
    }
    // Assign to constant
    DynArray& operator=(const double& s) {
      //std::cout << "Constant assignment" << std::endl;
      if (data == nullptr)
        data = new T[N];

      for (int i = 0; i < N; i++) data[i] = s;
      return *this;
    }

    // Element-wise addition
    DynArray& operator+=(const DynArray& rhs) {
      //std::cout << "Addition 1" << std::endl;
      for (int i = 0; i < N; i++) data[i] += rhs.data[i];
      return *this;
    }
    friend DynArray operator+(DynArray lhs, const DynArray& rhs) {
      //std::cout << "Addition 2" << std::endl;
      lhs += rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Add scalar to every element
    DynArray& operator+=(const double& rhs) {
      //std::cout << "Addition scalar 1" << std::endl;
      for (int i = 0; i < N; i++) data[i] += rhs;
      return *this;
    }
    friend DynArray operator+(DynArray lhs, const double& rhs) {
      //std::cout << "Addition scalar 2" << std::endl;
      lhs += rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise subtraction
    DynArray& operator-=(const DynArray& rhs) {
      for (int i = 0; i < N; i++) data[i] -= rhs.data[i];
      return *this;
    }
    friend DynArray operator-(DynArray lhs, const DynArray& rhs) {
      lhs -= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Subtract scalar from every element
    DynArray& operator-=(const double& rhs) {
      for (int i = 0; i < N; i++) data[i] -= rhs;
      return *this;
    }
    friend DynArray operator-(DynArray lhs, const double& rhs) {
      lhs -= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise multiplication
    DynArray& operator*=(const DynArray& rhs) {
      //std::cout << "Multiplication 1" << std::endl;
      for (int i = 0; i < N; i++) data[i] *= rhs.data[i];
      return *this;
    }
    friend DynArray operator*(DynArray lhs, const DynArray& rhs) {
      //std::cout << "Multiplication 2" << std::endl;
      lhs *= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Multiply every element by double
    DynArray& operator*=(const double& rhs) {
      //std::cout << "Multiplication scalar 1 " << N << " " << data << std::endl;
      for (int i = 0; i < N; i++) data[i] *= rhs;
      return *this;
    }
    friend DynArray operator*(DynArray lhs, const double& rhs) {
      lhs *= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // Element-wise division
    DynArray& operator/=(const DynArray& rhs) {
      for (int i = 0; i < N; i++) data[i] /= rhs.data[i];
      return *this;
    }
    friend DynArray operator/(DynArray lhs, const DynArray& rhs) {
      lhs /= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }
    // Divide all elements by double
    DynArray& operator/=(const double& rhs) {
      for (int i = 0; i < N; i++) data[i] /= rhs;
      return *this;
    }
    friend DynArray operator/(DynArray lhs, const double& rhs) {
      lhs /= rhs; // reuse compound assignment
      return lhs; // return the result by value (uses move constructor)
    }

    // DynArray subscript operator
    T& operator[](std::size_t idx) { return data[idx]; }
    const T& operator[](std::size_t idx) const { return data[idx]; }

    // Negation operator
    DynArray operator-() {
      DynArray s;
      for (int i = 0; i < N; i++) s.data[i] = -data[i];

      return s;
    }

  protected:
    int N;
    T *data = nullptr;
  };

  // Scalar addition, subtraction and multiplication from the left
  template<class T>
  DynArray<T>& operator+(double const& scalar, DynArray<T> rhs) {
    return rhs += scalar; // calls rhs.operator+=(scalar);
  }
  template<class T>
  DynArray<T>& operator-(double const& scalar, DynArray<T> rhs) {
    return (-rhs) += scalar;
  }
  template<class T>
  DynArray<T>& operator*(double const& scalar, DynArray<T> rhs) {
    //std::cout << "Multiplication scalar left" << std::endl;
     return rhs *= scalar; // calls rhs.operator*=(scalar);
  }
}

#endif // DG_DYNARRAY_HPP
