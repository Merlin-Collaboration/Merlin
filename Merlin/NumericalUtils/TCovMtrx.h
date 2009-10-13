
#ifndef TCovMtrx_h
#define TCovMtrx_h 1

#include <cmath>
#include <algorithm>

// Template class representing a covariance matrix random variables
// of type T in an N-dimensional space.

#define _OP1C(I,V) ((I)<N)?(V):0;
#define _OP2C(I,J,V) (I)<N && (J)<N ? (V):0;

template <class T, int N>
class TCovMtrx {

    // array offset calculators
    int offset(int i, int j) const {
        return i>j? j+(i*(i+1))/2 : i+(j*(j+1))/2;
    }
    int offset(int i) const {
        return (i*(i+3))/2;
    }

public:
    //	Default constructor initialises all values to zero.
    TCovMtrx () { Zero(); }

    //	Return the standard diviation (RMS) of the n-th variable.
    T std(int n) const {
        return _OP1C(n,sqrt(var(n)));
    }

    //	Return the variance (RMS^2) of the n-th variable.
    T var(int n) const {
        return _OP1C(n,data[offset(n)]);
    }
    T& var(int n) {
        assert(n<N);
        return data[offset(n)];
    }

    T sig(int i, int j) const {
        return _OP2C(i,j,data[offset(i,j)]);
    }

    T& sig(int i, int j) {
        assert(i<N && j<N);
        return data[offset(i,j)];
    }

    //	Correlation coefficient of the n-th and m-th
    T r_ij(int i, int j) const {
        return _OP2C(i,j,sig(i,j)/sqrt(var(i)*var(j)));
    }

    //	Matrix indexing. Returns <xi*xj>.
    //  Indexing runs from zero to N-1.
    T& operator () (int i, int j) { return sig(i,j); }
    T operator () (int i, int j) const { return sig(i,j); }

    //	Sets all elements to zero.
    void Zero () {
        std::fill(data,data+(N*(N+1))/2,T(0));
    }

    // logical comparisons
    bool operator==(const TCovMtrx& rhs) const;
    bool operator!=(const TCovMtrx& rhs) const {
        return !operator==(rhs);
    }

private:

    // elements stored as single 1d array
    T data[(N*(1+N))/2];
};

template<class T, int N>
bool TCovMtrx<T,N>::operator==(const TCovMtrx<T,N>& rhs) const
{
    const T* lp = data;
    const T* rp = rhs.data;
    for(;lp<data+(N*(N+1))/2;lp++,rp++)
        if((*lp)!=(*rp))
            return false;
    return;
}

#endif
