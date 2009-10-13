


#ifndef Space3D_h
#define Space3D_h 1

#include "merlin_config.h"


// VectorTags
#include "EuclideanGeometry/VectorTags.h"



//	Template class for generating arithmetic three vectors
//	(Euclidean) of the form (x,y,z). The template parameters
//	are the type of data storage (double, float etc.) and an
//	additional tag value. The tag class acts as a type-safe
//	mechanism for distinguising different types of three
//	vector (see Point3D and Vector3D).


template <class T, class tag>
class TVec3D
{
public:
    //	Default constructor. Data is not initialised.
    TVec3D ();

    //	Copy constructor.
    TVec3D (const TVec3D<T,tag>& v);

    //	Explicit constructor from the three vector components.
    TVec3D (const T& x1, const T& y1, const T& z1);

    bool operator==(const TVec3D< T,tag > &right) const;

    bool operator!=(const TVec3D< T,tag > &right) const;


    //	Copy assignment
    const TVec3D<T,tag>& operator = (const TVec3D<T,tag>& v);

    //	Arithmetic assignment.
    const TVec3D<T,tag>& operator += (const TVec3D<T,tag>& v);

    const TVec3D<T,tag>& operator -= (const TVec3D<T,tag>& v);

    const TVec3D<T,tag>& operator *= (const T& s);

    const TVec3D<T,tag>& operator /= (const T& s);

    TVec3D operator - () const;

    TVec3D operator - (const TVec3D<T,tag>& v) const;

    //	Dot (inner) product.
    T dot (const TVec3D<T,tag>& v) const;

    //	Return true if all components are zero.
    bool isZero () const;

    TVec3D operator + (const TVec3D<T,tag>& v) const;

    TVec3D operator * (T s) const;

    TVec3D operator / (T s) const;

    // Data Members for Class Attributes

    T x;

    T y;

    T z;

protected:
private:
private:
};

//	A 3-d vector (x,y,z).


typedef TVec3D< double,VectorTag > Vector3D;

//	A (x,y,z) point.


typedef TVec3D< double,PointTag > Point3D;

// Parameterized Class TVec3D

template <class T, class tag>
inline TVec3D<T,tag>::TVec3D ()
{
}

template <class T, class tag>
inline TVec3D<T,tag>::TVec3D (const TVec3D<T,tag>& v)
        :x(v.x),y(v.y),z(v.z)
{
}

template <class T, class tag>
inline TVec3D<T,tag>::TVec3D (const T& x1, const T& y1, const T& z1)
        :x(x1),y(y1),z(z1)
{
}


template <class T, class tag>
inline bool TVec3D<T,tag>::operator==(const TVec3D<T,tag> &right) const
{
    return x==right.x && y==right.y && z==right.z;
}

template <class T, class tag>
inline bool TVec3D<T,tag>::operator!=(const TVec3D<T,tag> &right) const
{
    return x!=right.x || y!=right.y || z!=right.z;
}



template <class T, class tag>
inline const TVec3D<T,tag>& TVec3D<T,tag>::operator = (const TVec3D<T,tag>& v)
{
    if(this!=&v){
        x=v.x;
        y=v.y;
        z=v.z;
    }
    return *this;
}

template <class T, class tag>
inline const TVec3D<T,tag>& TVec3D<T,tag>::operator += (const TVec3D<T,tag>& v)
{
    x+=v.x;
    y+=v.y;
    z+=v.z;
    return *this;
}

template <class T, class tag>
inline const TVec3D<T,tag>& TVec3D<T,tag>::operator -= (const TVec3D<T,tag>& v)
{
    x-=v.x;
    y-=v.y;
    z-=v.z;
    return *this;
}

template <class T, class tag>
inline const TVec3D<T,tag>& TVec3D<T,tag>::operator *= (const T& s)
{
    x*=s;
    y*=s;
    z*=s;
    return *this;
}

template <class T, class tag>
inline const TVec3D<T,tag>& TVec3D<T,tag>::operator /= (const T& s)
{
    x/=s;
    y/=s;
    z/=s;
    return *this;
}

template <class T, class tag>
inline TVec3D<T,tag> TVec3D<T,tag>::operator - () const
{
    return TVec3D<T,tag>(-x,-y,-z);
}

template <class T, class tag>
inline TVec3D<T,tag> TVec3D<T,tag>::operator - (const TVec3D<T,tag>& v) const
{
    return TVec3D<T,tag>(x-v.x,y-v.y,z-v.z);
}

template <class T, class tag>
inline T TVec3D<T,tag>::dot (const TVec3D<T,tag>& v) const
{
    return x*v.x+y*v.y+z*v.z;
}

template <class T, class tag>
inline bool TVec3D<T,tag>::isZero () const
{
    return x==0 && y==0 && z==0;
}

template <class T, class tag>
inline TVec3D<T,tag> TVec3D<T,tag>::operator + (const TVec3D<T,tag>& v) const
{
    return TVec3D<T,tag>(x+v.x,y+v.y,z+v.z);
}

template <class T, class tag>
inline TVec3D<T,tag> TVec3D<T,tag>::operator * (T s) const
{
    return TVec3D<T,tag>(s*x,s*y,s*z);
}

template <class T, class tag>
inline TVec3D<T,tag> TVec3D<T,tag>::operator / (T s) const
{
    return TVec3D<T,tag>(x/s,y/s,z/s);
}

// Parameterized Class TVec3D





template <class T, class tag>
TVec3D<T,tag> operator * (T s, const TVec3D<T,tag>& v)
{
    return TVec3D<T,tag>(s*v.x,s*v.y,s*v.z);
}

template<class T, class tag>
inline TVec3D<T,tag> cross(const TVec3D<T,tag>& a, const TVec3D<T,tag>& b)
{
    return TVec3D<T,tag>(a.y*b.z-a.z*b.y,b.x*a.z-b.z*a.x,a.x*b.y-a.y*b.x);
}

template<class T, class tag>
inline T operator*(const TVec3D<T,tag>& v1, const TVec3D<T,tag>& v2)
{
    return v1.dot(v2);
}


#endif
