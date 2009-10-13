


#ifndef Space2D_h
#define Space2D_h 1

#include "merlin_config.h"


// VectorTags
#include "EuclideanGeometry/VectorTags.h"



//	Template class for generating two-dimensional vectors
//	of the form (x,y). The template parameters
//	are the type of data storage (double, float etc.) and an
//	additional tag value. The tag class acts as a type-safe
//	mechanism for distinguising different types of three
//	vector (see Point2D and Vector2D).

template <class T, class tag>
class TVec2D
{
public:
    //	Default constructor. Data is not initialised.
    TVec2D ()
    {
    }

    //	Copy constructor.
    TVec2D (const TVec2D& v)
            : x(v.x),y(v.y)
    {
    }

    //	Explicit constructor from the three vector components.
    TVec2D (const T& x1, const T& y1)
            : x(x1),y(y1)
    {
    }


    //	Copy assignment
    const TVec2D& operator = (const TVec2D& v)
    {
        x=v.x;y=v.y;return *this;
    }

    //	Arithmetic assignment.
    const TVec2D& operator += (const TVec2D& v)
    {
        x+=v.x;y+=v.y;return *this;
    }

    const TVec2D& operator -= (const TVec2D& v)
    {
        x-=v.x;y-=v.y;return *this;
    }

    const TVec2D& operator *= (const T& s)
    {
        x*=s;y*=s;return *this;
    }

    const TVec2D& operator /= (const T& s)
    {
        x/=s;y/=s;return *this;
    }

    TVec2D operator - () const
    {
        return TVec2D(-x,-y);
    }

    TVec2D operator - (const TVec2D& v) const
    {
        return TVec2D(x-v.x,y-v.y);
    }

    //	Dot (inner) product.
    T dot (const TVec2D& v) const
    {
        return x*v.x+y*v.y;
    }

    //	Return true if all components are zero.
    bool isZero () const
    {
        return x==0&&y==0;
    }

    TVec2D operator + (const TVec2D& v) const
    {
        return TVec2D(x+v.x,y+v.y);
    }

    TVec2D operator * (T s) const
    {
        return TVec2D(s*x,s*y);
    }

    friend TVec2D operator * (T s, const TVec2D& v)
    {
        return TVec2D(s*v.x,s*v.y);
    }

    TVec2D operator / (T s) const
    {
        return TVec2D(x/s,y/s);
    }

    //	dot product.
    T operator * (const TVec2D& v) const
    {
        return dot(v);
    }

    // Data Members for Class Attributes

    T x;

    T y;

protected:
private:
private:
};



typedef TVec2D< double,PointTag > Point2D;



typedef TVec2D< double,VectorTag > Vector2D;

// Parameterized Class TVec2D



template class TVec2D< double,PointTag >;
template class TVec2D< double,VectorTag >;


#endif
