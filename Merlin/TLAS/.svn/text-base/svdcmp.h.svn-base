/////////////////////////////////////////////////////////////////////////
//
// Merlin C++ Class Library for Charged Particle Accelerator Simulations
//  
// Class library version 3 (2004)
// 
// Copyright: see Merlin/copyright.txt
//
// Last CVS revision:
// $Date: 2004/12/13 08:38:55 $
// $Revision: 1.2 $
// 
/////////////////////////////////////////////////////////////////////////

#ifndef  _MSC_VER
#define _MAX(a,b) std::max((a),(b))
#endif

template<class T,class V>
Vector<T>& svbksb(const Matrix<T>& u, const Vector<T>& w, const Matrix<T>& v, const V& b, Vector<T>& x)
{
    const int m = b.size();
    const int n = u.ncols();
    int jj,j,i;
    std::vector<T> tmp(n);

    for(j=0;j<n;j++) {
        T s=0.0;
        if (w[j]) {
            for (i=0;i<m;i++)
                s += u[i][j]*b[i];
            s /= w[j];
        }
        tmp[j]=s;
    }
    for(j=0;j<n;j++) {
        T s=0.0;
        for (jj=0;jj<n;jj++)
            s += v[j][jj]*tmp[jj];
        x[j]=s;
    }
    return x;
}

template<class T>
inline T PYTHAG(T a, T b)
{
    a = fabs(a);
    b = fabs(b);
    if(a>b) {
        T c=b/a;
        return a*sqrt(1+c*c);
    }
    else if(b!=0) {
        T c=a/b;
        return b*sqrt(1+c*c);
    }
    else
        return T(0);
}

template<class T> inline T SIGN(T a, T b) { return b>=T(0) ? fabs(a) : -fabs(a); }

template<class T>
void svdcmp(Matrix<T>& a, Vector<T>& w, Matrix<T>& v)
{
    const int m = a.nrows();
    const int n = a.ncols();
    assert(n<=m);

    int flag,i,its,j,jj,k,l,nm;
    T c,f,h,s,x,y,z;
    T anorm=0.0,g=0.0,scale=0.0;

    std::vector<T> rv1(n);

    for (i=0;i<n;i++) {
        l=i+1;
        rv1[i]=scale*g;
        g=s=scale=0.0;
        if (i < m) {
            for (k=i;k<m;k++) scale += fabs(a[k][i]);
            if (scale) {
                for (k=i;k<m;k++) {
                    a[k][i] /= scale;
                    s += a[k][i]*a[k][i];
                }
                f=a[i][i];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][i]=f-g;
                if (i!=n-1) {
                    for (j=l;j<n;j++) {
                        for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
                        f=s/h;
                        for (k=i;k<m;k++) a[k][j] += f*a[k][i];
                    }
                }
                for (k=i;k<m;k++) a[k][i] *= scale;
            }
        }
        w[i]=scale*g;
        g=s=scale=0.0;
        if (i < m && i != n-1) {
            for (k=l;k<n;k++) scale += fabs(a[i][k]);
            if (scale) {
                for (k=l;k<n;k++) {
                    a[i][k] /= scale;
                    s += a[i][k]*a[i][k];
                }
                f=a[i][l];
                g = -SIGN(sqrt(s),f);
                h=f*g-s;
                a[i][l]=f-g;
                for (k=l;k<n;k++) rv1[k]=a[i][k]/h;
                if (i != m-1) {
                    for (j=l;j<m;j++) {
                        for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
                        for (k=l;k<n;k++) a[j][k] += s*rv1[k];
                    }
                }
                for (k=l;k<n;k++) a[i][k] *= scale;
            }
        }
        anorm=_MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
    }
    for (i=n-1;i>=0;i--) {
        if (i < n-1) {
            if (g) {
                for (j=l;j<n;j++)
                    v[j][i]=(a[i][j]/a[i][l])/g;
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
                    for (k=l;k<n;k++) v[k][j] += s*v[k][i];
                }
            }
            for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
        }
        v[i][i]=1.0;
        g=rv1[i];
        l=i;
    }
    for (i=n-1;i>=0;i--) {
        l=i+1;
        g=w[i];
        if (i < n)
            for (j=l;j<n;j++) a[i][j]=0.0;
        if (g) {
            g=1.0/g;
            if (i != n-1) {
                for (j=l;j<n;j++) {
                    for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
                    f=(s/a[i][i])*g;
                    for (k=i;k<m;k++) a[k][j] += f*a[k][i];
                }
            }
            for (j=i;j<m;j++) a[j][i] *= g;
        } else {
            for (j=i;j<m;j++) a[j][i]=0.0;
        }
        ++a[i][i];
    }
    for (k=n-1;k>=0;k--) {
        for (its=1;its<=30;its++) {
            flag=1;
            for (l=k;l>=0;l--) {
                nm=l-1;
                if (fabs(rv1[l])+anorm == anorm) {
                    flag=0;
                    break;
                }
                if (fabs(w[nm])+anorm == anorm) break;
            }
            if (flag) {
                c=0.0;
                s=1.0;
                for (i=l;i<=k;i++) {
                    f=s*rv1[i];
                    if (fabs(f)+anorm != anorm) {
                        g=w[i];
                        h=PYTHAG(f,g);
                        w[i]=h;
                        h=1.0/h;
                        c=g*h;
                        s=(-f*h);
                        for (j=0;j<m;j++) {
                            y=a[j][nm];
                            z=a[j][i];
                            a[j][nm]=y*c+z*s;
                            a[j][i]=z*c-y*s;
                        }
                    }
                }
            }
            z=w[k];
            if (l == k) {
                if (z < 0.0) {
                    w[k] = -z;
                    for (j=0;j<n;j++) v[j][k]=(-v[j][k]);
                }
                break;
            }
            if (its == 30)
                throw  ConvergenceFailure();

            x=w[l];
            nm=k-1;
            y=w[nm];
            g=rv1[nm];
            h=rv1[k];
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g=PYTHAG(f,1.0);
            f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
            c=s=1.0;
            for (j=l;j<=nm;j++) {
                i=j+1;
                g=rv1[i];
                y=w[i];
                h=s*g;
                g=c*g;
                z=PYTHAG(f,h);
                rv1[j]=z;
                c=f/z;
                s=h/z;
                f=x*c+g*s;
                g=g*c-x*s;
                h=y*s;
                y=y*c;
                for (jj=0;jj<n;jj++) {
                    x=v[jj][j];
                    z=v[jj][i];
                    v[jj][j]=x*c+z*s;
                    v[jj][i]=z*c-x*s;
                }
                z=PYTHAG(f,h);
                w[j]=z;
                if (z) {
                    z=1.0/z;
                    c=f*z;
                    s=h*z;
                }
                f=(c*g)+(s*y);
                x=(c*y)-(s*g);
                for (jj=0;jj<m;jj++) {
                    y=a[jj][j];
                    z=a[jj][i];
                    a[jj][j]=y*c+z*s;
                    a[jj][i]=z*c-y*s;
                }
            }
            rv1[l]=0.0;
            rv1[k]=f;
            w[k]=x;
        }
    }
}
