#include "nonlin.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <iostream>
using namespace std;

namespace cpt {

static    int sign_of_func;
static    double (*user_func)(const double);
static    double func(const double x)
{
    return sign_of_func * (*user_func)(x);
}

inline void shft3(double &a, double &b, double &c, const double d)
{
    a=b; b=c; c=d;
}

inline double SIGN(const double& a, const double& b)
{
    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}

double find_maximum(const double a, const double b,
    double f(const double), double tol, double& xmin)
{
    sign_of_func = -1;
    user_func = f;
    double ax = a, bx = b, cx, fa, fb, fc;
    mnbrak(ax, bx, cx, fa, fb, fc, func);
    return -golden(ax, bx, cx, func, tol, xmin);
}

double find_minimum(const double a, const double b,
    double f(const double), double tol, double& xmin)
{
    double ax = a, bx = b, cx, fa, fb, fc;
    mnbrak(ax, bx, cx, fa, fb, fc, f);
    return golden(ax, bx, cx, f, tol, xmin);
}

void mnbrak(double &ax, double &bx, double &cx,
      double &fa, double &fb, double &fc,
      double func(const double))
{
    const double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
    double ulim,u,r,q,fu;

    fa=func(ax);
    fb=func(bx);
    if (fb > fa) {
        swap(ax,bx);
        swap(fb,fa);
    }
    cx=bx+GOLD*(bx-ax);
    fc=func(cx);
    while (fb > fc) {
        r=(bx-ax)*(fb-fc);
        q=(bx-cx)*(fb-fa);
        u=bx-((bx-cx)*q-(bx-ax)*r)/
            (2.0*SIGN(max(abs(q-r),TINY),q-r));
        ulim=bx+GLIMIT*(cx-bx);
        if ((bx-u)*(u-cx) > 0.0) {
            fu=func(u);
            if (fu < fc) {
                ax=bx;
                bx=u;
                fa=fb;
                fb=fu;
                return;
            } else if (fu > fb) {
                cx=u;
                fc=fu;
                return;
            }
            u=cx+GOLD*(cx-bx);
            fu=func(u);
        } else if ((cx-u)*(u-ulim) > 0.0) {
            fu=func(u);
            if (fu < fc) {
                shft3(bx,cx,u,cx+GOLD*(cx-bx));
                shft3(fb,fc,fu,func(u));
            }
        } else if ((u-ulim)*(ulim-cx) >= 0.0) {
            u=ulim;
            fu=func(u);
        } else {
            u=cx+GOLD*(cx-bx);
            fu=func(u);
        }
        shft3(ax,bx,cx,u);
        shft3(fa,fb,fc,fu);
    }

}

inline void shft2(double &a, double &b, const double c)
{
    a=b;  b=c;
}

double golden(const double ax, const double bx, const double cx,
              double f(const double), const double tol, double &xmin)
{
    const double R=0.61803399,C=1.0-R;
    double f1,f2,x0,x1,x2,x3;

    x0=ax;
    x3=cx;
    if (abs(cx-bx) > abs(bx-ax)) {
        x1=bx;
        x2=bx+C*(cx-bx);
    } else {
        x2=bx;
        x1=bx-C*(bx-ax);
    }
    f1=f(x1);
    f2=f(x2);
    while (abs(x3-x0) > tol*(abs(x1)+abs(x2))) {
        if (f2 < f1) {
            shft3(x0,x1,x2,R*x2+C*x3);
            shft2(f1,f2,f(x2));
        } else {
            shft3(x3,x2,x1,R*x1+C*x0);
            shft2(f2,f1,f(x1));
        }
    }
    if (f1 < f2) {
        xmin=x1;
        return f1;
    } else {
        xmin=x2;
        return f2;
    }
}

void line_search(
    Matrix<double,1>& xold, double fold, Matrix<double,1>& g,
    Matrix<double,1>& p, Matrix<double,1>& x, double& f, double stpmax,
    bool& check, double (*func)(Matrix<double,1>&))
{
    const double ALF=1.0e-4, TOLX=std::numeric_limits<double>::epsilon();
    int i;
    double a,alam,alam2=0.0,alamin,b,disc,f2=0.0;
    double rhs1,rhs2,slope,sum,temp,test,tmplam;

    int n = xold.size();
    check = false;
    sum=0.0;
    for (i=0;i<n;i++) sum += p[i]*p[i];
    sum=std::sqrt(sum);
    if (sum > stpmax)
        for (i=0;i<n;i++) p[i] *= stpmax/sum;
    slope = 0.0;
    for (i=0;i<n;i++)
        slope += g[i]*p[i];
    if (slope >= 0.0) error("nonlin: Roundoff problem in line_search");
    test=0.0;
    for (i=0;i<n;i++) {
        temp=std::abs(p[i])/std::max(std::abs(xold[i]),1.0);
        if (temp > test) test=temp;
    }
    alamin=TOLX/test;
    alam=1.0;
    for (;;) {
        for (i=0;i<n;i++) x[i]=xold[i]+alam*p[i];
        f=func(x);
        if (alam < alamin) {
            for (i=0;i<n;i++) x[i]=xold[i];
            check = true;
            return;
        } else if (f <= fold+ALF*alam*slope) return;
        else {
            if (alam == 1.0)
                tmplam = -slope/(2.0*(f-fold-slope));
            else {
                rhs1 = f-fold-alam*slope;
                rhs2=f2-fold-alam2*slope;
                a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
                b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
                if (a == 0.0) tmplam = -slope/(2.0*b);
                else {
                    disc=b*b-3.0*a*slope;
                    if (disc < 0.0) tmplam=0.5*alam;
                    else if (b <= 0.0) tmplam=(-b+std::sqrt(disc))/(3.0*a);
                    else tmplam=-slope/(b+std::sqrt(disc));
                }
                if (tmplam > 0.5*alam)
                    tmplam=0.5*alam;
            }
        }
        alam2=alam;
        f2 = f;
        alam=std::max(tmplam,0.1*alam);
    }
}

void minimize_BFGS(
    Matrix<double,1>& p, const double gtol, int& iter, double& fret,
    double (*func)(Matrix<double,1>&), void (*dfunc)(Matrix<double,1>&,
    Matrix<double,1>&))
{
    const int ITMAX = 200;
    const double EPS = std::numeric_limits<double>::epsilon();
    const double TOLX = 4*EPS, STPMX = 100.0;
    bool check;
    int i,its,j;
    double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;

    int n = p.size();
    Matrix<double,1> dg(n),g(n),hdg(n),pnew(n),xi(n);
    Matrix<double,2> hessin(n, n);

    fp=func(p);
    dfunc(p,g);
    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) hessin[i][j]=0.0;
        hessin[i][i]=1.0;
        xi[i] = -g[i];
        sum += p[i]*p[i];
    }
    stpmax=STPMX*std::max(std::sqrt(sum),double(n));
    for (its=0;its<ITMAX;its++) {
        iter=its;
        line_search(p,fp,g,xi,pnew,fret,stpmax,check,func);
        fp = fret;
        for (i=0;i<n;i++) {
            xi[i]=pnew[i]-p[i];
            p[i]=pnew[i];
        }
        test=0.0;
        for (i=0;i<n;i++) {
            temp=std::abs(xi[i])/std::max(std::abs(p[i]),1.0);
            if (temp > test) test=temp;
        }
        if (test < TOLX)
            return;
        for (i=0;i<n;i++) dg[i]=g[i];
        dfunc(p,g);
        test=0.0;
        den=std::max(fret,1.0);
        for (i=0;i<n;i++) {
            temp=std::abs(g[i])*std::max(std::abs(p[i]),1.0)/den;
            if (temp > test) test=temp;
        }
        if (test < gtol)
            return;
        for (i=0;i<n;i++) dg[i]=g[i]-dg[i];
        for (i=0;i<n;i++) {
            hdg[i]=0.0;
            for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
        }
        fac=fae=sumdg=sumxi=0.0;
        for (i=0;i<n;i++) {
            fac += dg[i]*xi[i];
            fae += dg[i]*hdg[i];
            sumdg += dg[i]*dg[i];
            sumxi += xi[i]*xi[i];
        }
        if (fac > std::sqrt(EPS*sumdg*sumxi)) {
            fac=1.0/fac;
            fad=1.0/fae;
            for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
            for (i=0;i<n;i++) {
                for (j=i;j<n;j++) {
                    hessin[i][j] += fac*xi[i]*xi[j]
                        -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
                    hessin[j][i]=hessin[i][j];
                }
            }
        }
        for (i=0;i<n;i++) {
            xi[i]=0.0;
            for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
        }
    }
    error("nonlin: Too many iterations in minimize_BFGS");
}

} /* end namespace cpt */
