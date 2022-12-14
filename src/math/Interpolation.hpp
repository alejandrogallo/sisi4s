#ifndef INTERPOLATION_DEFINED
#define INTERPOLATION_DEFINED

#include <math/MathFunctions.hpp>
#include <iostream>
#include <string>

// TODO: unify interpolation routines

namespace sisi4s {
template <typename F = double>
class Inter1D {
public:
  /*when you creat an instance of this class, you need to specify the size
    of the data set you are interpolating on, and the corresponding x and y
  */
  Inter1D(const int size, F *x, F *y) {
    N = size;
    xv = new F[N];
    yv = new F[N];
    y2d = new F[N];
    xv = x;
    yv = y;
  }

  /*x and y are arrays of the same length, containing the x data and
   the corresponding values of the function y. They two must be sorted
   from min to max according to x.N is the length of the two arrays. y1d0     is
   the 1st derivative of the function at the first point, and y1dn is the 1st
   derivative of the function at the end of the data. User can specify them two
   when the BoundaryCondition is set to "M" which stands for Manual, or it can
   be set to "N" which stands for "Natural".     This function returns an array
   of the second derivatives at each x points, which can be used in
   CubicSpline_getValue() function to inquire the value at certain point(s).
  */
  void cubicSpline(F y1d0 = 0.,
                   F y1dn = 0.,
                   const std::string boundaryCondition = "M") {
    double u[N];
    double qn(0.), un(0.), sig(0.), p(0.);
    if (boundaryCondition == "N") {
      LOG(1, "Interpolation::cubicSpline") << "Natural" << std::endl;
      y2d[0] = 0.;
      u[0] = 0.;
    } else if (boundaryCondition == "M") {
      LOG(1, "Interpolation::cubicSpline") << "Manual " << y1d0 << std::endl;
      y2d[0] = -0.5; // does not really matter, T.B.C
      u[0] =
          (3. / (xv[1] - xv[0])) * (((yv[1] - yv[0])) / (xv[1] - xv[0]) - y1d0);
      qn = 0.5;
      un = (3. / (xv[N - 1] - xv[N - 2]))
         * ((yv[N - 1] - yv[N - 2]) / (xv[N - 1] - xv[N - 2]) - y1dn);
    }
    for (int i(1); i < N - 1; i++) {
      sig = (xv[i] - xv[i - 1]) / (xv[i + 1] - xv[i - 1]);
      p = sig * y2d[i - 1] + 2.;
      y2d[i] = (sig - 1.) / p;
      u[i] = (6.
                  * ((yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i])
                     - (yv[i] - yv[i - 1]) / (xv[i] - xv[i - 1]))
                  / (xv[i + 1] - xv[i - 1])
              - sig * u[i - 1])
           / p;
    }
    y2d[N - 1] = (un - qn * u[N - 2]) / (qn * y2d[N - 2] + 1.);
    for (int t = N - 2; t >= 0; t--) y2d[t] = y2d[t] * y2d[t + 1] + u[t];
  }
  // The max(x) should not larger than max(xv).
  F getValue(F xinterp) {
    // LOG(1,"getValue") << "n= " << n << std::endl;
    // int n=1;
    double yinterp(0.);
    double h, a, b;
    int klow, khigh, k;
    // y= new double[n];
    // for (int i(0); i<n; ++i){
    klow = 0;
    khigh = N - 1;
    while ((khigh - klow) > 1) {
      k = (khigh + klow) / 2;
      if (xv[k - 1] > xinterp) khigh = k;
      else klow = k;
    }
    h = xv[khigh - 1] - xv[klow - 1];
    if (abs(h) < 1e-12) {
      LOG(1, "CubicSpline_getValue") << "Same x points in x data!"
                                     << "STOP!" << std::endl;
      exit(0);
    }
    a = (xv[khigh - 1] - xinterp) / h;
    b = (xinterp - xv[klow - 1]) / h;
    yinterp =
        a * yv[klow - 1] + b * yv[khigh - 1]
        + ((a * a * a - a) * y2d[klow - 1] + (b * b * b - b) * y2d[khigh - 1])
              * (h * h) / 6.;
    //}
    return yinterp;
  }
  int N;
  F *xv, *yv, *y2d;
};
} // namespace sisi4s

#endif
