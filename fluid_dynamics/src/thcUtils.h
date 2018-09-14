/*
  $Id: IncompressibleSAFlatPlate.cpp 26 2018-06-08 19:30:25Z ehowell $
  Some useful functions.  
  May want to push up to Camellia level.
 */

#ifndef THEACEAE_THCUTILS_H
#define THEACEAE_THCUTILS_H

#include "Camellia.h" 
#include <math.h>

using namespace Camellia;
using namespace std;

namespace Theaceae
{

/** \brief Function for time ramp of fields
 */
class RampFunction : public Function {
  public:
    /** \brief Actual function */
    RampFunction( int t0, int tf );
    /** \brief Update to step number   */
    void UpdateStep(int step);
    /** \brief values */
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache);
};


/** \brief Function for raising the fields to a power at each point
 */
class PowFunction : public Function {
  public:
    /** \brief Actual function */
    PowFunction(FunctionPtr function, double n);
    /** \brief values */
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache);
};

/** \brief Function for taking the absolute value of a field at each point
 */
class AbsFunction : public Function {
  public:
    /** \brief Actual function */
    AbsFunction(FunctionPtr function);
    /** \brief values */
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache);
};

/** \brief Function for taking field at each point
 */
class EFunction : public Function {
  public:
    /** \brief Actual function */
    EFunction(FunctionPtr function);
    /** \brief values */
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache);
};

/** \brief Function for taking field at each point
 */
class GtrZeroFunction : public Function {
  public:
    /** \brief Actual function */
    GtrZeroFunction(FunctionPtr function);
    /** \brief values */
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache);
};

/** \brief Calculate distance to plate
 */
class distanceToPlate : public SimpleFunction<double>
{
  public:
    /** \brief Actual function */
    distanceToPlate(double plateX0, double plateY);
    /** \brief values */
    double value(double x, double y);
};

/** \brief Calculate yvalue
 */
class yValue : public SimpleFunction<double>
{
public:
  /** \brief values */
  double value(double x, double y);
};

}
#endif
