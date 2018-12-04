/*
  $Id: IncompressibleSAFlatPlate.cpp 26 2018-06-08 19:30:25Z ehowell $
  Some useful functions.  
  May want to push up to Camellia level.
 */

#ifndef THEACEAE_THCUTILS_H
#define THEACEAE_THCUTILS_H

#include "Camellia.h" 
#include <math.h>
#include <functional>

// using namespace Camellia;
using namespace std;

namespace Theaceae
{

/** \brief Function for time ramp of fields
 */
template <typename Scalar>
class RampFunction : public Camellia::SimpleFunction<Scalar> {
  private:
    int _t0;
    int _tf;
    double _dStep;  
    Scalar _value;
  public:
    /** \brief Actual function */
    RampFunction( int t0, int tf );
    /** \brief Update to step number   */
    void UpdateStep(int step);
    /** \brief values */
    void values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache);
};


/** \brief Function for raising the fields to a power at each point
 */
template <typename Scalar>
class PowFunction : public Camellia::SimpleFunction<Scalar> {
  private:
    Camellia::TFunctionPtr<Scalar> _function{nullptr};
    double _n;
  public:
    /** \brief Actual function */
    PowFunction(Camellia::TFunctionPtr<Scalar> function, double n) {
      _function = function;
      _n = n;
    };
    /** \brief values */
    void values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache);
};

/** \brief Function for taking the absolute value of a field at each point
 */
template <typename Scalar>
class AbsFunction : public Camellia::SimpleFunction<Scalar> {
  private:
    Camellia::TFunctionPtr<Scalar> _function;
  public:
    /** \brief Actual function */
    AbsFunction(Camellia::TFunctionPtr<Scalar> function)  {
      _function = function;
    };
    /** \brief values */
    void values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache);
};

/** \brief Function for taking field at each point
 */
template <typename Scalar>
class EFunction : public Camellia::SimpleFunction<Scalar> {
  private:
    Camellia::TFunctionPtr<Scalar> _function;
  public:
    /** \brief Actual function */
    EFunction(Camellia::TFunctionPtr<Scalar> function) {
      _function = function;
    };
    /** \brief values */
    void values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache);
};

/** \brief Function for taking field at each point
 */
template <typename Scalar>
class GtrZeroFunction : public Camellia::SimpleFunction<Scalar> {
  private:
    Camellia::TFunctionPtr<Scalar> _function;
  public:
    /** \brief Actual function */
    GtrZeroFunction(Camellia::TFunctionPtr<Scalar> function) {
      _function = function;
    };
    /** \brief values */
    void values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache);
};

/** \brief Calculate distance to plate
 */
class distanceToPlate : public Camellia::SimpleFunction<double>
{
  double _plateX0; //x value of leading edge of the plate
  double _plateY; //y value of plate (assumes a horziontal plate) 
  public:
    /** \brief Actual function */
    distanceToPlate(double plateX0, double plateY) {
      _plateX0=plateX0;
      _plateY =plateY;
    };
    /** \brief values */
    double value(double x, double y);
};

/** \brief Calculate yvalue
 */
class yValue : public Camellia::SimpleFunction<double> {
public:
  /** \brief values */
  double value(double x, double y);
};

}
#endif
