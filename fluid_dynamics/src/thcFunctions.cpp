/*
  $Id: IncompressibleSAFlatPlate.cpp 26 2018-06-08 19:30:25Z ehowell $
  Some useful functions.  
  May want to push up to Camellia level.
 */

#include "Camellia.h" 
#include <math.h>
#include "thcFunctions.h" 

// using namespace Camellia;
using namespace std;
using namespace Theaceae;

/********************* Helper Functions ************************/
template <typename Scalar>
RampFunction<Scalar>::RampFunction( int t0, int tf ) : _t0(t0), _tf(tf) {
  _dStep = 1.0/(1.0 *(_tf - _t0) );
  if (_tf < 0)
  {
    _value = 1.0;
  }
  else if (_t0 < 0)
  {
    _value = -_dStep * _t0;
  }
  else
  {
    _value = 0.0;
  }
}

template <typename Scalar>
void RampFunction<Scalar>::UpdateStep(int step) {
  if (step >= RampFunction::_tf)
  {
    RampFunction::_value = 1.0;
  }
  else if ( step <= RampFunction::_t0)
  {
    RampFunction::_value = 0.0; 
  }
  else
  {
    RampFunction::_value = RampFunction::_dStep * (step - RampFunction::_t0);
  }
  cout << "Step: "<< step << " Value: " << _value << endl;
}

template <typename Scalar>
void RampFunction<Scalar>::values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
      values(cellIndex, ptIndex) = RampFunction::_value;
    }
  }
}


template <typename Scalar>
void PowFunction<Scalar>::values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  _function->values(values, basisCache);

  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
      values(cellIndex, ptIndex) = pow(values(cellIndex, ptIndex),PowFunction::_n);
    }
  }
}

template <typename Scalar>
void AbsFunction<Scalar>::values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  _function->values(values, basisCache);

  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
      values(cellIndex, ptIndex) = abs(values(cellIndex, ptIndex));
    }
  }
}

template <typename Scalar>
void EFunction<Scalar>::values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  EFunction::_function->values(values, basisCache);

  const Intrepid::FieldContainer<Scalar> *points = &(basisCache->getPhysicalCubaturePoints());
  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
      values(cellIndex, ptIndex) = exp(values(cellIndex, ptIndex));
    }
  }
}

template <typename Scalar>
void GtrZeroFunction<Scalar>::values(Intrepid::FieldContainer<Scalar> &values, Camellia::BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  GtrZeroFunction::_function->values(values, basisCache);

  const Intrepid::FieldContainer<Scalar> *points = &(basisCache->getPhysicalCubaturePoints());
  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
  if (values(cellIndex, ptIndex)>0)
    {
      values(cellIndex, ptIndex) = 1.0;
    }
  else
    {
      values(cellIndex, ptIndex) = 0.0;
    }
    }
  }
}

/********************* More Specific Functions ************************/
double distanceToPlate::value(double x, double y)
{
  if (x < _plateX0)
    {
  return sqrt( pow(x-_plateX0,2) + pow(y-_plateY,2)); 
    }
  else
    {
  return (y-_plateY);
    }
}


double yValue::value(double x, double y)
{
  return (y);
}

namespace Theaceae
{
template class RampFunction<double>;
template class PowFunction<double>;
template class AbsFunction<double>;
template class EFunction<double>;
template class GtrZeroFunction<double>;
}
