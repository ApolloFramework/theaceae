/*
  $Id: IncompressibleSAFlatPlate.cpp 26 2018-06-08 19:30:25Z ehowell $
  Some useful functions.  
  May want to push up to Camellia level.
 */

#include "Camellia.h" 
#include <math.h>
#include "thcUtils.h" 

using namespace Camellia;
using namespace std;
using namespace Theaceae;

/********************* Helper Functions ************************/
RampFunction::RampFunction( int t0, int tf ) : Function(0), _t0(t0), _tf(tf) {
  RampFunction::_dStep = 1.0/(1.0 *(RampFunction::_tf - RampFunction::_t0) );
  if (RampFunction::_tf < 0)
  {
    RampFunction::_value = 1.0;
  }
  else if (RampFunction::_t0 < 0)
  {
    RampFunction::_value = -RampFunction::_dStep * RampFunction::_t0;
  }
  else
  {
    RampFunction::_value = 0.0;
  }
}
void RampFunction::UpdateStep(int step) {
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
void RampFunction::values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  const Intrepid::FieldContainer<double> *points = &(basisCache->getPhysicalCubaturePoints());
  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
      values(cellIndex, ptIndex) = RampFunction::_value;
    }
  }
}


void PowFunction::values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  PowFunction::_function->values(values, basisCache);

  const Intrepid::FieldContainer<double> *points = &(basisCache->getPhysicalCubaturePoints());
  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
      values(cellIndex, ptIndex) = pow(values(cellIndex, ptIndex),PowFunction::_n);
    }
  }
}

void AbsFunction::values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  _function->values(values, basisCache);

  const Intrepid::FieldContainer<double> *points = &(basisCache->getPhysicalCubaturePoints());
  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
      values(cellIndex, ptIndex) = abs(values(cellIndex, ptIndex));
    }
  }
}

void EFunction::values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  EFunction::_function->values(values, basisCache);

  const Intrepid::FieldContainer<double> *points = &(basisCache->getPhysicalCubaturePoints());
  for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
    for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
      values(cellIndex, ptIndex) = exp(values(cellIndex, ptIndex));
    }
  }
}

void GtrZeroFunction::values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
  int numCells = values.dimension(0);
  int numPoints = values.dimension(1);

  GtrZeroFunction::_function->values(values, basisCache);

  const Intrepid::FieldContainer<double> *points = &(basisCache->getPhysicalCubaturePoints());
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
