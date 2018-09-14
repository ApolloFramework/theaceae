/*
  $Id: IncompressibleSAFlatPlate.cpp 26 2018-06-08 19:30:25Z ehowell $
  Some useful functions.  
  May want to push up to Camellia level.
 */

#include "Camellia.h" 
#include <math.h>

using namespace Camellia;
using namespace std;

/********************* Helper Functions ************************/
class RampFunction : public Function {
  private:
    int _t0;
    int _tf;
    double _dStep;  
    double _value;
  public:
  RampFunction( int t0, int tf ) : Function(0), _t0(t0), _tf(tf) {
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
  void UpdateStep(int step) {
    if (step >= _tf)
      {
	_value = 1.0;
      }
    else if ( step <= _t0)
      {
	_value = 0.0; 
      }
    else
      {
	_value = _dStep * (step - _t0);
      }
    cout << "Step: "<< step << " Value: " << _value << endl;
  }
  void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
    int numCells = values.dimension(0);
    int numPoints = values.dimension(1);

    

    const Intrepid::FieldContainer<double> *points = &(basisCache->getPhysicalCubaturePoints());
    for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
      for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
        values(cellIndex, ptIndex) = _value;
      }
    }
  }
};


class PowFunction : public Function {
  private:
    FunctionPtr _function;
    double _n;
  public:
    PowFunction(FunctionPtr function, double n) : Function(0), _function(function), _n(n) {}
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
      int numCells = values.dimension(0);
      int numPoints = values.dimension(1);

      _function->values(values, basisCache);

      const Intrepid::FieldContainer<double> *points = &(basisCache->getPhysicalCubaturePoints());
      for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
        for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
          values(cellIndex, ptIndex) = pow(values(cellIndex, ptIndex),_n);
        }
      }
    }
};

class AbsFunction : public Function {
  private:
    FunctionPtr _function;
  public:
    AbsFunction(FunctionPtr function) : Function(0), _function(function) {}
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
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
};

class EFunction : public Function {
  private:
    FunctionPtr _function;
  public:
    EFunction(FunctionPtr function) : Function(0), _function(function) {}
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
      int numCells = values.dimension(0);
      int numPoints = values.dimension(1);

      _function->values(values, basisCache);

      const Intrepid::FieldContainer<double> *points = &(basisCache->getPhysicalCubaturePoints());
      for (int cellIndex=0; cellIndex<numCells; cellIndex++) {
        for (int ptIndex=0; ptIndex<numPoints; ptIndex++) {
          values(cellIndex, ptIndex) = exp(values(cellIndex, ptIndex));
        }
      }
    }
};

class GtrZeroFunction : public Function {
  private:
    FunctionPtr _function;
  public:
    GtrZeroFunction(FunctionPtr function) : Function(0), _function(function) {}
    void values(Intrepid::FieldContainer<double> &values, BasisCachePtr basisCache) {
      int numCells = values.dimension(0);
      int numPoints = values.dimension(1);

      _function->values(values, basisCache);

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
};


class distanceToPlate : public SimpleFunction<double>
{
  double _plateX0; //x value of leading edge of the plate
  double _plateY; //y value of plate (assumes a horziontal plate) 
public:
  distanceToPlate(double plateX0, double plateY)
    {
      _plateX0=plateX0;
      _plateY =plateY;
    }
  double value(double x, double y)
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
};


class yValue : public SimpleFunction<double>
{
public:
  yValue() {}
  double value(double x, double y)
  {
    return (y);
  }
};
