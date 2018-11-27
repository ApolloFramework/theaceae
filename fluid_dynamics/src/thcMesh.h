/*
  $Id: IncompressibleSAFlatPlate.cpp 26 2018-06-08 19:30:25Z ehowell $
  This code is implaments the incompressible Navier-Stokes equation with
  a the Spalart-Allmaras Turbulence Model. 
  It then solves for the the flow over a flat plate.

  This model uses a conservative form for the navier stokes equations, and 
  sa-working variable equation. It also uses a modified inner product space.
 */

#ifndef THEACEAE_THCMESH_H
#define THEACEAE_THCMESH_H

#include "Camellia.h" 
#include "Teuchos_ParameterList.hpp" 
#include <math.h>

using namespace Camellia;
using namespace std;

namespace Theaceae
{

MeshTopologyPtr thcCreateRectMesh(ParameterList meshPL );

/** \brief Create packed rectangular mesh*/
MeshTopologyPtr RectMeshPack(vector<double> dimensions, vector<int> elementCounts, vector<double>x0, double xPack, double xTau, double yBeta);
/** \brief Create packed rectangular mesh*/
MeshTopologyPtr RectMeshPack2(vector<double> dimensions, vector<int> elementCounts, vector<double>x0, double xPack, double xBetaL, double xBetaR, int nXL, double yBeta);

}
#endif
