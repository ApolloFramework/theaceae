/*
  $Id: IncompressibleSAFlatPlate.cpp 26 2018-06-08 19:30:25Z ehowell $
  This code is implaments the incompressible Navier-Stokes equation with
  a the Spalart-Allmaras Turbulence Model. 
  It then solves for the the flow over a flat plate.

  This model uses a conservative form for the navier stokes equations, and 
  sa-working variable equation. It also uses a modified inner product space.
 */

#include "Camellia.h" 
#include "Teuchos_ParameterList.hpp" 
#include <math.h>

using namespace Camellia;
using namespace std;


/********************* Create a packed rectangular mesh ************************/
MeshTopologyPtr RectMeshPack(vector<double> dimensions, vector<int> elementCounts, vector<double>x0, double xPack, double xTau, double yBeta)
{
// assume 2D for now
  int spaceDim = 2;
  int numElements = elementCounts[0] * elementCounts[1];

// determine the boundary of the domain
  double xLeft = x0[0];
  double xRight = xLeft + dimensions[0];
  double yBottom = x0[1];
  double yTop = yBottom + dimensions[1];
  CellTopoPtr topo;
  topo = Camellia::CellTopology::quad();
  vector< CellTopoPtr > cellTopos(numElements, topo);

  Intrepid::FieldContainer<double> quadBoundaryPoints(4,spaceDim);
  quadBoundaryPoints(0,0) = xLeft;
  quadBoundaryPoints(0,1) = yBottom;
  quadBoundaryPoints(1,0) = xRight;
  quadBoundaryPoints(1,1) = yBottom;
  quadBoundaryPoints(2,0) = xRight;
  quadBoundaryPoints(2,1) = yTop;
  quadBoundaryPoints(3,0) = xLeft;
  quadBoundaryPoints(3,1) = yTop;

  vector<vector<double> > vertices;
  vector< vector<IndexType> > allElementVertices;
  vector< vector<int> > vertexIndices(elementCounts[0]+1, vector<int>(elementCounts[1]+1));

  double dxBar = 1.0/(elementCounts[0]);
  double dyBar = 1.0/(elementCounts[1]);
  double bLnNum = 1.0 + (std::exp( xTau)-1.0)*(xPack-xLeft)/(xRight-xLeft);
  double bLnDen = 1.0 + (std::exp(-xTau)-1.0)*(xPack-xLeft)/(xRight-xLeft);
  double xB = 1.0/(2.0*xTau) * std::log(bLnNum/bLnDen);
  double yBp1 = yBeta + 1.0;
  double yBm1 = yBeta - 1.0;

  double xTemp=0;
  for (int i = 0 ; i<=elementCounts[0]; i++)
  {
    if (i ==0)
      {
	xTemp= xLeft;
      }
    else if (i==elementCounts[0])
      {
	xTemp = xRight;
      }
    else
      {
	xTemp = (xPack-xLeft) * (1.0+std::sinh(xTau * (i*dxBar-xB))/std::sinh(xTau*xB))+xLeft;
      }
    for (int j = 0 ; j<=elementCounts[1]; j++)
      {
      vertexIndices[i][j] = vertices.size();
      vector<double> vertex(spaceDim);
      vertex[0] = xTemp;
      if (j ==0)
	{
	  vertex[1] = yBottom;
	}
      else if (j==elementCounts[1])
	{
	  vertex[1] = yTop;
	}
      else
	{
	  vertex[1] = dimensions[1] * (yBp1- yBm1*(std::pow(yBp1/yBm1,1-j*dyBar)))/
	    (1.0 + std::pow(yBp1/yBm1,1-j*dyBar))+yBottom;
	}
      vertices.push_back(vertex);
      }
  }
  for (int i = 0 ; i<elementCounts[0]; i++)
  {
    for (int j = 0 ; j<elementCounts[1]; j++)
      {

	vector<IndexType> elemVertices;
        elemVertices.push_back(vertexIndices[i][j]);
        elemVertices.push_back(vertexIndices[i+1][j]);
        elemVertices.push_back(vertexIndices[i+1][j+1]);
        elemVertices.push_back(vertexIndices[i][j+1]);
        allElementVertices.push_back(elemVertices);
      }
  }
  MeshGeometryPtr geometry = Teuchos::rcp( new MeshGeometry(vertices, allElementVertices, cellTopos));
  return Teuchos::rcp( new MeshTopology(geometry) );
}


MeshTopologyPtr RectMeshPack2(vector<double> dimensions, vector<int> elementCounts, vector<double>x0, double xPack, double xBetaL, double xBetaR, int nXL, double yBeta)
{
// assume 2D for now
  int spaceDim = 2;
  int numElements = elementCounts[0] * elementCounts[1];


// determine the boundary of the domain
  double xLeft = x0[0];
  double xRight = xLeft + dimensions[0];
  double yBottom = x0[1];
  double yTop = yBottom + dimensions[1];
  CellTopoPtr topo;
  topo = Camellia::CellTopology::quad();
  vector< CellTopoPtr > cellTopos(numElements, topo);

  Intrepid::FieldContainer<double> quadBoundaryPoints(4,spaceDim);
  quadBoundaryPoints(0,0) = xLeft;
  quadBoundaryPoints(0,1) = yBottom;
  quadBoundaryPoints(1,0) = xRight;
  quadBoundaryPoints(1,1) = yBottom;
  quadBoundaryPoints(2,0) = xRight;
  quadBoundaryPoints(2,1) = yTop;
  quadBoundaryPoints(3,0) = xLeft;
  quadBoundaryPoints(3,1) = yTop;

  vector<vector<double> > vertices;
  vector< vector<IndexType> > allElementVertices;
  vector< vector<int> > vertexIndices(elementCounts[0]+1, vector<int>(elementCounts[1]+1));

  double dyBar = 1.0/(elementCounts[1]);

  int nXR = elementCounts[0] - nXL;

  double dxBarL = 1.0 /nXL;
  double dxBarR = 1.0 /nXR;
  double xBp1L = xBetaL + 1.0;
  double xBm1L = xBetaL - 1.0;
  double xBp1R = xBetaR + 1.0;
  double xBm1R = xBetaR - 1.0;

  double yBp1 = yBeta + 1.0;
  double yBm1 = yBeta - 1.0;

  double xTemp=0;
  for (int i = 0 ; i<=elementCounts[0]; i++)
  {
    if (i ==0)
      {
	xTemp= xLeft;
      }
    else if (i==elementCounts[0])
      {
	xTemp = xRight;
      }

    else if (i == nXL)
      {
	xTemp = xPack;
      }
    else if (i < nXL)
      {
	xTemp = (xLeft-xPack) * (xBp1L- xBm1L*(std::pow(xBp1L/xBm1L,1-(nXL-i)*dxBarL)))/
	    (1.0 + std::pow(xBp1L/xBm1L,1-(nXL-i)*dxBarL))+xPack;
      }
    else
      {
	xTemp = (xRight-xPack) * (xBp1R- xBm1R*(std::pow(xBp1R/xBm1R,1-(i-nXL)*dxBarR)))/
	    (1.0 + std::pow(xBp1R/xBm1R,1-(i-nXL)*dxBarR))+xPack;
      }
    for (int j = 0 ; j<=elementCounts[1]; j++)
      {
      vertexIndices[i][j] = vertices.size();
      vector<double> vertex(spaceDim);
      vertex[0] = xTemp;
      if (j ==0)
	{
	  vertex[1] = yBottom;
	}
      else if (j==elementCounts[1])
	{
	  vertex[1] = yTop;
	}
      else
	{
	  vertex[1] = dimensions[1] * (yBp1- yBm1*(std::pow(yBp1/yBm1,1-j*dyBar)))/
	    (1.0 + std::pow(yBp1/yBm1,1-j*dyBar))+yBottom;
	}
      vertices.push_back(vertex);
      }
  }
  for (int i = 0 ; i<elementCounts[0]; i++)
  {
    for (int j = 0 ; j<elementCounts[1]; j++)
      {

	vector<IndexType> elemVertices;
        elemVertices.push_back(vertexIndices[i][j]);
        elemVertices.push_back(vertexIndices[i+1][j]);
        elemVertices.push_back(vertexIndices[i+1][j+1]);
        elemVertices.push_back(vertexIndices[i][j+1]);
        allElementVertices.push_back(elemVertices);
      }
  }
  MeshGeometryPtr geometry = Teuchos::rcp( new MeshGeometry(vertices, allElementVertices, cellTopos));
  return Teuchos::rcp( new MeshTopology(geometry) );
}
