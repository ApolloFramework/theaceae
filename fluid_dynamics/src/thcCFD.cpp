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
#include "thcCfdInput.h" 
#include "thcUtils.h" 
#include "thcMesh.h" 
#include <math.h>

using namespace std;
using namespace Teuchos;
using namespace Camellia;
using namespace Theaceae;

/********************* Main driver program ************************/

int main( int argc, char *argv[] )
{

  GlobalMPISession mpiSession(&argc, &argv); // initalize mpi
  int rank = MPIWrapper::CommWorld()->MyPID();

  ParameterList  thcPL = thcCFDInput(argc, argv );

/******************* Set up 2D Flat Plate ***************************/

#if 0
  double xDim = plateLength + plateStart;
  vector<double> dims = {xDim, yDim}; // dimensions of the mesh
  vector<int> meshDims = {mx,my}; 
  vector<double> x0 = {-plateStart, 0.0}; // mesh is s.t. plate begins at x=0
  MeshTopologyPtr meshTopo;
  if (packMesh)
    {
//      meshTopo = RectMeshPack(dims, meshDims, x0, 0.0, xTau, yBeta);
      meshTopo = RectMeshPack2(dims, meshDims, x0, 0.0, xBetaL, xBetaR, nXL, yBeta);
    }
  else
    {
      meshTopo = MeshFactory::rectilinearMeshTopology(dims, meshDims, x0);
    }


/******************* Define variables for calculation ***************************/
  VarFactoryPtr vf = VarFactory::varFactory();
  VarFactoryPtr tvf;

// continuity equation (div v = 0)
  VarPtr tauRho = vf ->testVar("tauRho", HGRAD);
  VarPtr rhoFlux = vf ->fluxVar("rhoFlux",L2); // currently v dot n

// momentum equation
  VarPtr v = vf->fieldVar("v", VECTOR_L2);
  VarPtr w = vf->fieldVar("w",L2); // pressure variable
  VarPtr vxFlux = vf ->fluxVar("vxFlux",L2);
  VarPtr vyFlux = vf ->fluxVar("vyFlux",L2);
  VarPtr tauVx = vf->testVar("tauVx", HGRAD);
  VarPtr tauVy = vf->testVar("tauVy", HGRAD);
  
// stress equation
  VarPtr sigmaxx = vf->fieldVar("sigmaxx", L2);
  VarPtr sigmaxy = vf->fieldVar("sigmaxy", L2);
  VarPtr sigmayy = vf->fieldVar("sigmayy", L2);
  VarPtr omega = vf->fieldVar("omega", L2);
  VarPtr vxHat = vf->traceVar("vxhat", HGRAD);
  VarPtr vyHat = vf->traceVar("vyhat", HGRAD);

  VarPtr tauSigmax = vf->testVar("tauSigmax", HDIV);
  VarPtr tauSigmay = vf->testVar("tauSigmay", HDIV);

  VarPtr chi; //chi is the normalized SA working varible
  VarPtr chiFlux;
  VarPtr tauChi;
  VarPtr psi; //psi = 1/Re 1/sigma_SA (1+chi_prev) grad_chi
  VarPtr chiHat;
  VarPtr tauPsi;


  if (useTurbModel)
    {
      if (splitTurbAdv)
	{
	  tvf = VarFactory::varFactory();
	  chi = tvf->fieldVar("chi", L2);
	  chiFlux = tvf->fluxVar("chiFlux", L2);
	  tauChi = tvf->testVar("tauChi", HGRAD);

	  psi = tvf->fieldVar("psi", VECTOR_L2);
	  chiHat = tvf->traceVar("chihat",HGRAD);
	  tauPsi = tvf->testVar("tauPsi", HDIV);
	}
      else
	{
	  chi = vf->fieldVar("chi", L2);
	  chiFlux = vf->fluxVar("chiFlux", L2);
	  tauChi = vf->testVar("tauChi", HGRAD);

	  psi = vf->fieldVar("psi", VECTOR_L2);
	  chiHat = vf->traceVar("chihat",HGRAD);
	  tauPsi = vf->testVar("tauPsi", HDIV);
	}
    }

/******************* Create the solver pointers ***************************/

/*
 Create pointers for the bilinear form, boundary conditions, rhs,
 inner product space, mesh, and the prevTime solution
*/

  BFPtr nonBF = BF::bf(vf);
  BCPtr bc = BC::bc();
  IPPtr ip = IP::ip(); //create an inner product space pointer
  RHSPtr rhs = RHS::rhs();

  BFPtr tBF; 
  BCPtr tbc; 
  IPPtr tip; 
  RHSPtr trhs;

  SolutionPtr tprevTime;
  MeshPtr mesh;
  MeshPtr tmesh;

  if (!readFile)
    {
      mesh = MeshFactory::minRuleMesh(meshTopo, nonBF, 
				      H1Order, delta_k);
    }
  else
    {
    mesh = MeshFactory::loadFromHDF5(nonBF, readFileName+".mesh");
    }
  
  
  SolutionPtr prevTime = Solution::solution(nonBF, mesh);
  mesh->registerSolution(prevTime);

  if (useTurbModel&&splitTurbAdv)
    {
      tBF = BF::bf(tvf);
      tbc = BC::bc();
      tip = IP::ip(); 
      trhs = RHS::rhs();

      if (!readFile)
	{
	  MeshTopologyPtr tMeshTopo = meshTopo->deepCopy();
	  tmesh = Teuchos::rcp( new Mesh(tMeshTopo, tBF, H1Order, delta_k) ) ;
	}
      else
	{
	  tmesh = MeshFactory::loadFromHDF5(tBF, readFileName+".mesh");
	}
      tprevTime = Solution::solution(tBF, tmesh);
      tmesh->registerSolution(tprevTime);
    }

 FunctionPtr sourceFactor = Teuchos::rcp(new RampFunction(sourceStep0,sourceStepf));
/***************** Define constant pointers ***********************************/
 
  FunctionPtr one = Function::constant(1.0);
  FunctionPtr zero = Function::constant(0);
  FunctionPtr invRe = Function::constant(1.0/Re);

/***************** Define Pointers for the prev time step *********************/

// Continuity
  FunctionPtr rhoFluxPrev = Function::solution(rhoFlux,prevTime,true); 

// momentum equation
  FunctionPtr vPrev =  Function::solution(v,prevTime);
  FunctionPtr wPrev =  Function::solution(w,prevTime);
  FunctionPtr vxFluxPrev = Function::solution(vxFlux,prevTime,true);
  FunctionPtr vyFluxPrev = Function::solution(vyFlux,prevTime,true);

// gradv equation
  FunctionPtr sigmaxxPrev = Function::solution(sigmaxx,prevTime);
  FunctionPtr sigmaxyPrev = Function::solution(sigmaxy,prevTime);
  FunctionPtr sigmayyPrev = Function::solution(sigmayy,prevTime);
  FunctionPtr omegaPrev  =  Function::solution(omega,prevTime);
  FunctionPtr vxHatPrev  =  Function::solution(vxHat,prevTime,false);
  FunctionPtr vyHatPrev  =  Function::solution(vyHat,prevTime,false);


  FunctionPtr chiPrev; 
  FunctionPtr chiFluxPrev;
  FunctionPtr tauChiPrev;
  FunctionPtr psiPrev; 
  FunctionPtr chiHatPrev;
  FunctionPtr tauPsiPrev;

  if(useTurbModel)
    {
      if (splitTurbAdv)
	{
	  chiPrev = Function::solution(chi,tprevTime); 
	  chiFluxPrev = Function::solution(chiFlux,tprevTime,true);
	  psiPrev = Function::solution(psi,tprevTime); 
	  chiHatPrev = Function::solution(chiHat,tprevTime,false);
	}
      else
	{
	  chiPrev = Function::solution(chi,prevTime); 
	  chiFluxPrev = Function::solution(chiFlux,prevTime,true);
	  psiPrev = Function::solution(psi,prevTime); 
	  chiHatPrev = Function::solution(chiHat,prevTime,false);
	}
    }
/*************** Define SA constants of function ptrs ********************/

  double const saCB1 = 0.1355;
  double const saCB2 = 0.622;
  double const saKappa = 0.41;
  double const saSigma = 2./3.;
  double const saCW1 = saCB1/(saKappa*saKappa)+(1.+saCB2)/(saSigma);
  double const saCW2 = 0.3;
  double const saCW3 = 2.0;
  double const saCV1 = 7.1;
  double const saCT3 = 1.2;
  double const saCT4 = 0.5;
  double const saRLim = 10.0;
  double const saCV2 = 0.7;
  double const saCV3 = 0.9;
  double const saCN1 =16.;
  double const smallNum =1.0e-20;


  FunctionPtr chiPrevSqr;
  FunctionPtr chiPrevCube;
  FunctionPtr chiPrevPos;
  FunctionPtr fv1Prev;
  FunctionPtr fv2Prev;
  FunctionPtr ft2Prev;

  FunctionPtr fnPrev;
  FunctionPtr sPrev;
  FunctionPtr sBarPrev;
  FunctionPtr sTildePrev;
  FunctionPtr sBarGtr;


  FunctionPtr saProd;
  FunctionPtr saDest;
  FunctionPtr saDiff;


  FunctionPtr dWall = Teuchos::rcp(new  distanceToPlate(0.0, 0.0));
  FunctionPtr dWallSqr =  Teuchos::rcp(new PowFunction( dWall, 2.0));
  FunctionPtr kappaDSqr = 
    Teuchos::rcp(new PowFunction(saKappa * dWall, 2.0));
  FunctionPtr smallNumPtr = Function::constant(smallNum);



  FunctionPtr rlim = Function::constant(saRLim);
  FunctionPtr saCV1Cube = Function::constant(pow(saCV1,3));
  FunctionPtr saCN1Ptr = Function::constant(saCN1);

  FunctionPtr prevVisc;





// define function pointers related to the destruction term
  FunctionPtr saR; 
  FunctionPtr saR6;
  FunctionPtr saG;
  FunctionPtr saCW3Six;
  FunctionPtr saG6;
  FunctionPtr fwFac;
  FunctionPtr fwPrev;


  if (useTurbModel)
    {
      chiPrevSqr = Teuchos::rcp(new PowFunction(chiPrev,2.));
      chiPrevCube = Teuchos::rcp(new PowFunction(chiPrev,3.));
      chiPrevPos = Teuchos::rcp(new GtrZeroFunction(chiPrev));
      fv1Prev = chiPrevCube/(chiPrevCube+saCV1Cube);
      fv2Prev = one - chiPrev/(one + chiPrev * fv1Prev);
 
      ft2Prev = Teuchos::rcp(new EFunction(-saCT4 * chiPrevSqr) );
      ft2Prev = saCT3 *  ft2Prev;

      sPrev = Teuchos::rcp(new AbsFunction(2.0 * Re * omegaPrev));   
      sBarPrev = invRe * chiPrev * fv2Prev / (kappaDSqr + smallNumPtr);
      sBarGtr = Teuchos::rcp(new GtrZeroFunction(sBarPrev + saCV2 * sPrev));

      if (useSBar)
	{
	  sTildePrev = sBarGtr * (sPrev + sBarPrev) +
	    (one-sBarGtr) * sPrev * (one + 
				     (saCV2 * saCV2 * sPrev + saCV3 * sBarPrev) / 
				     ( (saCV3 - 2 * saCV2) * sPrev - sBarPrev) );
	}
      else
	{
	  sTildePrev = sPrev + sBarPrev;
	}


      saR = Function::min(invRe * chiPrev / (sTildePrev * kappaDSqr + smallNumPtr), rlim);
      saR6 = Teuchos::rcp(new PowFunction(saR, 6) );
      saG = saR + saCW2 * (saR6 - saR);
      saCW3Six = Function::constant(pow(saCW3, 6.0) );
      saG6 = Teuchos::rcp(new PowFunction(saG, 6) );
      fwFac = Teuchos::rcp(new PowFunction( (one + saCW3Six) / (saG6 + saCW3Six), 1. / 6.0) );
      fwPrev = saG * fwFac;
     


      if (useSANeg)
	{
	  fnPrev = chiPrevPos + (one - chiPrevPos) * (saCN1Ptr + chiPrevCube) / ( saCN1Ptr - chiPrevCube );
	  saDiff = 1.0 / (saSigma * Re) * (one + fnPrev * chiPrev) ; 

	  saProd = saCB1 * chiPrev * (chiPrevPos *  ( one - ft2Prev ) * sTildePrev
	    + (1 - chiPrevPos) *  (one - saCT3) * sPrev );

	  saDest = (chiPrevPos * (saCW1 * fwPrev - saCB1 / pow(saKappa,2) * ft2Prev )
	    - (1 - chiPrevPos) *  saCW1) * invRe * chiPrevSqr / (dWallSqr + smallNumPtr) ; 

	  prevVisc = invRe * (one + chiPrevPos * fv1Prev * chiPrev);
	}
      else
	{
	  saDiff = 1.0 / (saSigma * Re) * (one + chiPrev);
	  saProd = saCB1 * ( one - ft2Prev ) * sTildePrev * chiPrev;
	  saDest = (saCW1 * fwPrev - saCB1 / pow(saKappa,2) * ft2Prev) * invRe * chiPrevSqr / (dWallSqr + smallNumPtr);
	  prevVisc = invRe * (one + fv1Prev * chiPrev);
	}
    }
  else
    {
      prevVisc = invRe;
    }


/************* Define function ptrs realted to dt and the mesh  ****************/
  FunctionPtr invDt = Function::constant(1.0/dt);
  FunctionPtr sqrtInvDt = Function::constant(sqrt(1.0/dt));
  FunctionPtr tinvDt = Function::constant(1.0*nTurbSteps/dt);
  FunctionPtr tsqrtInvDt = Function::constant(sqrt(1.0*nTurbSteps/dt));





/*************** Code up BF, RHS, and IP for each equation ********************/

  FunctionPtr e1 = Function::vectorize(one, zero);
  FunctionPtr e2 = Function::vectorize(zero, one);


// Continuity equation ( div v = 0 )
  nonBF->addTerm(rhoFlux, tauRho);
  nonBF->addTerm(-v, tauRho->grad());

  rhs->addTerm(-rhoFluxPrev * tauRho);
  rhs->addTerm(vPrev * tauRho->grad());

 
// grad Vx equation

  nonBF->addTerm(0.5/ prevVisc * (sigmaxx * e1 + sigmaxy * e2) , tauSigmax); 
  nonBF->addTerm(-Re * omega * e2, tauSigmax); 
  nonBF->addTerm(v->x(), tauSigmax->div());
  nonBF->addTerm(-vxHat, tauSigmax->dot_normal());

  rhs->addTerm(-0.5 *(e1 * sigmaxxPrev + e2 * sigmaxyPrev) / prevVisc * tauSigmax); 
  rhs->addTerm(Re * e2 * omegaPrev * tauSigmax); 
  rhs->addTerm(-vPrev->x() * tauSigmax->div());
  rhs->addTerm(vxHatPrev * tauSigmax->dot_normal());

// grad Vy equation
  nonBF->addTerm(0.5 / prevVisc * (sigmaxy * e1 + sigmayy * e2) , tauSigmay); 
  nonBF->addTerm(Re * omega * e1, tauSigmay); 
  nonBF->addTerm(v->y(), tauSigmay->div());
  nonBF->addTerm(-vyHat, tauSigmay->dot_normal());

  rhs->addTerm(-0.5 *(e1 * sigmaxyPrev + e2 * sigmayyPrev) / prevVisc * tauSigmay); 
  rhs->addTerm(-Re * e1 * omegaPrev * tauSigmay);  
  rhs->addTerm(-vPrev->y() * tauSigmay->div());
  rhs->addTerm(vyHatPrev * tauSigmay->dot_normal());


//momentum equation
// time derivative
  nonBF->addTerm(v->x() * invDt, tauVx);
  nonBF->addTerm(v->y() * invDt, tauVy);

// flux terms: n dot( vv+pI - nu gradV +...
  nonBF->addTerm(centerParam * vxFlux, tauVx);
  nonBF->addTerm(centerParam * vyFlux, tauVy);

  rhs->addTerm(-vxFluxPrev * tauVx);
  rhs->addTerm(-vyFluxPrev * tauVy);

// vv term
  nonBF->addTerm(-2.0 * centerParam * v->x() * vPrev->x(), tauVx->dx());

  nonBF->addTerm(-centerParam * vPrev->x() * v->y(), tauVy->dx());
  nonBF->addTerm(-centerParam * v->x() * vPrev->y(), tauVy->dx());

  nonBF->addTerm(-centerParam * vPrev->y() * v->x(), tauVx->dy());
  nonBF->addTerm(-centerParam * v->y() * vPrev->x(), tauVx->dy());

  nonBF->addTerm(-2.0 * centerParam * v->y() * vPrev->y(), tauVy->dy());

  rhs->addTerm(vPrev->x() * vPrev->x() * tauVx->dx());
  rhs->addTerm(vPrev->x() * vPrev->y() * tauVy->dx());
  rhs->addTerm(vPrev->y() * vPrev->x() * tauVx->dy());
  rhs->addTerm(vPrev->y() * vPrev->y() * tauVy->dy());

// grad p term
  nonBF->addTerm(-centerParam * w, tauVx->dx());
  nonBF->addTerm(-centerParam * w, tauVy->dy());

  rhs->addTerm(wPrev * tauVx->dx());
  rhs->addTerm(wPrev * tauVy->dy());

// viscosity terms

  nonBF->addTerm(centerParam * e1 * sigmaxx, tauVx->grad());
  nonBF->addTerm(centerParam * e2 * sigmaxy, tauVx->grad());
  nonBF->addTerm(centerParam * e1 * sigmaxy, tauVy->grad());
  nonBF->addTerm(centerParam * e2 * sigmayy, tauVy->grad());

  rhs->addTerm(-e1 * sigmaxxPrev * tauVx->grad());
  rhs->addTerm(-e2 * sigmaxyPrev * tauVx->grad());
  rhs->addTerm(-e1 * sigmaxyPrev * tauVy->grad());
  rhs->addTerm(-e2 * sigmayyPrev * tauVy->grad());


  if (useTurbModel)
    {
    if (splitTurbAdv)
            {
	//chi equation (chi is the normalized SA working variable)
	//dt term
	tBF->addTerm(chi * tinvDt, tauChi);

	// flux term: n dot V chi - psi
	tBF->addTerm(centerParam * chiFlux, tauChi);
	trhs->addTerm(-chiFluxPrev * tauChi);
	
	// advection term
	
	tBF->addTerm(-centerParam * vPrev * chi, tauChi->grad());
//      nonBF->addTerm(-centerParam * v * chiPrev, tauChi->grad());
	trhs->addTerm(vPrev * chiPrev * tauChi->grad());

	// add the diffusion term
	tBF->addTerm(centerParam * psi, tauChi->grad());
	trhs->addTerm(-psiPrev * tauChi->grad());

	// Production Term (Treated explicitly)
	trhs->addTerm(  sourceFactor * saProd * tauChi);

	// Destruction Term (Treated explicitly)
	trhs->addTerm( - sourceFactor * saDest * tauChi);

	// Grad Chi Squared Term (Treated explicitly) 
	// Grad Chi = psi/saDiff

	trhs->addTerm( saCB2/ (Re * saSigma) * psiPrev * psiPrev / (saDiff * saDiff) * tauChi);

	// gradchi eqn
	tBF->addTerm( psi / saDiff, tauPsi);
	tBF->addTerm(-chiHat, tauPsi->dot_normal());
	tBF->addTerm( chi, tauPsi->div());
	
	trhs->addTerm(-psiPrev / saDiff * tauPsi);
	trhs->addTerm( chiHatPrev * tauPsi->dot_normal());
	trhs->addTerm(-chiPrev * tauPsi->div());
      }
    else
      {
	//chi equation (chi is the normalized SA working variable)
	//dt term
	nonBF->addTerm(chi * invDt, tauChi);

	// flux term: n dot V chi - psi
	nonBF->addTerm(centerParam * chiFlux, tauChi);
	rhs->addTerm(-chiFluxPrev * tauChi);
	
	// advection term
	
	nonBF->addTerm(-centerParam * vPrev * chi, tauChi->grad());
//      nonBF->addTerm(-centerParam * v * chiPrev, tauChi->grad());
	rhs->addTerm(vPrev * chiPrev * tauChi->grad());

	// add the diffusion term
	nonBF->addTerm(centerParam * psi, tauChi->grad());
	rhs->addTerm(-psiPrev * tauChi->grad());

	// Production Term (Treated explicitly)
	rhs->addTerm(  sourceFactor * saProd * tauChi);

	// Destruction Term (Treated explicitly)
	rhs->addTerm( - sourceFactor * saDest * tauChi);

	// Grad Chi Squared Term (Treated explicitly) 
	// Grad Chi = psi/saDiff

	rhs->addTerm( saCB2/ (Re * saSigma) * psiPrev * psiPrev / (saDiff * saDiff) * tauChi);

	// gradchi eqn
	nonBF->addTerm( psi / saDiff, tauPsi);
	nonBF->addTerm(-chiHat, tauPsi->dot_normal());
	nonBF->addTerm( chi, tauPsi->div());
	
	rhs->addTerm(-psiPrev / saDiff * tauPsi);
	rhs->addTerm( chiHatPrev * tauPsi->dot_normal());
	rhs->addTerm(-chiPrev * tauPsi->div());
      }
    
    }

/******************* Define the inner product space ***************************/


  FunctionPtr invReSqrt =  Teuchos::rcp(new PowFunction( invRe, 0.5));

  FunctionPtr vecWork1 = Function::vectorize(one, -one);
  FunctionPtr vecWork2 = Function::vectorize(3.0*one, one);

  ipType desiredIp = IpStringEnum[ipSpace];
  switch ( desiredIp)
    {
    case ipType::simpleNorm:
     if (rank==0) cout <<"Using Simple Norm "  << endl;
     // This norm seems to work well for flows w/o SA model
 	ip->addTerm(invDt * tauVx);
	ip->addTerm(invDt * tauVy);
	ip->addTerm(tauRho);
	ip->addTerm(tauVx);
	ip->addTerm(tauVy);
	ip->addTerm(tauSigmax->div());
	ip->addTerm(tauSigmay->div());
	ip->addTerm(tauSigmax);
	ip->addTerm(tauSigmay);     
	if (useTurbModel)
	  {
	    if (splitTurbAdv)
	      {
		tip->addTerm(tinvDt * tauChi);
		tip->addTerm(tauChi);
		tip->addTerm(tauPsi);
		tip->addTerm(tauPsi->div());
		tip->addTerm(tauPsi->div() - vPrev * tauChi->grad());
		tip->addTerm(vPrev * tauChi->grad());	    
	      }
	    else
	      {
		ip->addTerm(invDt * tauChi);
		ip->addTerm(tauChi);
		ip->addTerm(tauPsi);
		ip->addTerm(tauPsi->div());
		ip->addTerm(tauPsi->div() - vPrev * tauChi->grad());
		ip->addTerm(vPrev * tauChi->grad());	    
	      }
	    }
      break;
    case ipType:: chanNorm:
     if (rank==0) cout <<"Using Jesse Chan's Norm "  << endl;
     if (rank==0) cout <<"Warning: this Norm is not working yet "  << endl;
     // Add time terms
	ip->addTerm(invDt * tauVx);
	ip->addTerm(invDt * tauVy);
	// Add ||V||^2 terms
	ip->addTerm(tauRho);
	ip->addTerm(tauVx);
	ip->addTerm(tauVy);
	// Add ||div W ||^2 Terms (not shown in paper but exist in NSDriver)
	ip->addTerm(tauSigmax->div());
	ip->addTerm(tauSigmay->div());
	// Add ||div W ||^2 Terms (not shown in paper but exist in NSDriver)
	// ip->addTerm(tauSigmax);
	//ip->addTerm(tauSigmay);

	// Add 1/Re ||A_visc^T Gard V||^2 Terms 
	ip->addTerm(invReSqrt * centerParam * tauVx->grad());
	ip->addTerm(invReSqrt * centerParam * tauVy->grad());
	// || A_euler^T Grad V||^2 Term
	ip->addTerm(tauRho->grad());
	ip->addTerm(centerParam * ((2*vPrev->x() + one) * e1 + (vPrev->x()+vPrev->y()) * e2) * tauVx->grad());
	ip->addTerm(centerParam * ((2*vPrev->y() + one) * e2 + (vPrev->x()+vPrev->y()) * e1) * tauVy->grad());
	// Add ||div W - A_euler^T grad V||^2 Terms
	ip->addTerm(tauSigmax->div() - centerParam * 
		   ((2*vPrev->x() + one) * e1 + (vPrev->x() + vPrev->y()) * e2) * tauVx->grad());
	ip->addTerm(tauSigmay->div() - centerParam *
		   ((2*vPrev->y() + one) * e2 + (vPrev->x() + vPrev->y()) * e1) * tauVy->grad());
	// Add min(Re, 1/K) ||E^T_visc W ||^2 Terms
	ip->addTerm(0.5 * vecWork1 / prevVisc * tauSigmax);
	ip->addTerm(0.5 * vecWork2 / prevVisc * tauSigmay);
      break;

    case ipType::graphNorm:
    default:
     if (rank==0) cout <<"Using Graph Norm "  << endl;

       ip = nonBF->graphNorm();
       if (useTurbModel&&splitTurbAdv)
	  {
	    tip = tBF->graphNorm();
	  }
    }



  
 

/******************* Apply the boundary conditions ***************************/
  cout << "Apply BC" << endl;

  SpatialFilterPtr top = SpatialFilter::matchingY(yDim);
  SpatialFilterPtr bottom = SpatialFilter::matchingY(0.0);
  SpatialFilterPtr inFlow = SpatialFilter::matchingX(-plateStart);
  SpatialFilterPtr outFlow = SpatialFilter::matchingX(plateLength);
  SpatialFilterPtr ltX = SpatialFilter::lessThanX(0.0);
  SpatialFilterPtr bottomGap = ltX & bottom;
  SpatialFilterPtr plate = !ltX & bottom;

  FunctionPtr flow = Function::constant(flowVelocity);
  FunctionPtr inflowVelocity = Function::vectorize(flow, zero);
  FunctionPtr nHat = Function::normal();


//  bc->addDirichlet(vxFlux, inFlow, e1 * flow * flow * nHat-vxFluxPrev);
  bc->addDirichlet(vxHat, inFlow, flow -vxHatPrev);
  bc->addDirichlet(vyFlux, inFlow, zero);
  bc->addDirichlet(vyHat, top, zero);
  bc->addDirichlet(vxFlux, top, zero);
  bc->addDirichlet(vxHat, plate, zero);
  bc->addDirichlet(vyHat, plate, zero);
  bc->addDirichlet(vyHat, bottomGap, zero);
  bc->addDirichlet(vxFlux, bottomGap, zero);
  bc->addZeroMeanConstraint(w);

  if (useTurbModel)
    {
      if (splitTurbAdv)
	{
	  tbc->addDirichlet(chiHat, inFlow, chiFsBc -chiHatPrev);
	  tbc->addDirichlet(chiFlux, top, zero);
	  tbc->addDirichlet(chiHat, plate, zero);
	  tbc->addDirichlet(chiFlux, bottomGap, zero);
	}
      else
	{
 //     bc->addDirichlet(chiFlux, inFlow, e1 * flow * chiFsBc * nHat-chiFluxPrev);
	  bc->addDirichlet(chiHat, inFlow, chiFsBc -chiHatPrev);
	  bc->addDirichlet(chiFlux, top, zero);
	  bc->addDirichlet(chiHat, plate, zero);
	  bc->addDirichlet(chiFlux, bottomGap, zero);
	}
    }

 

/******************* Define the soln are related pointers ***************************/

  cout << "Define Soln Pointers " << endl;
  SolutionPtr soln = Solution::solution(nonBF, mesh, bc, rhs, ip);
  SolutionPtr tsoln;

  cout << "Define Incr Pointers " << endl;
  FunctionPtr wIncr = Function::solution(w, soln);
  FunctionPtr vIncr = Function::solution(v, soln);
  FunctionPtr sigmaxxIncr = Function::solution(sigmaxx, soln);
  FunctionPtr sigmaxyIncr = Function::solution(sigmaxy, soln);
  FunctionPtr sigmayyIncr = Function::solution(sigmayy, soln);
  FunctionPtr omegaIncr = Function::solution(omega, soln);
  FunctionPtr l2WSqr = wIncr*wIncr;
  FunctionPtr l2VSqr = vIncr*vIncr;
  FunctionPtr l2SigmaxxSqr = sigmaxxIncr*sigmaxxIncr;
  FunctionPtr l2SigmaxySqr = sigmaxyIncr*sigmaxyIncr;
  FunctionPtr l2SigmayySqr = sigmayyIncr*sigmayyIncr;
  FunctionPtr l2OmegaSqr = omegaIncr*omegaIncr;


  FunctionPtr l2Sqr = l2WSqr + l2VSqr + l2SigmaxxSqr +
    l2SigmaxySqr +l2SigmayySqr + l2OmegaSqr;

  FunctionPtr l2PrevSolSqr =  wPrev * wPrev + vPrev * vPrev + sigmaxxPrev* sigmaxxPrev +
    sigmaxyPrev * sigmaxyPrev + sigmayyPrev * sigmayyPrev + omegaPrev * omegaPrev;

  FunctionPtr chiIncr;
  FunctionPtr psiIncr;
  FunctionPtr tl2Sqr;
  FunctionPtr l2ChiSqr;
  FunctionPtr l2PsiSqr;

  cout << "Use Turb Model If " << endl;
  if (useTurbModel)
    {
      if(splitTurbAdv)
	{
	  tsoln=Solution::solution(tBF, tmesh, tbc, trhs, tip);
	  chiIncr = Function::solution(chi, tsoln);
	  psiIncr = Function::solution(psi, tsoln);
	  l2ChiSqr = chiIncr*chiIncr;
	  l2PsiSqr = psiIncr*psiIncr;
	  tl2Sqr = l2ChiSqr + l2PsiSqr;
	}
      else
	{
	  chiIncr = Function::solution(chi, soln);
	  psiIncr = Function::solution(psi, soln);
	  l2ChiSqr = chiIncr*chiIncr;
	  l2PsiSqr = psiIncr*psiIncr;
	  l2Sqr = l2Sqr + l2ChiSqr + l2PsiSqr;
	}

      l2PrevSolSqr = l2PrevSolSqr + chiPrev * chiPrev + psiPrev * psiPrev;
    }

/******************* Add an initial guess ***************************/
  cout << "Add initial Guess" <<endl;
  const int solutionOrdinal = 0;
  if (!readFile)
    {
      map<int, FunctionPtr> initialGuess;
      map<int, FunctionPtr> tinitialGuess;

      FunctionPtr saYValue = Teuchos::rcp(new  yValue);
      
      FunctionPtr vxInitial = flowVelocity/yDim * saYValue;
      FunctionPtr gradVxyInitial = Function::constant(flowVelocity/yDim);
 
//  vxInitial = Function::constant(flowVelocity);
//  gradVxyInitial = zero;

      initialGuess[v->ID()] = Function::vectorize(vxInitial,zero);
      initialGuess[sigmaxx->ID()] = zero;
      initialGuess[sigmaxy->ID()] = gradVxyInitial * invRe;
      initialGuess[sigmayy->ID()] = zero;
      initialGuess[omega->ID()] = -0.5 * gradVxyInitial * invRe;
      initialGuess[w->ID()] = zero;
      
      FunctionPtr chiInitial; 
      FunctionPtr psiInitial;  

      if (useTurbModel)
	{
	  FunctionPtr normY = 1.0/yDim * saYValue;
	  FunctionPtr normY2 =  Teuchos::rcp(new PowFunction(normY , 2.0));
	  FunctionPtr normY3 =  Teuchos::rcp(new PowFunction(normY , 3.0));
	  FunctionPtr chi0 = chiFsBc * (-2.0 * normY3 + 3.0 * normY2);
	  FunctionPtr psi0y = 6.0*chiFsBc/(saSigma*yDim) * invRe * (one + chi0)*(normY-normY2); 
//      chiInitial = chiFsBc/yDim * saYValue;
//      psiInitial =  Function::vectorize(zero,zero); 
	  chiInitial = chi0;
	  psiInitial =  Function::vectorize(zero,psi0y); 
 //     chiInitial = Function::constant(chiFsBc);
 //     psiInitial =  Function::vectorize(zero,zero); 

	  if(splitTurbAdv)
	    {
	      tinitialGuess[chi->ID()] = chiInitial;
	      tinitialGuess[psi->ID()] = psiInitial;
	      tprevTime->projectOntoMesh(tinitialGuess, solutionOrdinal);
	    }
	  else
	    {
	      initialGuess[chi->ID()] = chiInitial;
	      initialGuess[psi->ID()] = psiInitial;
	    }
	}

      prevTime->projectOntoMesh(initialGuess, solutionOrdinal);
    }
  else
    {
      prevTime->loadFromHDF5(readFileName+".soln");
      if(useTurbModel&&splitTurbAdv)
	{
	  tprevTime->loadFromHDF5(readFileName+".tsoln");
	}
    }

  cout << "Here" << endl;
/******************* Solve the Navier Stokes Eqn ***************************/

  double l2Incr = 0;
  int refNumber = 0;
  int numRefinements = 2;
  HDF5Exporter exporter(mesh, "SA-FlatPlate");
  HDF5Exporter texporter(tmesh, "SA-Turb-FlatPlate");
  
//  HDF5Exporter fexporter(mesh, "SA-f");
  int numSubdivisions = 4;

  double newtonThreshold = 1.0e-3;
  double refThreshold = 1.0e-3;
  double energyThreshold = 0.2;
  RefinementStrategyPtr refStrategy = 
    RefinementStrategy::energyErrorRefinementStrategy(soln, energyThreshold);
  bool printToConsole = true;
  double energyError, l2Soln;
  
  soln->setUseCondensedSolve(true);
  exporter.exportSolution(prevTime, 0, numSubdivisions);
 
  if(useTurbModel&&splitTurbAdv)
    {
      tsoln->setUseCondensedSolve(true);
      texporter.exportSolution(tprevTime, 0, numSubdivisions);

    }
//  fexporter.exportFunction(saDest,"saDest", 0, numSubdivisions);

  int pltnum = 1;
  int step = 0;
  int tstep = 0;
  for (step = 0; step <= maxsteps; step++)
    {
      ((RampFunction*)sourceFactor.get())->UpdateStep(step);
      soln->solve();
      l2Incr = sqrt(l2Sqr->integrate(mesh));

      if (rank==0) cout <<"L^2(increment): " << l2Incr << endl;
//      if (!useTimeStep) 
//	{
//	  prevTime->addSolution(soln,newtonDamp);
//	}
//      else
//	{
	  prevTime->addSolution(soln,1.0);
//	}
	  if(useTurbModel&&splitTurbAdv)
	    {
	      for (tstep = 0; tstep <= nTurbSteps; tstep++)
		{
		  tsoln->solve();
		  l2Incr = l2Incr+sqrt(tl2Sqr->integrate(tmesh));
		  tprevTime->addSolution(tsoln,1.0);
		}
	    }
      if (step % ndum ==0)
      {
	if (rank==0) cout << "Writing solution to file at step: " << step  << endl;
	exporter.exportSolution(prevTime, pltnum, numSubdivisions);
//	fexporter.exportFunction(saDest,"saDest", pltnum, numSubdivisions);
	if(useTurbModel&&splitTurbAdv)
	  {
	    texporter.exportSolution(tprevTime, pltnum, numSubdivisions);
	  }
	pltnum++;
      }
      if (false)
//      if (step % nsave ==0)

	{
	  if (rank==0)
	    {
	      cout<<"Saving solution to file at step: "<< step <<endl;
	    }
	  string oss;
	  oss="";
	  oss = filename+std::to_string(step);
	  
	  prevTime->mesh()->saveToHDF5(oss+".mesh");
	  prevTime->saveToHDF5(oss+".soln");
	  if(useTurbModel&&splitTurbAdv)
	    {
	      tprevTime->mesh()->saveToHDF5(oss+".tmesh");
	      tprevTime->saveToHDF5(oss+".tsoln");
	    }
	}
    }

#endif
  return 0;

}


/******************************* End of Code **********************************/



/* skip refinement until done debugging
  for (refNumber =0; refNumber <= numRefinements; refNumber++)
    {
      do
	{
	  soln->solve();
	  l2_incr = sqrt(l2_squared->integrate(mesh));
	  if (rank==0) cout <<"L^2(increment): " << l2_incr << endl;
	  backFlow->addSolution(soln,1.0);
	}
      while (l2_incr > newtonThreshold);

      exporter.exportSolution(backFlow, refNumber, numSubdivisions);
      
      refStrategy->refine(printToConsole);
      energyError = refStrategy->getEnergyError(refNumber);
      l2_soln = sqrt(l2_backFlow_squared->integrate(mesh));
      
      newtonThreshold = refThreshold * energyError/l2_soln;
      cout << "L2(soln): " << l2_soln << endl;
      cout << "Newton threshold: " << newtonThreshold << endl;
    }
  energyError = soln->energyErrorTotal();
  cout << "Final energy error: " << energyError << endl;

  return 0;
}
*/
