/*
  $Id: SemiImpModelUW 26 2018-06-08 19:30:25Z ehowell $

  This code modeling the propogation of sound waves with a background flow.
  The purpose is to test the behavior of the semi-implicit operator with
  large flows, and to see if the DPG methods allows us to take large time
  steps.
 */

#include "Camellia.h" 
#include <math.h>

using namespace Camellia;
using namespace std;

/********************* Helper Functions ************************/


class coswtkx : public SimpleFunction<double>
{
// return the cos(omega t - kx) omega t is given as an input
  double _k;
  double _wt;//
public:
  coswtkx(double k, double wt)
    {
      _k=k;
      _wt=wt;
    }
  double value(double x, double y)
  {
    return cos(_wt - _k * x);
  }
};

class sinwtkx : public SimpleFunction<double>
{
// return the cos(omega t - kx) omega t is given as an input
  double _k;
  double _wt;//
public:
  sinwtkx(double k, double wt)
    {
      _k=k;
      _wt=wt;
    }
  double value(double x, double y)
  {
    return sin(_wt - _k * x);
  }
};

/********************* Main driver program ************************/

int main(int argc, char *argv[]) 
{

  Teuchos::GlobalMPISession mpiSession(&argc, &argv); // initalize mpi
  int rank = MPIWrapper::CommWorld()->MyPID();

  Teuchos::CommandLineProcessor CLP;


/******************** Mesh related inputs *********************/
  int mx = 8; 
  int my = 8;

//  double xDim = M_PI;
//  double yDim = M_PI;
  double xDim = 1.;
  double yDim = 1.;

  int H1Order = 2;
  int delta_k = 2;


/******************** Time step inputs *********************/
  double dt = 100;
  double center = 0.5;
  double semiCo = 0.51;
  int nstep = 1000;
  int ndump = 10;


  bool addAdvectIp = true;
  bool addFluxes = false;
/******************** Physics inputs *********************/

  bool nonlinear = false;
  double Re = 1.0e8;
  double Mach = 1.;
  double gamma = 5.0/3.0;

/******************** Equilibrium inputs *********************/
  double v0 = 1; //intial flow veloicty 
  double theta = M_PI/4.0; // direction of v0 (0 -> x, pi/2 ->y)
    
  double amp = 1.0e-4; // amplitude of pressure perturbation
  int nx = 1; // number of periods of sound waves in x direction


/******************* Set up Command Line Processor  ***************************/
  CLP.setOption("mx", &mx, "Number of Elements in x direction");
  CLP.setOption("my", &my, "Number of Elements in y direction");
  CLP.setOption("H1Order", &H1Order, "Polynomial Degree");
  CLP.setOption("delta_k", &delta_k, "Polynomial Degree Enhancement");
  CLP.setOption("xDim", &xDim, "Domain Height");
  CLP.setOption("yDim", &yDim, "Domain Height");


  CLP.setOption("dt", &dt, "Time Step");
  CLP.setOption("center", &center, "centering parameter");
  CLP.setOption("semiCo", &semiCo, "semiimplicit operator coeff");

  CLP.setOption("Re", &Re, "Reynolds number");
  CLP.setOption("Mach", &Mach, "Mach Number");
  CLP.setOption("gamma", &gamma, "Ratio of specific heats");

  CLP.parse( argc, argv );


/******************* Begin Code ***************************/
/******************* Define variables for calculation *************************/
  VarFactoryPtr vVf  = VarFactory::varFactory();
  VarFactoryPtr pVf  = VarFactory::varFactory();
  

  VarPtr p = pVf->fieldVar("p", L2);
  VarPtr pTau = pVf->testVar("pTau", HGRAD);
  VarPtr pFlux = pVf->fluxVar("pFlux", L2);

  VarPtr gradP = pVf->fieldVar("gradP", VECTOR_L2);
  VarPtr gradPTau = pVf->testVar("gradPTau", HDIV);
  VarPtr pHat = pVf->traceVar("pHat", L2);
  
  VarPtr vx = vVf->fieldVar("vx", L2);
  VarPtr vy = vVf->fieldVar("vy", L2);
  VarPtr vxFlux = vVf->fluxVar("vxFLux", L2);
  VarPtr vyFlux = vVf->fluxVar("vyFLux", L2);
  VarPtr vxTau = vVf->testVar("vxTau", HGRAD);
  VarPtr vyTau = vVf->testVar("vyTau", HGRAD);


  VarPtr gradVx = vVf->fieldVar("gradVx", VECTOR_L2);
  VarPtr gradVy = vVf->fieldVar("gradVy", VECTOR_L2);
  VarPtr gradVxTau = vVf->testVar("gradVxTau", HDIV);
  VarPtr gradVyTau = vVf->testVar("gradVyTau", HDIV);
  VarPtr vxHat = vVf->traceVar("vxHat", L2);
  VarPtr vyHat = vVf->traceVar("vyHat", L2);


/******************* Create the solver pointers ***************************/

  BFPtr pBF = BF::bf(pVf);
  BCPtr pBc = BC::bc();
  IPPtr pIp = IP::ip(); 
  RHSPtr pRHS = RHS::rhs();

  BFPtr vBF = BF::bf(vVf);
  BCPtr vBc = BC::bc();
  IPPtr vIp = IP::ip(); 
  RHSPtr vRHS = RHS::rhs();


/******************* Set up periodic BCs ***************************/
  std:vector<PeriodicBCPtr> periodicBCs;
  auto pBCx = PeriodicBC::xIdentification(0,xDim);
  auto pBCy = PeriodicBC::yIdentification(0,yDim);
  periodicBCs.push_back(pBCx);  
  periodicBCs.push_back(pBCy);  



/******************* Set up Mesh ***************************/

  vector<double> dims = {xDim, yDim}; // dimensions of the mesh
  vector<int> meshDims = {mx,my}; 
  vector<double> x0 = {0.0, 0.0}; 
  MeshTopologyPtr meshTopo =
    MeshFactory::rectilinearMeshTopology(dims, meshDims, x0, periodicBCs);

  MeshPtr pMesh = MeshFactory::minRuleMesh(meshTopo, pBF, H1Order, delta_k);
  MeshPtr vMesh = MeshFactory::minRuleMesh(meshTopo, vBF, H1Order, delta_k);


  SolutionPtr pPrevTime = Solution::solution(pBF, pMesh);
  pMesh->registerSolution(pPrevTime);
  
  SolutionPtr vPrevTime = Solution::solution(vBF, vMesh);
  vMesh->registerSolution(vPrevTime);

/**************** Define useful function pointer *********************/
  FunctionPtr zero = Function::constant(0.0); 

  FunctionPtr invDt = Function::constant(1.0/dt); 
  FunctionPtr invRe = Function::constant(1.0/Re); 
  FunctionPtr invGammaM2 = Function::constant(1.0/(gamma*Mach*Mach) ); 

  FunctionPtr centerPtr = Function::constant(center); 
  FunctionPtr gammaPtr = Function::constant(gamma); 
  FunctionPtr semiPtr = Function::constant(dt * semiCo / (Mach * Mach));

// set up equilibrium
  FunctionPtr vxEq = cos(theta) * Function::constant(v0);
  FunctionPtr vyEq = sin(theta) * Function::constant(v0);
  FunctionPtr pEq = Function::constant(1.0); //presure has been normalize

  FunctionPtr vxPrev = Function::solution(vx,vPrevTime);
  FunctionPtr vyPrev = Function::solution(vy,vPrevTime);

  FunctionPtr vxFluxPrev = Function::solution(vxFlux,vPrevTime,true);
  FunctionPtr vyFluxPrev = Function::solution(vyFlux,vPrevTime,true);

  FunctionPtr pPrev  = Function::solution(p, pPrevTime);
  FunctionPtr gradPPrev  = Function::solution(gradP, pPrevTime);
  FunctionPtr pHatPrev  = Function::solution(pHat, pPrevTime,false);

  FunctionPtr gradVxPrev  = Function::solution(gradVx, vPrevTime);
  FunctionPtr vxHatPrev  = Function::solution(vxHat, vPrevTime,false);
  FunctionPtr gradVyPrev  = Function::solution(gradVy, vPrevTime);
  FunctionPtr vyHatPrev  = Function::solution(vyHat, vPrevTime,false);


/**************** Define BF, ip, and RHS *********************/


  pBF->addTerm(invDt * p,pTau);

  pBF->addTerm(centerPtr * vxEq * gradP->x(), pTau);
  pBF->addTerm(centerPtr * vyEq * gradP->y(), pTau);

// assume eq flow is incompressible
/*
  pBf->addTerm(centerPtr * gammaPtr * vxEq->dx() * p, pTau);
  pBf->addTerm(centerPtr * gammaPtr * vyEq->dy() * p, pTau);
*/

  pRHS->addTerm(-vxEq * gradPPrev->x() * pTau);
  pRHS->addTerm(-vyEq * gradPPrev->y() * pTau);

  pRHS->addTerm(-vxPrev * pEq->dx() * pTau);
  pRHS->addTerm(-vyPrev * pEq->dy() * pTau);

  pRHS->addTerm(-gammaPtr * gradVxPrev->x() * pEq * pTau);
  pRHS->addTerm(-gammaPtr * gradVyPrev->y() * pEq * pTau);
// assume Eq flow is incompressible
/*
  pRHS->addTerm(-gammaPtr * vxEq->dx() * pPrev * pTau);
  pRHS->addTerm(-gammaPtr * vyEq->dy() * pPrev * pTau);
*/

  // add semi implicit operator to pressure equation
  pBF->addTerm(semiPtr * gradP, pTau->grad());
  pBF->addTerm(-semiPtr * pFlux, pTau);

  if (nonlinear)
    {
      pBF->addTerm(centerPtr * vxPrev * gradP->x(), pTau);
      pBF->addTerm(centerPtr * vyPrev * gradP->y(), pTau);

      pBF->addTerm(centerPtr * gammaPtr * gradVxPrev->x() * p, pTau);
      pBF->addTerm(centerPtr * gammaPtr * gradVyPrev->y() * p, pTau);

      pRHS->addTerm(-vxPrev * gradPPrev->x() * pTau);
      pRHS->addTerm(-vyPrev * gradPPrev->y() * pTau);
      pRHS->addTerm(-gammaPtr * gradVxPrev->x() * pPrev * pTau);
      pRHS->addTerm(-gammaPtr * gradVyPrev->y() * pPrev * pTau);
    }

// grad P equation
  pBF->addTerm( centerPtr * gradP, gradPTau);
  pBF->addTerm(-centerPtr * pHat, gradPTau->dot_normal());
  pBF->addTerm( centerPtr * p, gradPTau->div());

  pRHS->addTerm(-gradPPrev * gradPTau);
  pRHS->addTerm( pHatPrev * gradPTau->dot_normal());
  pRHS->addTerm(-pPrev *  gradPTau->div());
 


// velocity equation
  vBF->addTerm(invDt * vx, vxTau);
  vBF->addTerm(invDt * vy, vyTau);

// pressure term
  vRHS->addTerm(-invGammaM2 * gradPPrev->x() * vxTau);
  vRHS->addTerm(-invGammaM2 * gradPPrev->y() * vyTau);


// Equilibrium pressure term should cancel other eq terms
//  vRHS->addTerm(-invGammaM2 * pEq->dx() * vxTau);
//  vRHS->addTerm(-invGammaM2 * pEq->dy() * vyTau);

// viscoisty term


  vBF->addTerm(invRe * centerPtr * gradVx->x(), vxTau->dx());
  vBF->addTerm(invRe * centerPtr * gradVx->y(), vxTau->dy());
  vBF->addTerm(invRe * centerPtr * gradVy->x(), vyTau->dx());
  vBF->addTerm(invRe * centerPtr * gradVy->y(), vyTau->dy());

  vRHS->addTerm(-invRe * gradVxPrev->x() * vxTau->dx());
  vRHS->addTerm(-invRe * gradVxPrev->y() * vxTau->dy());
  vRHS->addTerm(-invRe * gradVyPrev->x() * vyTau->dx());
  vRHS->addTerm(-invRe * gradVyPrev->y() * vyTau->dy());

// no equilibrium viscoisty to due equilbrium assumption

// add flux term 
  vBF->addTerm(-centerPtr * vxFlux, vxTau);
  vBF->addTerm(-centerPtr * vyFlux, vyTau);

  vRHS->addTerm(vxFluxPrev * vxTau);
  vRHS->addTerm(vyFluxPrev * vyTau);


// advection terms  
  vBF->addTerm(centerPtr * vxEq * gradVx->x(), vxTau);
  vBF->addTerm(centerPtr * vyEq * gradVx->y(), vxTau);
  vBF->addTerm(centerPtr * vxEq * gradVy->x(), vyTau);
  vBF->addTerm(centerPtr * vyEq * gradVy->y(), vyTau);

  vBF->addTerm(centerPtr * vx * vxEq->dx(), vxTau);
  vBF->addTerm(centerPtr * vy * vxEq->dy(), vxTau);
  vBF->addTerm(centerPtr * vx * vyEq->dx(), vyTau);
  vBF->addTerm(centerPtr * vy * vyEq->dy(), vyTau);

  vRHS->addTerm(-vxEq * gradVxPrev->x() * vxTau);
  vRHS->addTerm(-vyEq * gradVxPrev->y() * vxTau);
  vRHS->addTerm(-vxEq * gradVyPrev->x() * vyTau);
  vRHS->addTerm(-vyEq * gradVyPrev->y() * vyTau);

  vRHS->addTerm(-vxPrev * vxEq->dx() * vxTau);
  vRHS->addTerm(-vyPrev * vxEq->dy() * vxTau);
  vRHS->addTerm(-vxPrev * vyEq->dx() * vyTau);
  vRHS->addTerm(-vyPrev * vyEq->dy() * vyTau);


  if (nonlinear)
    {
      vBF->addTerm(centerPtr * vxPrev * gradVx->x(), vxTau);
      vBF->addTerm(centerPtr * vyPrev * gradVx->y(), vxTau);
      vBF->addTerm(centerPtr * vxPrev * gradVy->x(), vyTau);
      vBF->addTerm(centerPtr * vyPrev * gradVy->y(), vyTau);

      vBF->addTerm(centerPtr * vx * gradVxPrev->x(), vxTau);
      vBF->addTerm(centerPtr * vy * gradVxPrev->y(), vxTau);
      vBF->addTerm(centerPtr * vx * gradVyPrev->x(), vyTau);
      vBF->addTerm(centerPtr * vy * gradVyPrev->y(), vyTau);

      vRHS->addTerm(-vxPrev * gradVxPrev->x() * vxTau);
      vRHS->addTerm(-vyPrev * gradVxPrev->y() * vxTau);
      vRHS->addTerm(-vxPrev * gradVyPrev->x() * vyTau);
      vRHS->addTerm(-vyPrev * gradVyPrev->y() * vyTau);

      // assume equilibrium so neglect veq cdot nabla veq
  
    }

// grad vx equation
  vBF->addTerm( centerPtr * gradVx, gradVxTau);
  vBF->addTerm(-centerPtr * vxHat, gradVxTau->dot_normal());
  vBF->addTerm( centerPtr * vx, gradVxTau->div());

  vRHS->addTerm(-gradVxPrev * gradVxTau);
  vRHS->addTerm( vxHatPrev * gradVxTau->dot_normal());
  vRHS->addTerm(-vxPrev *  gradVxTau->div());

// grad vy equation
  vBF->addTerm( centerPtr * gradVy, gradVyTau);
  vBF->addTerm(-centerPtr * vyHat, gradVyTau->dot_normal());
  vBF->addTerm( centerPtr * vy, gradVyTau->div());

  vRHS->addTerm(-gradVyPrev * gradVyTau);
  vRHS->addTerm( vyHatPrev * gradVyTau->dot_normal());
  vRHS->addTerm(-vyPrev *  gradVyTau->div());


/******************* Define the inner product ***************************/

/*
  vIp = vBF->graphNorm();
  pIp = pBF->graphNorm();
*/

  vIp->addTerm(invDt * vxTau);
  vIp->addTerm(invDt * vyTau);
  vIp->addTerm(vxTau->grad());
  vIp->addTerm(vyTau->grad());
  pIp->addTerm(invDt * pTau);
  pIp->addTerm(pTau->grad());

  if (addAdvectIp)
    {
      pIp->addTerm(vxEq * pTau->dx() + vyEq * pTau->dy());
      vIp->addTerm(vxEq * vxTau->dx() + vyEq * vxTau->dy());
      vIp->addTerm(vxEq * vyTau->dx() + vyEq * vyTau->dy());
    }

  pIp->addTerm(gradPTau);
  pIp->addTerm(gradPTau->div());
  vIp->addTerm(gradVxTau);
  vIp->addTerm(gradVyTau);
  vIp->addTerm(gradVxTau->div());
  vIp->addTerm(gradVyTau->div());

/******************* Initialize perturbed fields ***************************/
  int solutionOrdinal = 0;

  map<int, FunctionPtr> vPert;
  map<int, FunctionPtr> pPert;

  double kxPert = 2 * M_PI/(xDim * nx) ;

  FunctionPtr vxPertFunc = Teuchos::rcp(new coswtkx(kxPert, 0));
  FunctionPtr gradVxPertFunc = Teuchos::rcp(new sinwtkx(kxPert, 0));

  vPert[vx->ID()] = amp /(gamma * Mach) * vxPertFunc ;
  vPert[vy->ID()] = zero;

  vPert[gradVx->ID()] = 
    Function::vectorize(-kxPert * amp /(gamma * Mach) * gradVxPertFunc,
			zero);

  vPert[gradVy->ID()] = Function::vectorize(zero,zero);

  double wdt = Mach * dt / (2 * kxPert); // pressire is at 1/2 time step
  FunctionPtr pPertFunc = Teuchos::rcp(new coswtkx(kxPert, wdt));
  FunctionPtr gradPPertFunc = Teuchos::rcp(new sinwtkx(kxPert, wdt));

  pPert[p->ID()] = amp * pPertFunc;
  pPert[gradP->ID()] = Function::vectorize(-kxPert * amp * gradPPertFunc,zero);

  vPrevTime->projectOntoMesh(vPert, solutionOrdinal);
  pPrevTime->projectOntoMesh(pPert, solutionOrdinal);


/******************* Apply the boundary conditions ***************************/
/*
  I have a doubly periodic box, so I don't have to do anything here

  SpatialFilterPtr top = SpatialFilter::matchingY(yDim);
  SpatialFilterPtr bottom = SpatialFilter::matchingY(0.0);
  SpatialFilterPtr left = SpatialFilter::matchingX(0.0);
  SpatialFilterPtr right = SpatialFilter::matchingX(xDim);


  FunctionPtr zero = Function::constant(0.0);



  cout << "before add dirichlet" << endl;
  bc->addDirichlet(phi, top, one);
  bc->addDirichlet(phi, bottom, one);
  bc->addDirichlet(phi, left, one);
  bc->addDirichlet(phi, right, one);
  cout << "After add dirichlet" << endl;
*/


/******************* Define the soln  pointers ***************************/


  SolutionPtr pSoln = Solution::solution(pBF, pMesh, pBc, pRHS, pIp);
  SolutionPtr vSoln = Solution::solution(vBF, vMesh, vBc, vRHS, vIp);
  
  FunctionPtr vxIncr = Function::solution(vx, vSoln);
  FunctionPtr vyIncr = Function::solution(vy, vSoln);
  FunctionPtr gradVxIncr = Function::solution(gradVx, vSoln);
  FunctionPtr gradVyIncr = Function::solution(gradVy, vSoln);

  FunctionPtr pIncr = Function::solution(p, pSoln);
  FunctionPtr gradPIncr = Function::solution(gradP, pSoln);


  FunctionPtr vL2Sqr = vxIncr * vxIncr + vyIncr * vyIncr 
                 + gradVxIncr * gradVxIncr + gradVyIncr * gradVyIncr;
  FunctionPtr pL2Sqr = pIncr * pIncr + gradPIncr * gradPIncr;

/******************* Solve  ***************************/

  double vL2Incr = 0;
  double pL2Incr = 0;
  HDF5Exporter vExporter(vMesh, "SoundWavesFlow");
  HDF5Exporter pExporter(pMesh, "SoundWavesPressure");
  int numSubdivisions = 4;
  
//  pSoln->setUseCondensedSolve(true);
//  vSoln->setUseCondensedSolve(true);
  vExporter.exportSolution(vPrevTime, 0, numSubdivisions);
  pExporter.exportSolution(pPrevTime, 0, numSubdivisions);

  int pltnum = 1;
  int step = 0;

  for (step = 0; step <= nstep; step++)
    {
      vSoln->solve();
      vL2Incr = sqrt(vL2Sqr->integrate(vMesh));
      
      if (rank==0) cout <<"V L^2(increment): " << vL2Incr << endl;
      vPrevTime->addSolution(vSoln,1.0);

      pSoln->solve();
      pL2Incr = sqrt(pL2Sqr->integrate(pMesh));
      
      if (rank==0) cout <<"P L^2(increment): " << pL2Incr << endl;
      pPrevTime->addSolution(pSoln,1.0);


      if (step % ndump ==0)
	{
	  if (rank==0) cout << "Writing dumpfile at step: " << step << endl;
	  vExporter.exportSolution(vPrevTime, pltnum, numSubdivisions);
	  pExporter.exportSolution(pPrevTime, pltnum, numSubdivisions);
	  pltnum++;
	}
    }


  return 0;

}



