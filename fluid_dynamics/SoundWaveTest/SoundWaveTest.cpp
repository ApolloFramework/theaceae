/*
  $Id: SoundWaveTest ehowell$

  This code models the propogation of sound waves with a background flow.
  The purpose is to test the behavior of the semi-implicit operator with
  large flows, and to see if the DPG methods allows us to take large time
  steps than other methods. Thre diffrenet methods will be implamented
  cg, primal dpg, and ultraweak dpg.
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
// return the sin(omega t - kx) omega t is given as an input
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
  double xDim = 1.;
  double yDim = 1.;
  int H1Order = 2;
  int delta_k = 2; 

/******************** Finite element type inputs *********************/
  enum class feType
  {
    CG,
    pDPG,
    uwDPG
  };

  std::map<std::string, feType> feTypeEnum = 
    {
      {"continuous", feType::CG},
      {"primal", feType::pDPG},
      {"ultraweak", feType::uwDPG},
    };

  string feMethod = "ultraweak";

/******************** Time step inputs *********************/
  double dt = 100;
  double center = 0.5;
  double semiCo = 0.25;
  int nstep = 1000;
  int ndump = 100;
  bool addAdvectIp = true;
  bool addFluxes = false;

/******************** Physics inputs *********************/
  bool nonlinear = false;
  double Re = 1.0e8;
  double Mach = 1.;
  double gamma = 5.0/3.0;

/******************** Equilibrium and perturbation inputs *********************/
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
  CLP.setOption("nstep", &nstep, "number of time steps");
  CLP.setOption("ndump", &ndump, "frequency to write plot data");
  CLP.setOption("Re", &Re, "Reynolds number");
  CLP.setOption("Mach", &Mach, "Mach Number");
  CLP.setOption("gamma", &gamma, "Ratio of specific heats");
  CLP.setOption("feMethod", &feMethod, "Finite Element Method");

  CLP.parse( argc, argv );

/******************* Define variables for calculation *************************/
  feType thisFeMethod = feTypeEnum[feMethod];
  VarFactoryPtr vVf  = VarFactory::varFactory();
  VarFactoryPtr pVf  = VarFactory::varFactory();  
  VarPtr p, pTau, pFlux, gradP, gradPTau, pHat;
  VarPtr vx, vy, vxFlux, vyFlux, vxTau, vyTau;
  VarPtr gradVx, gradVy, gradVxTau, gradVyTau, vxHat, vyHat;

  switch (thisFeMethod)
    {
      case feType::CG:
	p = pVf->fieldVar("p", HGRAD);
	pTau = pVf->testVar("pTau", HGRAD);
	vx = vVf->fieldVar("vx", HGRAD);
	vy = vVf->fieldVar("vy", HGRAD);
	vxTau = vVf->testVar("vxTau", HGRAD);
	vyTau = vVf->testVar("vyTau", HGRAD);
      break;
      case feType::pDPG:
// todo
      break;
      case feType::uwDPG:
	p = pVf->fieldVar("p", L2);
	pTau = pVf->testVar("pTau", HGRAD);
	pFlux = pVf->fluxVar("pFlux", L2);
	gradP = pVf->fieldVar("gradP", VECTOR_L2);
	gradPTau = pVf->testVar("gradPTau", HDIV);
	pHat = pVf->traceVar("pHat", L2);
	vx = vVf->fieldVar("vx", L2);
	vy = vVf->fieldVar("vy", L2);
	vxFlux = vVf->fluxVar("vxFLux", L2);
	vyFlux = vVf->fluxVar("vyFLux", L2);
	vxTau = vVf->testVar("vxTau", HGRAD);
	vyTau = vVf->testVar("vyTau", HGRAD);
	gradVx = vVf->fieldVar("gradVx", VECTOR_L2);
	gradVy = vVf->fieldVar("gradVy", VECTOR_L2);
	gradVxTau = vVf->testVar("gradVxTau", HDIV);
	gradVyTau = vVf->testVar("gradVyTau", HDIV);
	vxHat = vVf->traceVar("vxHat", L2);
	vyHat = vVf->traceVar("vyHat", L2);
      break;
    };

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
  vector<double> dims = {xDim, yDim}; 
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
  TFunctionPtr<double> n = TFunction<double>::normal(); //todo
  FunctionPtr boundaryIndicator = Function::meshBoundaryCharacteristic();
  FunctionPtr invDt = Function::constant(1.0/dt); 
  FunctionPtr invRe = Function::constant(1.0/Re); 
  FunctionPtr invGammaM2 = Function::constant(1.0/(gamma*Mach*Mach) ); 
  FunctionPtr centerPtr = Function::constant(center); 
  FunctionPtr gammaPtr = Function::constant(gamma); 
  FunctionPtr semiPtr = Function::constant(dt * semiCo / (Mach * Mach));

/**************** Set Up Equilibrium *********************/
  FunctionPtr vxEq = cos(theta) * Function::constant(v0);
  FunctionPtr vyEq = sin(theta) * Function::constant(v0);
  FunctionPtr pEq = Function::constant(1.0); //presure has been normalized

/**************** Define pointers to prev soln *********************/
  FunctionPtr pPrev  = Function::solution(p, pPrevTime);
  FunctionPtr vxPrev = Function::solution(vx,vPrevTime);
  FunctionPtr vyPrev = Function::solution(vy,vPrevTime);
  FunctionPtr gradPPrev, pHatPrev, vxFluxPrev, vyFluxPrev;
  FunctionPtr gradVxPrev, vxHatPrev, gradVyPrev, vyHatPrev; 

  switch (thisFeMethod)
    {
      case feType::CG:
// do nothing
      break;
      case feType::pDPG:
	vxFluxPrev = Function::solution(vxFlux,vPrevTime,true);
	vyFluxPrev = Function::solution(vyFlux,vPrevTime,true);
      break;
      case feType::uwDPG:
	gradPPrev  = Function::solution(gradP, pPrevTime);
	pHatPrev   = Function::solution(pHat, pPrevTime,false);
	vxFluxPrev = Function::solution(vxFlux,vPrevTime,true);
	vyFluxPrev = Function::solution(vyFlux,vPrevTime,true);
	gradVxPrev = Function::solution(gradVx, vPrevTime);
	gradVyPrev = Function::solution(gradVy, vPrevTime);
	vxHatPrev  = Function::solution(vxHat, vPrevTime,false);
	vyHatPrev  = Function::solution(vyHat, vPrevTime,false);
      break;
    };
/**************** Define BF, ip, and RHS *********************/
  switch (thisFeMethod)
    {
      case feType::CG:
// todo check me
	// pressure equation
	pBF->addTerm(invDt * p,pTau);
	pBF->addTerm(centerPtr * vxEq * p->dx(), pTau);
	pBF->addTerm(centerPtr * vyEq * p->dy(), pTau);
	pRHS->addTerm(-vxEq * pPrev->dx() * pTau);
	pRHS->addTerm(-vyEq * pPrev->dy() * pTau);
	pRHS->addTerm(-vxPrev * pEq->dx() * pTau);
	pRHS->addTerm(-vyPrev * pEq->dy() * pTau);
	pRHS->addTerm(-gammaPtr * vxPrev->dx() * pEq * pTau);
	pRHS->addTerm(-gammaPtr * vyPrev->dy() * pEq * pTau);
//	pBF->addTerm(semiPtr * p->grad(), pTau->grad());
	// add boundary flux term
//	pBF->addTerm(-semiPtr * p->grad() * n, boundaryIndicator * pTau);
// ech I think there's a bug related to calcualting the dimensionality of grad()
	pBF->addTerm(semiPtr * p->dx(), pTau->dx());
	pBF->addTerm(semiPtr * p->dy(), pTau->dy());
	// add boundary flux term
	pBF->addTerm(-semiPtr * p->dx() * n->x(), boundaryIndicator * pTau);
	pBF->addTerm(-semiPtr * p->dy() * n->y(), boundaryIndicator * pTau);
	if (nonlinear)
	  {
	    pBF->addTerm(centerPtr * vxPrev * p->dx(), pTau);
	    pBF->addTerm(centerPtr * vyPrev * p->dy(), pTau);
	    pBF->addTerm(centerPtr * gammaPtr * vxPrev->dx() * p, pTau);
	    pBF->addTerm(centerPtr * gammaPtr * vyPrev->dy() * p, pTau);
	    pRHS->addTerm(-vxPrev * pPrev->dx() * pTau);
	    pRHS->addTerm(-vyPrev * pPrev->dy() * pTau);
	    pRHS->addTerm(-gammaPtr * vxPrev->dx() * pPrev * pTau);
	    pRHS->addTerm(-gammaPtr * vyPrev->dy() * pPrev * pTau);
	  }
	// velocity equation
	vBF->addTerm(invDt * vx, vxTau);
	vBF->addTerm(invDt * vy, vyTau);
	// pressure term
	vRHS->addTerm(-invGammaM2 * pPrev->dx() * vxTau);
	vRHS->addTerm(-invGammaM2 * pPrev->dy() * vyTau);
	// viscoisty term
/*
	vBF->addTerm(invRe * centerPtr * vx->grad(), vxTau->grad());
	vBF->addTerm(invRe * centerPtr * vy->grad(), vyTau->grad());
	vRHS->addTerm(-invRe * vxPrev->grad() * vxTau->grad());
	vRHS->addTerm(-invRe * vyPrev->grad() * vyTau->grad());
	// add boundary flux terms
	vBF->addTerm(centerPtr * invRe * vx->grad() * n,
                     boundaryIndicator * vxTau);
	vBF->addTerm(centerPtr * invRe * vy->grad() * n, 
                     boundaryIndicator * vyTau);
	vRHS->addTerm(-invRe * vxPrev->grad() * n * boundaryIndicator * vxTau);
	vRHS->addTerm(-invRe * vyPrev->grad() * n * boundaryIndicator * vyTau);

*/
	vBF->addTerm(invRe * centerPtr * vx->dx(), vxTau->dx());
	vBF->addTerm(invRe * centerPtr * vx->dy(), vxTau->dy());
	vBF->addTerm(invRe * centerPtr * vy->dx(), vyTau->dx());
	vBF->addTerm(invRe * centerPtr * vy->dy(), vyTau->dy());
	vRHS->addTerm(-invRe * vxPrev->dx() * vxTau->dx());
	vRHS->addTerm(-invRe * vxPrev->dy() * vxTau->dy());
	vRHS->addTerm(-invRe * vyPrev->dx() * vyTau->dx());
	vRHS->addTerm(-invRe * vyPrev->dy() * vyTau->dy());
	// add boundary flux terms
	vBF->addTerm(centerPtr * invRe * vx->dx() * n->x(),
                     boundaryIndicator * vxTau);
	vBF->addTerm(centerPtr * invRe * vx->dy() * n->y(),
                     boundaryIndicator * vxTau);
	vBF->addTerm(centerPtr * invRe * vy->dx() * n->x(), 
                     boundaryIndicator * vyTau);
	vBF->addTerm(centerPtr * invRe * vy->dy() * n->y(), 
                     boundaryIndicator * vyTau);

	vRHS->addTerm(-invRe * vxPrev->dx() * n->x() * boundaryIndicator * vxTau);
	vRHS->addTerm(-invRe * vxPrev->dy() * n->y() * boundaryIndicator * vxTau);
	vRHS->addTerm(-invRe * vyPrev->dx() * n->x() * boundaryIndicator * vyTau);
	vRHS->addTerm(-invRe * vyPrev->dy() * n->y() * boundaryIndicator * vyTau);

	// advection terms  
	vBF->addTerm(centerPtr * vxEq * vx->dx(), vxTau);
	vBF->addTerm(centerPtr * vyEq * vx->dy(), vxTau);
	vBF->addTerm(centerPtr * vxEq * vy->dx(), vyTau);
	vBF->addTerm(centerPtr * vyEq * vy->dy(), vyTau);
	vBF->addTerm(centerPtr * vx * vxEq->dx(), vxTau);
	vBF->addTerm(centerPtr * vy * vxEq->dy(), vxTau);
	vBF->addTerm(centerPtr * vx * vyEq->dx(), vyTau);
	vBF->addTerm(centerPtr * vy * vyEq->dy(), vyTau);
	vRHS->addTerm(-vxEq * vxPrev->dx() * vxTau);
	vRHS->addTerm(-vyEq * vxPrev->dy() * vxTau);
	vRHS->addTerm(-vxEq * vyPrev->dx() * vyTau);
	vRHS->addTerm(-vyEq * vyPrev->dy() * vyTau);
	vRHS->addTerm(-vxPrev * vxEq->dx() * vxTau);
	vRHS->addTerm(-vyPrev * vxEq->dy() * vxTau);
	vRHS->addTerm(-vxPrev * vyEq->dx() * vyTau);
	vRHS->addTerm(-vyPrev * vyEq->dy() * vyTau);
	if (nonlinear)
	  {
	    vBF->addTerm(centerPtr * vxPrev * vx->dx(), vxTau);
	    vBF->addTerm(centerPtr * vyPrev * vx->dy(), vxTau);
	    vBF->addTerm(centerPtr * vxPrev * vy->dx(), vyTau);
	    vBF->addTerm(centerPtr * vyPrev * vy->dy(), vyTau);
	    vBF->addTerm(centerPtr * vx * vxPrev->dx(), vxTau);
	    vBF->addTerm(centerPtr * vy * vxPrev->dy(), vxTau);
	    vBF->addTerm(centerPtr * vx * vyPrev->dx(), vyTau);
	    vBF->addTerm(centerPtr * vy * vyPrev->dy(), vyTau);
	    vRHS->addTerm(-vxPrev * vxPrev->dx() * vxTau);
	    vRHS->addTerm(-vyPrev * vxPrev->dy() * vxTau);
	    vRHS->addTerm(-vxPrev * vyPrev->dx() * vyTau);
	    vRHS->addTerm(-vyPrev * vyPrev->dy() * vyTau);
	  }
      break;
      case feType::pDPG:
// todo
      break;
      case feType::uwDPG:
	pBF->addTerm(invDt * p,pTau);
	pBF->addTerm(centerPtr * vxEq * gradP->x(), pTau);
	pBF->addTerm(centerPtr * vyEq * gradP->y(), pTau);
	pRHS->addTerm(-vxEq * gradPPrev->x() * pTau);
	pRHS->addTerm(-vyEq * gradPPrev->y() * pTau);
	pRHS->addTerm(-vxPrev * pEq->dx() * pTau);
	pRHS->addTerm(-vyPrev * pEq->dy() * pTau);
	pRHS->addTerm(-gammaPtr * gradVxPrev->x() * pEq * pTau);
	pRHS->addTerm(-gammaPtr * gradVyPrev->y() * pEq * pTau);
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
	// viscoisty term
	vBF->addTerm(invRe * centerPtr * gradVx->x(), vxTau->dx());
	vBF->addTerm(invRe * centerPtr * gradVx->y(), vxTau->dy());
	vBF->addTerm(invRe * centerPtr * gradVy->x(), vyTau->dx());
	vBF->addTerm(invRe * centerPtr * gradVy->y(), vyTau->dy());
	vRHS->addTerm(-invRe * gradVxPrev->x() * vxTau->dx());
	vRHS->addTerm(-invRe * gradVxPrev->y() * vxTau->dy());
	vRHS->addTerm(-invRe * gradVyPrev->x() * vyTau->dx());
	vRHS->addTerm(-invRe * gradVyPrev->y() * vyTau->dy());
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
      break;
    };

/******************* Define the inner product ***************************/
  switch (thisFeMethod)
    {
      case feType::CG:
	vIp = Teuchos::null;
	pIp = Teuchos::null;
      break;
      case feType::pDPG:
// todo
      break;
      case feType::uwDPG:
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
      break;
    };

/******************* Initialize perturbed fields ***************************/
  int solutionOrdinal = 0;
  map<int, FunctionPtr> vPert;
  map<int, FunctionPtr> pPert;
  double kxPert = 2 * M_PI/(xDim * nx) ;
  double wdt = Mach * dt / (2 * kxPert); // pressire is at 1/2 time step
  FunctionPtr vxPertFunc = Teuchos::rcp(new coswtkx(kxPert, 0));
  FunctionPtr gradVxPertFunc = Teuchos::rcp(new sinwtkx(kxPert, 0));
  FunctionPtr pPertFunc = Teuchos::rcp(new coswtkx(kxPert, wdt));
  FunctionPtr gradPPertFunc = Teuchos::rcp(new sinwtkx(kxPert, wdt));
  vPert[vx->ID()] = amp /(gamma * Mach) * vxPertFunc ;
  vPert[vy->ID()] = zero;
  pPert[p->ID()] = amp * pPertFunc;
  switch (thisFeMethod)
    {
      case feType::CG:
// do nothing
      break;
      case feType::pDPG:
// todo
      break;
      case feType::uwDPG:
	vPert[gradVx->ID()] = 
	  Function::vectorize(kxPert * amp /(gamma * Mach) * gradVxPertFunc,
			      zero);
	vPert[gradVy->ID()] = Function::vectorize(zero,zero);
	pPert[gradP->ID()] =
          Function::vectorize(kxPert * amp * gradPPertFunc,zero);
      break;
    };
  vPrevTime->projectOntoMesh(vPert, solutionOrdinal);
  pPrevTime->projectOntoMesh(pPert, solutionOrdinal);

/******************* Apply the boundary conditions ***************************/
//  I have a doubly periodic box, so I don't have to do anything here

/******************* Define the soln  pointers ***************************/
  SolutionPtr pSoln = Solution::solution(pBF, pMesh, pBc, pRHS, pIp);
  SolutionPtr vSoln = Solution::solution(vBF, vMesh, vBc, vRHS, vIp);  
  FunctionPtr vxIncr = Function::solution(vx, vSoln);
  FunctionPtr vyIncr = Function::solution(vy, vSoln);
  FunctionPtr pIncr = Function::solution(p, pSoln);
  FunctionPtr gradVxIncr, gradVyIncr, gradPIncr;
  FunctionPtr vL2Sqr, pL2Sqr;

  switch (thisFeMethod)
    {
      case feType::CG:
	vL2Sqr = vxIncr * vxIncr + vyIncr * vyIncr;
	pL2Sqr = pIncr * pIncr;
      break;
      case feType::pDPG:
// todo
      break;
      case feType::uwDPG:
	gradVxIncr = Function::solution(gradVx, vSoln);
	gradVyIncr = Function::solution(gradVy, vSoln);
	gradPIncr = Function::solution(gradP, pSoln);
	vL2Sqr = vxIncr * vxIncr + vyIncr * vyIncr 
	       + gradVxIncr * gradVxIncr + gradVyIncr * gradVyIncr;
	pL2Sqr = pIncr * pIncr + gradPIncr * gradPIncr;
      break;
    };

/******************* Solve  ***************************/
  double vL2Incr = 0;
  double pL2Incr = 0;
  HDF5Exporter vExporter(vMesh, "SoundWavesFlow");
  HDF5Exporter pExporter(pMesh, "SoundWavesPressure");
  int numSubdivisions = 4;
  int pltnum = 1;
  int step = 0;

/* satic condensation does not work with continuous g  
  pSoln->setUseCondensedSolve(true);
  vSoln->setUseCondensedSolve(true);
*/
  vExporter.exportSolution(vPrevTime, 0, numSubdivisions);
  pExporter.exportSolution(pPrevTime, 0, numSubdivisions);
  cout << "before solve loop" << endl;
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
