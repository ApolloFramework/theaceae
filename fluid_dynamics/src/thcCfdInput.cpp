/*
  $Id: thcCFDInput.cpp 26 2018-06-08 19:30:25Z ehowell $
  Use Teuchos routines to configure all of the RANS cases
  This is the entire parameter list.  Nalu puts the individual
  control with the classes for better localization and an argument 
  could be made that this is a better appraoch, but lacking documentation
  it is easier to have it in a single file and use sublists.
 */

#include "thcCfdInput.h"

using namespace std;
using namespace Teuchos;

//  Using Teuchos utilities, enable reading in an input file (assumed to be
//  yaml, but could be XML) and putting into a ParameterList.  Enable
//  command-line options to specify input filename and certain parameter
//  overrides.  Finally, the final ParameterList is written out for provenance
//  purposes.
//
namespace Theaceae
{

ParameterList thcCFDInput(int argc, char* argv[] )
{
  using std::cout;
  using std::endl;


  try {

    /******************* Input validators *************************************/
    // Create a validator for the "ispace" parameter
    /*  TODO: 
    RCP<StringToIntegralParameterEntryValidator<string>
      ipSpaceValidator = rcp(
      new StringToIntegralParameterEntryValidator<std::string>(
         tuple<std::string>("simple","chan","graph")
         "ipSpace", 
        )
      );
  
    */

    /******************* Set up Parameters in Parameter List*******************/
    //  Convention: explicitly type
    
    ParameterList thcPL;

    // Mesh options
    ParameterList& PLmesh = thcPL.sublist("mesh");
    PLmesh.set<int>("mx", 21, "Number of Elements in x direction");
    PLmesh.set<int>("my", 20, "Number of Elements in y direction");
    PLmesh.set<int>("H1Order", 4, "Polynomial Degree");
    PLmesh.set<int>("delta_k", 2, "Polynomial Degree Enhancement");

    PLmesh.set<double>("plateLength", 2.0, "Length of Plate");
    PLmesh.set<double>("plateStart", 0.33333333333333, "Distance before the start of th plate");
    PLmesh.set<int>("yDim", 1., "Domain Height");
    PLmesh.set<bool>("packMesh", true, "Option to turn on Mesh Packing");
    PLmesh.set<double>("xTau", 4.01, "x Direction backing factor");   // was 0.01
    PLmesh.set<double>("yBeta", 1.0025, "y direction Mesh Packing Factor");
                                                //was 1.025 1.000025 is need for y1 ~8e-6 
    PLmesh.set<double>("xBetaL", 1.009, "x Mesh Packing Factor left of plate");
    PLmesh.set<double>("xBetaR", 1.025, "x Mesh Packing Factor right of plate");
    PLmesh.set<int>("nXL", 3, "Number of Elements to pack left of plate");

    // time step options
    ParameterList& PLtime = thcPL.sublist("timeControl");
    PLtime.set<double>("centerParam", 1.0, "Centering parameter for advection, diffusion, and pressure terms");
    PLtime.set<double>("sourceCenter", 0.0, "Centering parameter for source terms");
    PLtime.set<double>("dt", 1e-2, "time step"); // do something fancier in the future CFL ~ u_xdt/dx~ 1/2
    PLtime.set<int>("maxsteps", 20, "number of steps to take");

    // check pointing options 
    ParameterList& PLdump = thcPL.sublist("dump");
    PLdump.set("filename", "flatplate", "Base name of out file");
    PLdump.set<int>("ndum", 5, "Frequency of output plotting files");
    PLdump.set<int>("nsave", 5000, "Frequency of outputing check pointing files");
    PLdump.set<bool>("readFile", false, "Flag to restart from file");
    PLdump.set<string>("readFileName", "../18061401/flateplate3000", "File name to read from");

     // physics input
    ParameterList& PLphys = thcPL.sublist("physics");
    PLphys.set<double>("Re", 1.e6, "Reynolds Number");
    PLphys.set<int>("sourceStep0", -1000, "First step to apply turbulent source term");
    PLphys.set<int>("sourceStepf", -1100, "Step after which the turbulent source terms have their full magnitude");
    PLphys.set<string>("ipSpace", "simple", "Inner product space: simple, chan, graph");
    //PLphys.set<string>("ipSpace", "simple", "Inner product space: simple, chan, graph",ipSpaceValidator);
     // initialization input but put into physics
    PLphys.set<double>("flowVelocity", 1.0, "input flow Velocity");

    // turbulent model inputs
    ParameterList& PLturb = thcPL.sublist("turbulence");
    PLturb.set<bool>("useTurbModel", true, "Use Turbulance Model");
    PLturb.set<bool>("useSANeg", true, "Use SA-neg model"); //When true use the SA-neg model (see NASA website)
    PLturb.set<bool>("useSBar",  true, "Use S bar"); // When true use SBar not STilde (see NASA website)
    PLturb.set<double>("chiFsBc", 3.0, "chi freesteam bc");
    PLturb.set<bool>("splitTurbAdv", true, "split turbulence advance");
    PLturb.set<int>("nTurbSteps", 100, "number of steps in turbulant advance per fluid advance");

    /******************* Command-line parsing *********************************/
    // Two things:  
    //    1) get filename to parse for parameter list, 
    //    2) allow some variables to be overriden by command-line arguments
    
    std::string    thcInFileName = "theaceae.yaml";
    std::string    thcOutFileName = "theaceae.xml";

    CommandLineProcessor  clp(false); // Don't throw exceptions
    clp.setOption("inputFile",  &thcInFileName,  "Theaceae input file (yaml format)");
    clp.setOption("outputFile", &thcOutFileName, "The XML file to debug yaml file");
    clp.setDocString("Input file for theaceae");

    // Enable options that can override parameters in the file
    //const ParameterEntry& lpl_ms = getConst(PLtime).getEntryPtr("maxsteps");
    const ParameterEntry& lpl_ms = getConst(thcPL).sublist("timeControl").getEntry("maxsteps");
    int clp_maxsteps = getValue<int>(lpl_ms);
    string maxsteps_doc = lpl_ms.docString();

    // clp.setOption("maxsteps", &clp_maxsteps);
    clp.setOption("maxsteps", &clp_maxsteps, maxsteps_doc.c_str());

    // Now parse the command line arguments
    CommandLineProcessor::EParseCommandLineReturn
      parse_return = clp.parse(argc,argv);
    if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
      std::cout << "\nParse command-line arguments failed " << std::endl;
      //TODO return parse_return;
    }

    /******************* Parse file and update from defaults ******************/
    updateParametersFromXmlFile(thcInFileName, inoutArg(thcPL));

    // Enable options that can override parameters in the file
    //PLtime.set<int>("maxsteps", &clp_maxsteps, &maxsteps_doc);
    PLtime.set<int>("maxsteps", clp_maxsteps, "number of steps to take");
    //PLtime.set<int>("maxsteps", &clp_maxsteps, "number of steps to take");

    // State the norm
    string lpl_ips = getConst(thcPL).sublist("physics").get<string>("ipSpace");
    string ipType = lpl_ips + "Norm";
    PLphys.set<string>("ipType", ipType, "Inner-product norm to use");

    // Finally write the final ParameterList to a file for good provenance
    // practices
    writeParameterListToXmlFile(thcPL,thcOutFileName);

    return thcPL;
  }
  catch(const std::exception& e)
  {
      std::cerr << "caught std::exception:\n\n";
      OSTab tab(std::cerr);
      std::cerr << e.what() << std::endl;
  }

};
}
