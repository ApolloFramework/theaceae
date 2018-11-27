/*
  $Id: thcCFDInput.cpp 26 2018-06-08 19:30:25Z ehowell $
  Use Teuchos routines to configure all of the RANS cases
  This is the entire parameter list.  Nalu puts the individual
  control with the classes for better localization and an argument 
  could be made that this is a better appraoch, but lacking documentation
  it is easier to have it in a single file and use sublists.
 */

#ifndef THEACEAE_THCCFDINPUT_H
#define THEACEAE_THCCFDINPUT_H

/*
 *  thcCFDInput.h
 *
 *  Created by Scott Kruger on 9/10/18.
 *
 */

#include "Teuchos_getConst.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ParameterList.hpp" 
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#ifdef HAVE_TEUCHOS_EXTENDED
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_YamlParser_decl.hpp"
#include "Teuchos_YamlParameterListCoreHelpers.hpp"
#endif


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
  /** \brief Yaml and CLI parsing for RANS Turbulence model runs
   */
  int thcCFDInput(int argc, char* argv[], ParameterList &thcPL );
}
#endif
