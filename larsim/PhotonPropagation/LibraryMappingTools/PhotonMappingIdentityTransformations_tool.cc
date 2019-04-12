/**
 * @file   larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations_tool.cc
 * @brief  A photon mapping identity transformation: toolification.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   March 20, 2019
 * @see    `larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations.h`
 *         `larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations.cxx`
 * 
 * Note that the separation of the code into a implementation and a tool part
 * is not because of factorization to have art-independent code, but because
 * we want this object to be customizable by inheritance.
 * An _art_ tool derived from `phot::PhotonMappingIdentityTransformations` class
 * must link with the library where its implementation is defined. When that
 * library also contains the registration of the tool into the _art_ framework
 * (`DEFINE_ART_CLASS_TOOL()` macro), the functions registering the base tool
 * will conflict with the ones for the registration of the derived tool.
 * 
 */


// LArSoft libraries
#include "larsim/PhotonPropagation/LibraryMappingTools/PhotonMappingIdentityTransformations.h"

// framework libraries
#include "art/Utilities/ToolMacros.h"


//------------------------------------------------------------------------------
DEFINE_ART_CLASS_TOOL(phot::PhotonMappingIdentityTransformations)


//------------------------------------------------------------------------------

