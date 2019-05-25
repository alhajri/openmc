#include "openmc/tallies/filter_meshborn.h"

#include <sstream>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/mesh.h"
#include "openmc/xml_interface.h"

namespace openmc {

void
MeshbornFilter::get_all_bins(const Particle* p, int estimator, FilterMatch& match)
const
{
 auto bin = model::meshes[mesh_]->get_bin(p->r_born_);
 if (bin >= 0) {
   match.bins_.push_back(bin);
   match.weights_.push_back(1.0);
 }
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_meshborn_filter_get_mesh(int32_t index, int32_t* index_mesh)
{return openmc_mesh_filter_get_mesh(index, index_mesh);}

extern "C" int
openmc_meshborn_filter_set_mesh(int32_t index, int32_t index_mesh)
{return openmc_mesh_filter_set_mesh(index, index_mesh);}

} // namespace openmc
