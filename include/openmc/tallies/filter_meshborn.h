#ifndef OPENMC_TALLIES_FILTER_MESHBORN_H
#define OPENMC_TALLIES_FILTER_MESHBORN_H

#include <cstdint>

#include "openmc/tallies/filter_mesh.h"

namespace openmc {

//==============================================================================
//! Indexes the location of particle events to a regular mesh.  For tracklength
//! tallies, it will produce multiple valid bins and the bin weight will
//! correspond to the fraction of the track length that lies in that bin.
//==============================================================================

class MeshbornFilter : public MeshFilter
{
public:
  std::string type() const override {return "meshborn";}

  void get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match)
  const override;

  std::string text_label(int bin) const override;

};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESH_H
