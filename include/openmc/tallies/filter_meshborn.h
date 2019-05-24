#ifndef OPENMC_TALLIES_FILTER_MESHBORN_H
#define OPENMC_TALLIES_FILTER_MESHBORN_H

#include <cstdint>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Indexes the location of particle events to a regular mesh.  For tracklength
//! tallies, it will produce multiple valid bins and the bin weight will
//! correspond to the fraction of the track length that lies in that bin.
//==============================================================================

class MeshbornFilter : public Filter
{
public:
  ~MeshbornFilter() = default;

  std::string type() const override {return "meshborn";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle* p, int estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  virtual int32_t mesh() const {return mesh_;}

  virtual void set_mesh(int32_t mesh);

protected:
  int32_t mesh_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_MESH_H
