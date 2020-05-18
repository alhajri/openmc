#ifndef OPENMC_TALLIES_FILTER_DERIVATIVE_H
#define OPENMC_TALLIES_FILTER_DERIVATIVE_H

#include <vector>
#include <cstdint>

#include <gsl/gsl>

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Multiplies each tally by the importance corresponding to the mesh index 
//==============================================================================

class DerivativeFilter : public Filter
{
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ~DerivativeFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type() const override {return "derivative";}

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match)
  const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;

  //----------------------------------------------------------------------------
  // Accessors=
  void set_score(std::vector<int> scores)

protected:
  //----------------------------------------------------------------------------
  // Data members

  std::vectot<int> nuclide_;

  std::vectot<int> score_;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_IMPORTANCE_H