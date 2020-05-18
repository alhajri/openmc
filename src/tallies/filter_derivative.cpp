#include "openmc/tallies/filter_derivative.h"

#include <fmt/core.h>

#include "openmc/capi.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/mesh.h"
#include "openmc/xml_interface.h"
#include "openmc/nuclide.h"

namespace openmc {

//==============================================================================
// DerivativeFilter implementation
//==============================================================================

void
DerivativeFilter::from_xml(pugi::xml_node node)
{
  // Set the importance (adjoint flux)
  auto nuclides = get_node_array<std::string>(node, "nuclides");
  nuclides_.clear();

  for (const auto& nuc : nuclides) {
    if (nuc == "total") {
      nuclides_.push_back(-1);
    } else {
      auto search = data::nuclide_map.find(nuc);
      if (search == data::nuclide_map.end())
        fatal_error(fmt::format("Could not find the nuclide {} specified in "
          "tally {} in any material", nuc, id_));
      nuclides_.push_back(search->second);
      
    }
  }

  // Set the size
  const auto& nuc {*data::nuclides[nuclides_[0]]};
  n_bins_ = nuc.multipole_->data_.shape()[0] * nuc.multipole_->data_.shape()[1] * 2;
}

void DerivativeFilter::set_score(std::vector<int> scores)
{
    score_.clear();
    for (int i = 0; i < scores.size(); ++i){
        score_.push_back(scores[i]);
    }
}

void
DerivativeFilter::get_all_bins(const Particle& p, TallyEstimator estimator, FilterMatch& match)
const
{

  if (nuclide_[0] >= 0) {
    auto j = model::materials[p.material_]->mat_nuclide_index_[nuclide_[0]];
    if (j == C_NONE) return;
  }
  const auto& nuc {*data::nuclides[nuclide_[0]]};
  if (!multipole_in_range(nuc, p.E_)) return;
  double sig_s, sig_a, sig_f;
  std::tie(sig_s, sig_a, sig_f)
    = nuc.multipole_->evaluate(p.E_, p.sqrtkT_);
  double sigma; //microscopic cross section
  std::pair<int, std::vector<double>> derivative;

  // Get the correct cross section
  switch (score_[0]) {
  case SCORE_TOTAL:
    sigma = sig_s + sig_a;
    derivative = nuc.multipole_->evaluate_pole_deriv_total(p.E_, p.sqrtkT_);
    break;
  case SCORE_SCATTER:
    sigma = sig_s;
    derivative = nuc.multipole_->evaluate_pole_deriv_scatter(p.E_, p.sqrtkT_);
    break;
  case ELASTIC:
    sigma = sig_s;
    derivative = nuc.multipole_->evaluate_pole_deriv_scatter(p.E_, p.sqrtkT_);
    break;
  case SCORE_ABSORPTION:
    sigma = sig_a;
    derivative = nuc.multipole_->evaluate_pole_deriv_absorption(p.E_, p.sqrtkT_);
    break;
  case SCORE_FISSION:
    sigma = sig_f;
    derivative = nuc.multipole_->evaluate_pole_deriv_fission(p.E_, p.sqrtkT_);
    break;
  }

  // the score is atom_density * derivative_total * distance
  int start = derivative.first;
  int size  = derivative.second.size();

  for (int deriv_idx = start; deriv_idx < start + size ; deriv_idx++){
    match.bins_.push_back(deriv_idx);
    match.weights_.push_back(derivative.second[deriv_idx - start] / sigma);
  }
}

void
DerivativeFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  write_dataset(filter_group, "nuclide", nuclide_);
  write_dataset(filter_group, "score", score_);
}

std::string
DerivativeFilter::text_label(int bin) const
{
  return "derivative";
}

}// namespace openmc