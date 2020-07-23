#include "openmc/tallies/sensitivity.h"

#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/source.h"
#include "openmc/search.h"
#include "openmc/tallies/tally.h"
#include "openmc/xml_interface.h"
#include "openmc/message_passing.h"
#include "openmc/tallies/filter_importance.h"

#include <fmt/core.h>
#include "xtensor/xadapt.hpp"
#include "xtensor/xbuilder.hpp" // for empty_like
#include "xtensor/xview.hpp"

template class std::vector<openmc::TallySensitivity>;

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  std::vector<TallySensitivity> tally_sens;
  std::unordered_map<int, int> tally_sens_map;
}

//==============================================================================
// SensitivityTally implementation
//==============================================================================

SensitivityTally::SensitivityTally(int32_t id)
{
  index_ = model::tallies.size(); // Avoids warning about narrowing
  this->set_id(id);
  this->set_filters({});
}

SensitivityTally::SensitivityTally(pugi::xml_node node)
{
  index_ = model::tallies.size(); // Avoids warning about narrowing

  // Copy and set tally id
  if (!check_for_node(node, "id")) {
    throw std::runtime_error{"Must specify id for tally in tally XML file."};
  }
  int32_t id = std::stoi(get_node_value(node, "id"));
  this->set_id(id);

  if (check_for_node(node, "name")) name_ = get_node_value(node, "name");

  // =======================================================================
  // READ DATA FOR FILTERS

  // Check if user is using old XML format and throw an error if so
  if (check_for_node(node, "filter")) {
    throw std::runtime_error{"Tally filters must be specified independently of "
      "tallies in a <filter> element. The <tally> element itself should "
      "have a list of filters that apply, e.g., <filters>1 2</filters> "
      "where 1 and 2 are the IDs of filters specified outside of "
      "<tally>."};
  }

  // Determine number of filters
  std::vector<int> filter_ids;
  if (check_for_node(node, "filters")) {
    filter_ids = get_node_array<int>(node, "filters");
  }

  // Allocate and store filter user ids
  std::vector<Filter*> filters;
  for (int filter_id : filter_ids) {
    // Determine if filter ID is valid
    auto it = model::filter_map.find(filter_id);
    if (it == model::filter_map.end()) {
      throw std::runtime_error{fmt::format(
        "Could not find filter {} specified on tally {}", filter_id, id_)};
    }

    // Store the index of the filter
    filters.push_back(model::tally_filters[it->second].get());
  }

  // Set the filters
  this->set_filters(filters);

  // =======================================================================
  // READ DATA FOR NUCLIDES

  // Sensitivity only allows total material rates.

  nuclides_.clear();

  // By default, we tally just the total material rates.
  if (!check_for_node(node, "nuclides")) {
    nuclides_.push_back(-1);
  } else {
    this->set_nuclides(node);
  }

  // =======================================================================
  // READ DATA FOR SCORES
  // Sensitivity only allows nu-fission score

  // Reset state and prepare for the new scores.
  if (!check_for_node(node, "scores")){
    scores_.clear();
    depletion_rx_ = false;
    scores_.reserve(1);

    auto score = SCORE_NU_FISSION;

    scores_.push_back(score);
  } else {
    this->set_scores(node);
  }

  // =======================================================================
  // READ DATA FOR SENSITIVITIES

  if (!check_for_node(node, "sensitivity")) {
    fatal_error(fmt::format("No sensitivity specified on tally {}.", id_));
  }

   // Check for a tally sensitivity.
  if (check_for_node(node, "sensitivity")) {
    int sens_id = std::stoi(get_node_value(node, "sensitivity"));

    // Find the sensitivity with the given id, and store it's index.
    auto it = model::tally_sens_map.find(sens_id);
    if (it == model::tally_sens_map.end()) {
      fatal_error(fmt::format(
        "Could not find sensitivity {} specified on tally {}", sens_id, id_));
    }

    sens_ = it->second;
    //sens_ = it->first;

    // Only analog or collision estimators are supported for sensitivity
    // tallies.
    if (estimator_ == TallyEstimator::TRACKLENGTH) {
      estimator_ = TallyEstimator::COLLISION;
    }
  }


  // If settings.xml trigger is turned on, create tally triggers
  if (settings::trigger_on) {
    this->init_triggers(node);
  }

  // check for gpt
  if (check_for_node(node, "gpt")) {
    gpt_ = true;
  }

  // =======================================================================
  // SET TALLY ESTIMATOR

  // Check if user specified estimator
  if (check_for_node(node, "estimator")) {
    std::string est = get_node_value(node, "estimator");
    if (est == "analog") {
      estimator_ = TallyEstimator::ANALOG;
    } else if (est == "collision") {
      // If the estimator was set to an analog estimator, this means the
      // tally needs post-collision information
      if (estimator_ == TallyEstimator::ANALOG) {
        throw std::runtime_error{fmt::format("Cannot use collision estimator "
          "for tally ", id_)};
      }

      // Set estimator to collision estimator
      estimator_ = TallyEstimator::COLLISION;

    } else {
      throw std::runtime_error{fmt::format(
        "Invalid estimator '{}' on tally {}", est, id_)};
    }
  }
}

SensitivityTally::~SensitivityTally()
{
  model::tally_map.erase(id_);
}

//SensitivityTally*
//SensitivityTally::create(int32_t id)
//{
//  model::tallies.push_back(std::make_unique<SensitivityTally>(id));
//  return model::tallies.back().get();
//}

void
SensitivityTally::set_filters(gsl::span<Filter*> filters)
{
  // Clear old data.
  filters_.clear();
  strides_.clear();

  // Copy in the given filter indices.
  auto n = filters.size();
  if (n != 1) {
      throw std::runtime_error{fmt::format("Cannot use more than one filter for sensitivity")};
  }
  filters_.reserve(n);

  for (int i = 0; i < n; ++i) {
    // Add index to vector of filters
    auto& f {filters[i]};
    filters_.push_back(model::filter_map.at(f->id()));

    // Keep track of indices for special filters.
    if (!dynamic_cast<const ImportanceFilter*>(f)) {
      throw std::runtime_error{fmt::format("Must use an importance filter for sensitivity")};
    }
  }

  // Set the strides.  Filters are traversed in reverse so that the last filter
  // has the shortest stride in memory and the first filter has the longest
  // stride.
  strides_.resize(n, 0);
  int stride = 1;
  for (int i = n-1; i >= 0; --i) {
    strides_[i] = stride;
    stride *= model::tally_filters[filters_[i]]->n_bins();
  }
  n_filter_bins_ = stride;
}

void SensitivityTally::init_results()
{
  int n_scores = 1;
  int n_filter_bins = model::tally_sens[sens_].n_bins_;
  results_ = xt::empty<double>({n_filter_bins, n_scores, 3});
  previous_results_ = xt::empty<double>({n_filter_bins, n_scores, 1});
}

void SensitivityTally::reset()
{
  n_realizations_ = 0;
  if (results_.size() != 0) {
    xt::view(results_, xt::all()) = 0.0;
    xt::view(previous_results_, xt::all()) = 0.0;
  }
}

void SensitivityTally::accumulate()
{
  // Increment number of realizations
  n_realizations_ += settings::reduce_tallies ? 1 : mpi::n_procs;

  if (mpi::master || !settings::reduce_tallies) {
    // Calculate total source strength for normalization
    double total_source = 0.0;
    if (settings::run_mode == RunMode::FIXED_SOURCE) {
      for (const auto& s : model::external_sources) {
        total_source += s.strength();
      }
    } else {
      total_source = 1.0;
    }

    // Account for number of source particles in normalization
    if (gpt_) {
      double norm = total_source / (settings::n_particles * settings::gen_per_batch);
      denominator_ = norm;
    }

    // Accumulate each result
    // TODO: ignore the first realization
    for (int i = 0; i < results_.shape()[0]; ++i) {
      for (int j = 0; j < results_.shape()[1]; ++j) {
        double val = previous_results_(i, j, SensitivityTallyResult::VALUE)/ denominator_;
        results_(i, j, SensitivityTallyResult::SUM) += val;
        results_(i, j, SensitivityTallyResult::SUM_SQ) += val*val;
      }
    }

    denominator_ = 0.0; //<! zero out denom

    // Move current value to previous value and zero out each result
    for (int i = 0; i < results_.shape()[0]; ++i) {
      for (int j = 0; j < results_.shape()[1]; ++j) {
        double val = results_(i, j, SensitivityTallyResult::VALUE);
        results_(i, j, SensitivityTallyResult::VALUE) = 0.0;
        previous_results_(i, j, SensitivityTallyResult::VALUE) = val;
      }
    }
  }
}

//==============================================================================
// TallySensitivity implementation
//==============================================================================

TallySensitivity::TallySensitivity(pugi::xml_node node)
{
  if (check_for_node(node, "id")) {
    id = std::stoi(get_node_value(node, "id"));
  } else {
    fatal_error("Must specify an ID for <sensitivity> elements in the tally "
                "XML file");
  }

  if (id <= 0)
    fatal_error("<sensitivity> IDs must be an integer greater than zero");

  std::string variable_str = get_node_value(node, "variable");

  std::string nuclide_name = get_node_value(node, "nuclide");
  bool found = false;
  for (auto i = 0; i < data::nuclides.size(); ++i) {
    if (data::nuclides[i]->name_ == nuclide_name) {
      found = true;
      sens_nuclide = i;
    }
  }
  if (!found) {
    fatal_error(fmt::format("Could not find the nuclide \"{}\" specified in "
      "derivative {} in any material.", nuclide_name, id));
  }

  if (variable_str == "cross_section") {
    variable = SensitivityVariable::CROSS_SECTION;

    // ADD LOGIC TO LOOK FOR ENERGY BINS OTHERWISE DEFAULT TO MAX AND MIN NEUTRON ENERGY
    auto bins = get_node_array<double>(node, "energy");
    this->set_bins(bins);

    // ADD LOGIC TO SET THE SCORE TYPE
    auto reaction = get_node_value(node, "reaction");
    sens_reaction = score_str_to_int(reaction);

  } else if (variable_str == "multipole") {
    variable = SensitivityVariable::MULTIPOLE;

    // check if curvefit sensitivities were asked for, maybe move to a different variable

    // ADD LOGIC TO SET THE BINS TO SIZE OF MULTIPOLE PARAMETERS
    // set n_bins_ 
    const auto& nuc {*data::nuclides[sens_nuclide]};
    n_bins_ = nuc.multipole_->data_.shape()[0] * nuc.multipole_->data_.shape()[1] * 2;

  } else if (variable_str == "curve_fit") {
    variable = SensitivityVariable::CURVE_FIT;

    // check if curvefit sensitivities were asked for, maybe move to a different variable

    // ADD LOGIC TO SET THE BINS TO SIZE OF MULTIPOLE PARAMETERS
    // set n_bins_ 
    const auto& nuc {*data::nuclides[sens_nuclide]};
    n_bins_ = nuc.multipole_->curvefit_.shape()[0] * nuc.multipole_->curvefit_.shape()[1] * nuc.multipole_->curvefit_.shape()[2];

  } else {
    fatal_error(fmt::format("Unrecognized variable \"{}\" on derivative {}",
      variable_str, id));
  }

}

void
TallySensitivity::set_bins(gsl::span<const double> bins)
{
  // Clear existing bins
  energy_bins_.clear();
  energy_bins_.reserve(bins.size());

  // Copy bins, ensuring they are valid
  for (gsl::index i = 0; i < bins.size(); ++i) {
    if (i > 0 && bins[i] <= bins[i-1]) {
      throw std::runtime_error{"Energy bins must be monotonically increasing."};
    }
    energy_bins_.push_back(bins[i]);
  }

  n_bins_ = energy_bins_.size() - 1;
}

//==============================================================================
// Non-method functions
//==============================================================================

void
read_tally_sensitivities(pugi::xml_node node)
{
  // Populate the sensitivities array.
  for (auto sens_node : node.children("sensitivity"))
    model::tally_sens.emplace_back(sens_node);

  // Fill the sensitivity map.
  for (auto i = 0; i < model::tally_sens.size(); ++i) {
    auto id = model::tally_sens[i].id;
    auto search = model::tally_sens_map.find(id);
    if (search == model::tally_sens_map.end()) {
      model::tally_sens_map[id] = i;
    } else {
      fatal_error("Two or more sensitivities use the same unique ID: "
                  + std::to_string(id));
    }
  }

  // Make sure sensitivities were not requested for an MG run.
  if (!settings::run_CE && !model::tally_sens.empty())
    fatal_error("Sensitivities not supported in multi-group mode");
}

void
score_track_sensitivity(Particle& p, double distance)
{
  // A void material cannot be perturbed so it will not affect flux sensitivities.
  if (p.material_ == MATERIAL_VOID) return;
  const Material& material {*model::materials[p.material_]};

  for (auto idx = 0; idx < model::tally_sens.size(); idx++) {
    const auto& sens = model::tally_sens[idx];
    auto& cumulative_sensitivities = p.cumulative_sensitivities_[idx];

    double atom_density = 0.;
        if (sens.sens_nuclide >= 0) {
          auto j = model::materials[p.material_]->mat_nuclide_index_[sens.sens_nuclide];
          if (j == C_NONE) continue;
          atom_density = model::materials[p.material_]->atom_density_(j);
        }

    switch (sens.variable) {

    case SensitivityVariable::CROSS_SECTION:
    {
      // Calculate the sensitivity with respect to the cross section
      // at this energy

      // Get the post-collision energy of the particle.
      auto E = p.E_;


      // Get the correct cross section
      double macro_xs;
      switch (sens.sens_reaction) {
      case SCORE_TOTAL:
        if (sens.sens_nuclide >=0){
            macro_xs = p.neutron_xs_[sens.sens_nuclide].total * atom_density;
        } else {
            macro_xs = p.macro_xs_.total;
        }
        break;
      case SCORE_SCATTER:
        if (sens.sens_nuclide >=0){
            macro_xs = (p.neutron_xs_[sens.sens_nuclide].total 
            - p.neutron_xs_[sens.sens_nuclide].absorption) * atom_density;
        } else {
            macro_xs = p.macro_xs_.total - p.macro_xs_.absorption;
        }
        break;
      case ELASTIC:
        if (sens.sens_nuclide >= 0) {
            if (p.neutron_xs_[sens.sens_nuclide].elastic == CACHE_INVALID)
              data::nuclides[sens.sens_nuclide]->calculate_elastic_xs(p);
            macro_xs = p.neutron_xs_[sens.sens_nuclide].elastic * atom_density;
          } 
        break;
      case SCORE_ABSORPTION: 
        if (sens.sens_nuclide >=0){
            macro_xs = p.neutron_xs_[sens.sens_nuclide].absorption * atom_density;
        } else {
            macro_xs = p.macro_xs_.absorption;
        }
        break;
      case SCORE_FISSION:
        if (p.macro_xs_.absorption == 0) continue;

        if (sens.sens_nuclide >= 0) {
          macro_xs = p.neutron_xs_[sens.sens_nuclide].fission * atom_density;
        } else {
          macro_xs = p.macro_xs_.fission;
        }
        
        break;
      }
  
      // Bin the energy.
      if (E >= sens.energy_bins_.front() && E <= sens.energy_bins_.back()) {
        auto bin = lower_bound_index(sens.energy_bins_.begin(), sens.energy_bins_.end(), E);
        cumulative_sensitivities[bin] -= distance * macro_xs;
      }

    }
    break;

    case SensitivityVariable::MULTIPOLE:    
    {
      // check if in resonance range
      const auto& nuc {*data::nuclides[sens.sens_nuclide]};
      if (multipole_in_range(nuc, p.E_)){
        // Calculate derivative of the total cross section at p->E_
        auto derivative = nuc.multipole_->evaluate_pole_deriv_total(p.E_, p.sqrtkT_);

        // the score is atom_density * derivative_total * distance
        int start = derivative.first;
        int size  = derivative.second.size();

        double score = atom_density*distance;

        for (int deriv_idx = start; deriv_idx < start + size ; deriv_idx++){
          cumulative_sensitivities[deriv_idx] -= score*derivative.second[deriv_idx - start];
        }
      }
    }
      break;

    case SensitivityVariable::CURVE_FIT:
    {
      // check if in resonance range
      const auto& nuc {*data::nuclides[sens.sens_nuclide]};
      if (multipole_in_range(nuc, p.E_)){
        // Calculate derivative of the total cross section at p->E_
        auto derivative = nuc.multipole_->evaluate_fit_deriv_total(p.E_, p.sqrtkT_);

        // the score is atom_density * derivative_total * distance
        int start = derivative.first;
        int size  = derivative.second.size();

        double score = atom_density*distance;

        for (int deriv_idx = start; deriv_idx < start + size ; deriv_idx++){
          cumulative_sensitivities[deriv_idx] -= score*derivative.second[deriv_idx - start];
        }
      }
    }
      break;
    }
  }
}

void score_collision_sensitivity(Particle& p)
{
  // A void material cannot be perturbed so it will not affect flux derivatives.
  if (p.material_ == MATERIAL_VOID) return;

  // only scattering events effect the cumulative tallies
  if (p.event_ != TallyEvent::SCATTER) return;

  const Material& material {*model::materials[p.material_]};

  for (auto idx = 0; idx < model::tally_sens.size(); idx++) {
    const auto& sens = model::tally_sens[idx];
    auto& cumulative_sensitivities = p.cumulative_sensitivities_[idx];

    if (p.event_nuclide_ != sens.sens_nuclide) continue;
    // Find the index in this material for the diff_nuclide.
    int i;
    for (i = 0; i < material.nuclide_.size(); ++i)
      if (material.nuclide_[i] == sens.sens_nuclide) break;
    // Make sure we found the nuclide.
    if (material.nuclide_[i] != sens.sens_nuclide) {
      fatal_error(fmt::format(
        "Could not find nuclide {} in material {} for tally sensitivity {}",
        data::nuclides[sens.sens_nuclide]->name_, material.id_, sens.id));
    }

    switch (sens.variable) {

    case SensitivityVariable::CROSS_SECTION:
    {

      // Get the pre-collision energy of the particle.
      double E = p.E_last_;
      
      // Get the correct cross section
      double score;
      switch (sens.sens_reaction) {
      case SCORE_TOTAL:
        score = 0.0;
        break;
      case SCORE_SCATTER:
        score = 1.0;
        break;
      case ELASTIC:
        if (p.event_mt_ != ELASTIC) continue;
        score = 1.0;
        break;
      case SCORE_ABSORPTION: 
        score = 0.0;
        break;
      case SCORE_FISSION:
        score = 0.0;
        break;
      }

  
      // Bin the energy.
      if (E >= sens.energy_bins_.front() && E <= sens.energy_bins_.back()) {
        auto bin = lower_bound_index(sens.energy_bins_.begin(), sens.energy_bins_.end(), E);
        cumulative_sensitivities[bin] += score;
      }
    }
      break;

    case SensitivityVariable::MULTIPOLE:
    {
      // check if in resonance range
      const auto& nuc {*data::nuclides[sens.sens_nuclide]};
      if (multipole_in_range(nuc, p.E_last_)){
        // Calculate derivative of the scattering cross section at p->E_last_
        const auto& micro_xs {p.neutron_xs_[sens.sens_nuclide]};
        auto derivative = nuc.multipole_->evaluate_pole_deriv_scatter(p.E_last_, p.sqrtkT_);

        // sum/bin 1/micro_sigma_scatter * derivative
        int start = derivative.first;
        int size  = derivative.second.size();

        double scatter = (micro_xs.total - micro_xs.absorption);

        for (int deriv_idx = start; deriv_idx < start + size ; deriv_idx++){
          cumulative_sensitivities[deriv_idx] += derivative.second[deriv_idx - start]/scatter;
        }
      }
    }
      break;

    case SensitivityVariable::CURVE_FIT:
    {
      // check if in resonance range
      const auto& nuc {*data::nuclides[sens.sens_nuclide]};
      if (multipole_in_range(nuc, p.E_last_)){
        // Calculate derivative of the scattering cross section at p->E_last_
        const auto& micro_xs {p.neutron_xs_[sens.sens_nuclide]};
        auto derivative = nuc.multipole_->evaluate_fit_deriv_scatter(p.E_last_, p.sqrtkT_);

        // sum/bin 1/micro_sigma_scatter * derivative
        int start = derivative.first;
        int size  = derivative.second.size();

        double scatter = (micro_xs.total - micro_xs.absorption);

        for (int deriv_idx = start; deriv_idx < start + size ; deriv_idx++){
          cumulative_sensitivities[deriv_idx] += derivative.second[deriv_idx - start]/scatter;
        }
      }
    }
      break;
    }
  }
}

}// namespace openmc
