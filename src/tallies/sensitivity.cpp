#include "openmc/tallies/sensitivity.h"

#include "openmc/error.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/tallies/tally.h"
#include "openmc/xml_interface.h"

#include <fmt/core.h>

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
    return;
  }

  // =======================================================================
  // READ DATA FOR SCORES
  // Sensitivity only allows nu-fission score

  // Reset state and prepare for the new scores.

  scores_.clear();
  depletion_rx_ = false;
  scores_.reserve(scores.size());

  auto score = SCORE_NU_FISSION;

  scores_.push_back(score);

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

SensitivityTally*
SensitivityTally::create(int32_t id)
{
  model::tallies.push_back(std::make_unique<SensitivityTally>(id));
  return model::tallies.back().get();
}

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
  results_ = xt::empty<double>({n_filter_bins, n_scores, 4});
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
    double norm = total_source / (settings::n_particles * settings::gen_per_batch);

    // Accumulate each result
    for (int i = 0; i < results_.shape()[0]; ++i) {
      for (int j = 0; j < results_.shape()[1]; ++j) {
        double val = results_(i, j, SensitivityTally::PREVIOUS_VALUE) * norm;
        results_(i, j, SensitivityTally::SUM) += val;
        results_(i, j, SensitivityTally::SUM_SQ) += val*val;
      }
    }

    // Move current value to previous value and zero out each result
    for (int i = 0; i < results_.shape()[0]; ++i) {
      for (int j = 0; j < results_.shape()[1]; ++j) {
        double val = results_(i, j, SensitivityTally::VALUE);
        results_(i, j, SensitivityTally::VALUE) = 0.0;
        results_(i, j, SensitivityTally::PREVIOUS_VALUE) = val;
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

  if (variable_str == "cross_section") {
    variable = SensitivityVariable::CROSS_SECTION;

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

    // ADD LOGIC TO LOOK FOR ENERGY BINS OTHERWISE DEFAULT TO MAX AND MIN NEUTRON ENERGY
    auto bins = get_node_array<double>(node, "energy");
    this->set_bins(bins);

    // ADD LOGIC TO SET THE SCORE TYPE
    auto reaction = get_node_value(node, "reaction");
    sens_reaction = score_str_to_int(reaction);

  } else if (variable_str == "multipole") {
    variable = SensitivityVariable::MULTIPOLE;

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

    // ADD LOGIC TO SET THE BINS TO SIZE OF MULTIPOLE PARAMETERS

    // ADD LOGIC TO SET SCORE TO SCORE_TOTAL

  }  else {
    fatal_error(fmt::format("Unrecognized variable \"{}\" on derivative {}",
      variable_str, id));
  }

  sens_material = std::stoi(get_node_value(node, "material"));
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

  n_bins_ = bins_.size() - 1;
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
apply_sensitivity_to_score(const Particle* p, int i_tally, int i_nuclide,
  double atom_density, int score_bin, double& score)
{
  const Tally& tally {*model::tallies[i_tally]};

  if (score == 0.0) return;

  // If our score was previously c then the new score is
  // c * (1/f * d_f/d_p + 1/c * d_c/d_p)
  // where (1/f * d_f/d_p) is the (logarithmic) flux derivative and p is the
  // perturbated variable.

  const auto& deriv {model::tally_derivs[tally.deriv_]};
  const auto flux_deriv = p->flux_derivs_[tally.deriv_];

  // Handle special cases where we know that d_c/d_p must be zero.
  if (score_bin == SCORE_FLUX) {
    score *= flux_deriv;
    return;
  } else if (p->material_ == MATERIAL_VOID) {
    score *= flux_deriv;
    return;
  }
  const Material& material {*model::materials[p->material_]};
  if (material.id_ != deriv.diff_material) {
    score *= flux_deriv;
    return;
  }

  switch (deriv.variable) {

  //============================================================================
  // Density derivative:
  // c = Sigma_MT
  // c = sigma_MT * N
  // c = sigma_MT * rho * const
  // d_c / d_rho = sigma_MT * const
  // (1 / c) * (d_c / d_rho) = 1 / rho

  case DerivativeVariable::DENSITY:
    switch (tally.estimator_) {

    case TallyEstimator::ANALOG:
    case TallyEstimator::COLLISION:
      switch (score_bin) {

      case SCORE_TOTAL:
      case SCORE_SCATTER:
      case SCORE_ABSORPTION:
      case SCORE_FISSION:
      case SCORE_NU_FISSION:
        score *= flux_deriv + 1. / material.density_gpcc_;
        break;

      default:
        fatal_error("Tally derivative not defined for a score on tally "
          + std::to_string(tally.id_));
      }
      break;

    default:
      fatal_error("Differential tallies are only implemented for analog and "
        "collision estimators.");
    }
    break;

  //============================================================================
  // Nuclide density derivative:
  // If we are scoring a reaction rate for a single nuclide then
  // c = Sigma_MT_i
  // c = sigma_MT_i * N_i
  // d_c / d_N_i = sigma_MT_i
  // (1 / c) * (d_c / d_N_i) = 1 / N_i
  // If the score is for the total material (i_nuclide = -1)
  // c = Sum_i(Sigma_MT_i)
  // d_c / d_N_i = sigma_MT_i
  // (1 / c) * (d_c / d_N) = sigma_MT_i / Sigma_MT
  // where i is the perturbed nuclide.

  case DerivativeVariable::NUCLIDE_DENSITY:
    switch (tally.estimator_) {

    case TallyEstimator::ANALOG:
      if (p->event_nuclide_ != deriv.diff_nuclide) {
        score *= flux_deriv;
        return;
      }

      switch (score_bin) {

      case SCORE_TOTAL:
      case SCORE_SCATTER:
      case SCORE_ABSORPTION:
      case SCORE_FISSION:
      case SCORE_NU_FISSION:
        {
          // Find the index of the perturbed nuclide.
          int i;
          for (i = 0; i < material.nuclide_.size(); ++i)
            if (material.nuclide_[i] == deriv.diff_nuclide) break;
          score *= flux_deriv + 1. / material.atom_density_(i);
        }
        break;

      default:
        fatal_error("Tally derivative not defined for a score on tally "
          + std::to_string(tally.id_));
      }
      break;

    case TallyEstimator::COLLISION:
      switch (score_bin) {

      case SCORE_TOTAL:
        if (i_nuclide == -1 && p->macro_xs_.total > 0.0) {
          score *= flux_deriv
            + p->neutron_xs_[deriv.diff_nuclide].total
            / p->macro_xs_.total;
        } else if (i_nuclide == deriv.diff_nuclide
                   && p->neutron_xs_[i_nuclide].total) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_SCATTER:
        if (i_nuclide == -1 && (p->macro_xs_.total
                                - p->macro_xs_.absorption) > 0.0) {
          score *= flux_deriv
            + (p->neutron_xs_[deriv.diff_nuclide].total
            - p->neutron_xs_[deriv.diff_nuclide].absorption)
            / (p->macro_xs_.total
            - p->macro_xs_.absorption);
        } else if (i_nuclide == deriv.diff_nuclide) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_ABSORPTION:
        if (i_nuclide == -1 && p->macro_xs_.absorption > 0.0) {
          score *= flux_deriv
            + p->neutron_xs_[deriv.diff_nuclide].absorption
            / p->macro_xs_.absorption;
        } else if (i_nuclide == deriv.diff_nuclide
                   && p->neutron_xs_[i_nuclide].absorption) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_FISSION:
        if (i_nuclide == -1 && p->macro_xs_.fission > 0.0) {
          score *= flux_deriv
            + p->neutron_xs_[deriv.diff_nuclide].fission
            / p->macro_xs_.fission;
        } else if (i_nuclide == deriv.diff_nuclide
                   && p->neutron_xs_[i_nuclide].fission) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_NU_FISSION:
        if (i_nuclide == -1 && p->macro_xs_.nu_fission > 0.0) {
          score *= flux_deriv
            + p->neutron_xs_[deriv.diff_nuclide].nu_fission
            / p->macro_xs_.nu_fission;
        } else if (i_nuclide == deriv.diff_nuclide
                   && p->neutron_xs_[i_nuclide].nu_fission) {
          score *= flux_deriv + 1. / atom_density;
        } else {
          score *= flux_deriv;
        }
        break;

      default:
        fatal_error("Tally derivative not defined for a score on tally "
          + std::to_string(tally.id_));
      }
      break;

    default:
      fatal_error("Differential tallies are only implemented for analog and "
        "collision estimators.");
    }
    break;

  //============================================================================
  // Temperature derivative:
  // If we are scoring a reaction rate for a single nuclide then
  // c = Sigma_MT_i
  // c = sigma_MT_i * N_i
  // d_c / d_T = (d_sigma_Mt_i / d_T) * N_i
  // (1 / c) * (d_c / d_T) = (d_sigma_MT_i / d_T) / sigma_MT_i
  // If the score is for the total material (i_nuclide = -1)
  // (1 / c) * (d_c / d_T) = Sum_i((d_sigma_MT_i / d_T) * N_i) / Sigma_MT_i
  // where i is the perturbed nuclide.  The d_sigma_MT_i / d_T term is
  // computed by multipole_deriv_eval.  It only works for the resolved
  // resonance range and requires multipole data.

  case DerivativeVariable::TEMPERATURE:
    switch (tally.estimator_) {

    case TallyEstimator::ANALOG:
      {
        // Find the index of the event nuclide.
        int i;
        for (i = 0; i < material.nuclide_.size(); ++i)
          if (material.nuclide_[i] == p->event_nuclide_) break;

        const auto& nuc {*data::nuclides[p->event_nuclide_]};
        if (!multipole_in_range(&nuc, p->E_last_)) {
          score *= flux_deriv;
          break;
        }

        switch (score_bin) {

        case SCORE_TOTAL:
          if (p->neutron_xs_[p->event_nuclide_].total) {
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + (dsig_s + dsig_a) * material.atom_density_(i)
              / p->macro_xs_.total;
          } else {
            score *= flux_deriv;
          }
          break;

        case SCORE_SCATTER:
          if (p->neutron_xs_[p->event_nuclide_].total
              - p->neutron_xs_[p->event_nuclide_].absorption) {
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + dsig_s * material.atom_density_(i)
              / (p->macro_xs_.total
              - p->macro_xs_.absorption);
          } else {
            score *= flux_deriv;
          }
          break;

        case SCORE_ABSORPTION:
          if (p->neutron_xs_[p->event_nuclide_].absorption) {
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + dsig_a * material.atom_density_(i)
              / p->macro_xs_.absorption;
          } else {
            score *= flux_deriv;
          }
          break;

        case SCORE_FISSION:
          if (p->neutron_xs_[p->event_nuclide_].fission) {
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + dsig_f * material.atom_density_(i)
              / p->macro_xs_.fission;
          } else {
            score *= flux_deriv;
          }
          break;

        case SCORE_NU_FISSION:
          if (p->neutron_xs_[p->event_nuclide_].fission) {
            double nu = p->neutron_xs_[p->event_nuclide_].nu_fission
              / p->neutron_xs_[p->event_nuclide_].fission;
            double dsig_s, dsig_a, dsig_f;
            std::tie(dsig_s, dsig_a, dsig_f)
              = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
            score *= flux_deriv + nu * dsig_f * material.atom_density_(i)
              / p->macro_xs_.nu_fission;
          } else {
            score *= flux_deriv;
          }
          break;

        default:
          fatal_error("Tally derivative not defined for a score on tally "
            + std::to_string(tally.id_));
        }
      }
      break;

    case TallyEstimator::COLLISION:
      if (i_nuclide != -1) {
        const auto& nuc {data::nuclides[i_nuclide]};
        if (!multipole_in_range(nuc.get(), p->E_last_)) {
          score *= flux_deriv;
          return;
        }
      }

      switch (score_bin) {

      case SCORE_TOTAL:
        if (i_nuclide == -1 && p->macro_xs_.total > 0.0) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && p->neutron_xs_[i_nuc].total) {
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += (dsig_s + dsig_a) * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / p->macro_xs_.total;
        } else if (p->neutron_xs_[i_nuclide].total) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv
            + (dsig_s + dsig_a) / p->neutron_xs_[i_nuclide].total;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_SCATTER:
        if (i_nuclide == -1 && (p->macro_xs_.total
            - p->macro_xs_.absorption)) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && (p->neutron_xs_[i_nuc].total
                - p->neutron_xs_[i_nuc].absorption)) {
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += dsig_s * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / (p->macro_xs_.total
            - p->macro_xs_.absorption);
        } else if (p->neutron_xs_[i_nuclide].total
                   - p->neutron_xs_[i_nuclide].absorption) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv + dsig_s / (p->neutron_xs_[i_nuclide].total
            - p->neutron_xs_[i_nuclide].absorption);
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_ABSORPTION:
        if (i_nuclide == -1 && p->macro_xs_.absorption > 0.0) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && p->neutron_xs_[i_nuc].absorption) {
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += dsig_a * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / p->macro_xs_.absorption;
        } else if (p->neutron_xs_[i_nuclide].absorption) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv
            + dsig_a / p->neutron_xs_[i_nuclide].absorption;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_FISSION:
        if (i_nuclide == -1 && p->macro_xs_.fission > 0.0) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && p->neutron_xs_[i_nuc].fission) {
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += dsig_f * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / p->macro_xs_.fission;
        } else if (p->neutron_xs_[i_nuclide].fission) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv
            + dsig_f / p->neutron_xs_[i_nuclide].fission;
        } else {
          score *= flux_deriv;
        }
        break;

      case SCORE_NU_FISSION:
        if (i_nuclide == -1 && p->macro_xs_.nu_fission > 0.0) {
          double cum_dsig = 0;
          for (auto i = 0; i < material.nuclide_.size(); ++i) {
            auto i_nuc = material.nuclide_[i];
            const auto& nuc {*data::nuclides[i_nuc]};
            if (multipole_in_range(&nuc, p->E_last_)
                && p->neutron_xs_[i_nuc].fission) {
              double nu = p->neutron_xs_[i_nuc].nu_fission
                / p->neutron_xs_[i_nuc].fission;
              double dsig_s, dsig_a, dsig_f;
              std::tie(dsig_s, dsig_a, dsig_f)
                = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
              cum_dsig += nu * dsig_f * material.atom_density_(i);
            }
          }
          score *= flux_deriv + cum_dsig / p->macro_xs_.nu_fission;
        } else if (p->neutron_xs_[i_nuclide].fission) {
          const auto& nuc {*data::nuclides[i_nuclide]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          score *= flux_deriv
            + dsig_f / p->neutron_xs_[i_nuclide].fission;
        } else {
          score *= flux_deriv;
        }
        break;

      default:
        break;
      }
      break;

    default:
      fatal_error("Differential tallies are only implemented for analog and "
        "collision estimators.");
    }
    break;
  }
}

void
score_track_sensitivity(Particle* p, double distance)
{
  // A void material cannot be perturbed so it will not affect flux sensitivities.
  if (p->material_ == MATERIAL_VOID) return;
  const Material& material {*model::materials[p->material_]};

  for (auto idx = 0; idx < model::tally_sens.size(); idx++) {
    const auto& sens = model::tally_sens[idx];
    auto& cumulative_sensitivities = p->cumulative_sensitivities_[idx];
    if (sens.sens_material != material.id_) continue;

    switch (sens.variable) {

    case SensitivityVariable::CROSS_SECTION:
      // Calculate the sensitivity with respect to the cross section
      // at this energy

      // Get the pre-collision energy of the particle.
      auto E = p->E_last_;
  
      // Bin the energy.
      if (E >= sens.energy_bins_.front() && E <= sens.energy_bins_.back()) {
        auto bin = lower_bound_index(sens.energy_bins_.begin(), sens.energy_bins_.end(), E);
      }

      double atom_density = 0.;
          if (sens.sens_nuclide >= 0) {
            auto j = model::materials[p->material_]->mat_nuclide_index_[sens.sens_nuclide];
            if (j == C_NONE) continue;
            atom_density = model::materials[p->material_]->atom_density_(j);
          }


      // Get the correct cross section
      double macro_xs;
      switch (sens.sens_reaction) {
      case SCORE_TOTAL:
        if (sens.sens_nuclide >=0){
            macro_xs = p->neutron_xs_[sens.sens_nuclide].total * atom_density;
        } else {
            macro_xs = p->macro_xs_.total;
        }
        break;
      case SCORE_SCATTER:
        if (sens.sens_nuclide >=0){
            macro_xs = (p->neutron_xs_[sens.sens_nuclide].total 
            - p->neutron_xs_[sens.sens_nuclide].absorption) * atom_density;
        } else {
            macro_xs = p->macro_xs_.total - p->macro_xs_.absorption;
        }
        break;
      case SCORE_ABSORPTION: 
        if (sens.sens_nuclide >=0){
            macro_xs = p->neutron_xs_[sens.sens_nuclide].absorption * atom_density;
        } else {
            macro_xs = p->macro_xs_.absorption;
        }
        break;
      case SCORE_FISSION:
        if (sens.sens_nuclide >=0){
            macro_xs = (p->neutron_xs_[sens.sens_nuclide].total 
            - p->neutron_xs_[sens.sens_nuclide].absorption) * atom_density;
        } else {
            macro_xs = p->macro_xs_.total - p->macro_xs_.absorption;
        }
        break;
      }


      
      cumulative_sensitivities[bin] -= distance * macro_xs;


      break;

    case SensitivityVariable::MULTIPOLE:
      for (auto i = 0; i < material.nuclide_.size(); ++i) {
        const auto& nuc {*data::nuclides[material.nuclide_[i]]};
        if (multipole_in_range(&nuc, p->E_last_)) {
          // phi is proportional to e^(-Sigma_tot * dist)
          // (1 / phi) * (d_phi / d_T) = - (d_Sigma_tot / d_T) * dist
          // (1 / phi) * (d_phi / d_T) = - N (d_sigma_tot / d_T) * dist
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_, p->sqrtkT_);
          flux_deriv -= distance * (dsig_s + dsig_a)
            * material.atom_density_(i);
        }
      }
      break;
    }
  }
}

void score_collision_sensitivity(Particle* p)
{
  // A void material cannot be perturbed so it will not affect flux derivatives.
  if (p->material_ == MATERIAL_VOID) return;

  const Material& material {*model::materials[p->material_]};

  for (auto idx = 0; idx < model::tally_sens.size(); idx++) {
    const auto& sens = model::tally_sens[idx];
    auto& cumulative_sensitivities = p->cumulative_sensitivities_[idx];

    if (sens.sens_material != material.id_) continue;

    switch (sens.variable) {

    case SensitivityVariable::CROSS_SECTION:
      if (p->event_nuclide_ != deriv.diff_nuclide) continue;
      // Find the index in this material for the diff_nuclide.
      int i;
      for (i = 0; i < material.nuclide_.size(); ++i)
        if (material.nuclide_[i] == deriv.diff_nuclide) break;
      // Make sure we found the nuclide.
      if (material.nuclide_[i] != deriv.diff_nuclide) {
        fatal_error(fmt::format(
          "Could not find nuclide {} in material {} for tally derivative {}",
          data::nuclides[deriv.diff_nuclide]->name_, material.id_, deriv.id));
      }
      // phi is proportional to Sigma_s
      // (1 / phi) * (d_phi / d_N) = (d_Sigma_s / d_N) / Sigma_s
      // (1 / phi) * (d_phi / d_N) = sigma_s / Sigma_s
      // (1 / phi) * (d_phi / d_N) = 1 / N
      flux_deriv += 1. / material.atom_density_(i);
      break;

    case SensitivityVariable::MULTIPOLE:
      // Loop over the material's nuclides until we find the event nuclide.
      for (auto i_nuc : material.nuclide_) {
        const auto& nuc {*data::nuclides[i_nuc]};
        if (i_nuc == p->event_nuclide_ && multipole_in_range(&nuc, p->E_last_)) {
          // phi is proportional to Sigma_s
          // (1 / phi) * (d_phi / d_T) = (d_Sigma_s / d_T) / Sigma_s
          // (1 / phi) * (d_phi / d_T) = (d_sigma_s / d_T) / sigma_s
          const auto& micro_xs {p->neutron_xs_[i_nuc]};
          double dsig_s, dsig_a, dsig_f;
          std::tie(dsig_s, dsig_a, dsig_f)
            = nuc.multipole_->evaluate_deriv(p->E_last_, p->sqrtkT_);
          flux_deriv += dsig_s / (micro_xs.total - micro_xs.absorption);
          // Note that this is an approximation!  The real scattering cross
          // section is
          // Sigma_s(E'->E, u'->u) = Sigma_s(E') * P(E'->E, u'->u).
          // We are assuming that d_P(E'->E, u'->u) / d_T = 0 and only
          // computing d_S(E') / d_T.  Using this approximation in the vicinity
          // of low-energy resonances causes errors (~2-5% for PWR pincell
          // eigenvalue derivatives).
        }
      }
      break;
    }
  }
}

}// namespace openmc
