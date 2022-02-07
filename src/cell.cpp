
#include "openmc/cell.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iterator>
#include <set>
#include <sstream>
#include <string>

#include <fmt/core.h>
#include <gsl/gsl>

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/dagmc.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/nuclide.h"
#include "openmc/settings.h"
#include "openmc/surface.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  std::unordered_map<int32_t, int32_t> cell_map;
  vector<unique_ptr<CSGCell>> cells;

  std::unordered_map<int32_t, int32_t> universe_map;
  vector<unique_ptr<Universe>> universes;
} // namespace model

#ifdef __CUDACC__
namespace gpu {
__constant__ unique_ptr<CSGCell>* cells;
__constant__ unique_ptr<Universe>* universes;
} // namespace gpu
#endif

//==============================================================================
// Universe implementation
//==============================================================================

void
Universe::to_hdf5(hid_t universes_group) const
{
  // Create a group for this universe.
  auto group = create_group(universes_group, fmt::format("universe {}", id_));

  // Write the contained cells.
  if (cells_.size() > 0) {
    vector<int32_t> cell_ids;
    for (auto i_cell : cells_) cell_ids.push_back(model::cells[i_cell]->id_);
    write_dataset(group, "cells", cell_ids);
  }

  close_group(group);
}

BoundingBox Universe::bounding_box() const {
  BoundingBox bbox = {INFTY, -INFTY, INFTY, -INFTY, INFTY, -INFTY};
  if (cells_.size() == 0) {
    return {};
  } else {
    for (const auto& cell : cells_) {
      auto& c = model::cells[cell];
      bbox |= c->bounding_box();
    }
  }
  return bbox;
}

//==============================================================================
// Cell implementation
//==============================================================================

double HD Cell::temperature(int32_t instance) const
{
  if (sqrtkT_.size() < 1) {
#ifdef __CUDA_ARCH__
    asm("trap;");
#else
    throw std::runtime_error{"Cell temperature has not yet been set."};
#endif
  }

  if (instance >= 0) {
    double sqrtkT = sqrtkT_.size() == 1 ?
      sqrtkT_.at(0) :
      sqrtkT_.at(instance);
    return sqrtkT * sqrtkT / K_BOLTZMANN;
  } else {
    return sqrtkT_[0] * sqrtkT_[0] / K_BOLTZMANN;
  }
}

void
Cell::set_temperature(double T, int32_t instance, bool set_contained)
{
  if (settings::temperature_method == TemperatureMethod::INTERPOLATION) {
    if (T < data::temperature_min) {
      throw std::runtime_error{"Temperature is below minimum temperature at "
        "which data is available."};
    } else if (T > data::temperature_max) {
      throw std::runtime_error{"Temperature is above maximum temperature at "
        "which data is available."};
    }
  }

  if (type_ == Fill::MATERIAL) {
    if (instance >= 0) {
      // If temperature vector is not big enough, resize it first
      if (sqrtkT_.size() != n_instances_) sqrtkT_.resize(n_instances_, sqrtkT_[0]);

      // Set temperature for the corresponding instance
      sqrtkT_.at(instance) = std::sqrt(K_BOLTZMANN * T);
    } else {
      // Set temperature for all instances
      for (auto& T_ : sqrtkT_) {
        T_ = std::sqrt(K_BOLTZMANN * T);
      }
    }
  } else {
    if (!set_contained) {
      throw std::runtime_error{fmt::format("Attempted to set the temperature of cell {} "
                                           "which is not filled by a material.", id_)};
    }

    auto contained_cells = this->get_contained_cells();
    for (const auto& entry : contained_cells) {
      auto& cell = model::cells[entry.first];
      Expects(cell->type_ == Fill::MATERIAL);
      auto& instances =  entry.second;
      for (auto instance : instances) {
        cell->set_temperature(T, instance);
      }
    }
  }
}

//==============================================================================
// CSGCell implementation
//==============================================================================

CSGCell::CSGCell(pugi::xml_node cell_node)
{
  if (check_for_node(cell_node, "id")) {
    id_ = std::stoi(get_node_value(cell_node, "id"));
  } else {
    fatal_error("Must specify id of cell in geometry XML file.");
  }

  if (check_for_node(cell_node, "name")) {
    name_ = get_node_value(cell_node, "name");
  }

  if (check_for_node(cell_node, "universe")) {
    universe_ = std::stoi(get_node_value(cell_node, "universe"));
  } else {
    universe_ = 0;
  }

  // Make sure that either material or fill was specified, but not both.
  bool fill_present = check_for_node(cell_node, "fill");
  bool material_present = check_for_node(cell_node, "material");
  if (!(fill_present || material_present)) {
    fatal_error(fmt::format(
      "Neither material nor fill was specified for cell {}", id_));
  }
  if (fill_present && material_present) {
    fatal_error(fmt::format("Cell {} has both a material and a fill specified; "
      "only one can be specified per cell", id_));
  }

  if (fill_present) {
    fill_ = std::stoi(get_node_value(cell_node, "fill"));
    if (fill_ == universe_) {
      fatal_error(fmt::format("Cell {} is filled with the same universe that"
        "it is contained in.", id_));
    }
  } else {
    fill_ = C_NONE;
  }

  // Read the material element.  There can be zero materials (filled with a
  // universe), more than one material (distribmats), and some materials may
  // be "void".
  if (material_present) {
    vector<string> mats {get_node_array<string>(cell_node, "material", true)};
    if (mats.size() > 0) {
      material_.reserve(mats.size());
      for (std::string mat : mats) {
        if (mat.compare("void") == 0) {
          // THis is to avoid a branch in the XS lookup kernel. Because most of the time,
          // voids are indeed filled with gas, the solution here is to just define a material
          // with no nuclides in it.
          fatal_error("Void materials not treated correctly in GPU mode!");
          // material_.push_back(MATERIAL_VOID);
        } else {
          material_.push_back(std::stoi(mat));
        }
      }
    } else {
      fatal_error(fmt::format("An empty material element was specified for cell {}",
        id_));
    }
  }

  // Read the temperature element which may be distributed like materials.
  if (check_for_node(cell_node, "temperature")) {
    sqrtkT_ = get_node_array<double>(cell_node, "temperature");
    sqrtkT_.shrink_to_fit();

    // Make sure this is a material-filled cell.
    if (material_.size() == 0) {
      fatal_error(fmt::format(
        "Cell {} was specified with a temperature but no material. Temperature"
        "specification is only valid for cells filled with a material.", id_));
    }

    // Make sure all temperatures are non-negative.
    for (auto T : sqrtkT_) {
      if (T < 0) {
        fatal_error(fmt::format(
          "Cell {} was specified with a negative temperature", id_));
      }
    }

    // Convert to sqrt(k*T).
    for (auto& T : sqrtkT_) {
      T = std::sqrt(K_BOLTZMANN * T);
    }
  }

  // Read the region specification.
  std::string region_spec;
  if (check_for_node(cell_node, "region")) {
    region_spec = get_node_value(cell_node, "region");
  }

  // Get a tokenized representation of the region specification.
  // Should throw if user attempts to use "complex" cell.
  std::stringstream firstpass(region_spec);
  std::string value;
  int region_size = 0;
  while (firstpass >> value) region_size++;
  std::stringstream secondpass(region_spec);
  region_.reserve(region_size);
  while (secondpass >> value) region_.push_back(std::stoi(value));

  // Convert user IDs to surface indices.
  for (auto& r : region_) {
    if (r < OP_UNION) {
      const auto& it {model::surface_map.find(abs(r))};
      if (it == model::surface_map.end()) {
        throw std::runtime_error{"Invalid surface ID " + std::to_string(abs(r))
          + " specified in region for cell " + std::to_string(id_) + "."};
      }
      r = (r > 0) ? it->second + 1 : -(it->second + 1);
    }
  }

  // Read the translation vector.
  if (check_for_node(cell_node, "translation")) {
    if (fill_ == C_NONE) {
      fatal_error(fmt::format("Cannot apply a translation to cell {}"
        " because it is not filled with another universe", id_));
    }

    auto xyz {get_node_array<double>(cell_node, "translation")};
    if (xyz.size() != 3) {
      fatal_error(fmt::format(
        "Non-3D translation vector applied to cell {}", id_));
    }
    translation_ = xyz;
  }

  // Read the rotation transform.
  if (check_for_node(cell_node, "rotation")) {
    if (fill_ == C_NONE) {
      fatal_error(fmt::format("Cannot apply a rotation to cell {}"
        " because it is not filled with another universe", id_));
    }

    auto rot {get_node_array<double>(cell_node, "rotation")};
    if (rot.size() != 3 && rot.size() != 9) {
      fatal_error(fmt::format(
        "Non-3D rotation vector applied to cell {}", id_));
    }

    // Compute and store the rotation matrix.
    rotation_.reserve(rot.size() == 9 ? 9 : 12);
    if (rot.size() == 3) {
      double phi = -rot[0] * PI / 180.0;
      double theta = -rot[1] * PI / 180.0;
      double psi = -rot[2] * PI / 180.0;
      rotation_.push_back(std::cos(theta) * std::cos(psi));
      rotation_.push_back(-std::cos(phi) * std::sin(psi)
                          + std::sin(phi) * std::sin(theta) * std::cos(psi));
      rotation_.push_back(std::sin(phi) * std::sin(psi)
                          + std::cos(phi) * std::sin(theta) * std::cos(psi));
      rotation_.push_back(std::cos(theta) * std::sin(psi));
      rotation_.push_back(std::cos(phi) * std::cos(psi)
                          + std::sin(phi) * std::sin(theta) * std::sin(psi));
      rotation_.push_back(-std::sin(phi) * std::cos(psi)
                          + std::cos(phi) * std::sin(theta) * std::sin(psi));
      rotation_.push_back(-std::sin(theta));
      rotation_.push_back(std::sin(phi) * std::cos(theta));
      rotation_.push_back(std::cos(phi) * std::cos(theta));

      // When user specifies angles, write them at end of vector
      rotation_.push_back(rot[0]);
      rotation_.push_back(rot[1]);
      rotation_.push_back(rot[2]);
    } else {
      rotation_.reserve(rot.size());
      std::copy(rot.begin(), rot.end(), std::back_inserter(rotation_));
    }
  }
}

//==============================================================================

void
CSGCell::to_hdf5(hid_t cell_group) const
{
  // Create a group for this cell.
  auto group = create_group(cell_group, fmt::format("cell {}", id_));

  if (!name_.empty()) {
    write_string(group, "name", name_, false);
  }

  write_dataset(group, "universe", model::universes[universe_]->id_);

  // Write the region specification.
  if (!region_.empty()) {
    std::stringstream region_spec {};
    for (int32_t token : region_) {
      if (token == OP_LEFT_PAREN) {
        region_spec << " (";
      } else if (token == OP_RIGHT_PAREN) {
        region_spec << " )";
      } else if (token == OP_COMPLEMENT) {
        region_spec << " ~";
      } else if (token == OP_INTERSECTION) {
      } else if (token == OP_UNION) {
        region_spec << " |";
      } else {
        // Note the off-by-one indexing
        const auto surf_id = model::surfaces[abs(token)-1].id_;
        region_spec << " " << ((token > 0) ? surf_id : -surf_id);
      }
    }
    write_string(group, "region", region_spec.str(), false);
  }

  // Write fill information.
  if (type_ == Fill::MATERIAL) {
    write_dataset(group, "fill_type", "material");
    vector<int32_t> mat_ids;
    for (auto i_mat : material_) {
      if (i_mat != MATERIAL_VOID) {
        mat_ids.push_back(model::materials[i_mat]->id_);
      } else {
        mat_ids.push_back(MATERIAL_VOID);
      }
    }
    if (mat_ids.size() == 1) {
      write_dataset(group, "material", mat_ids[0]);
    } else {
      write_dataset(group, "material", mat_ids);
    }

    vector<double> temps;
    for (auto sqrtkT_val : sqrtkT_)
      temps.push_back(sqrtkT_val * sqrtkT_val / K_BOLTZMANN);
    write_dataset(group, "temperature", temps);

  } else if (type_ == Fill::UNIVERSE) {
    write_dataset(group, "fill_type", "universe");
    write_dataset(group, "fill", model::universes[fill_]->id_);
    if (translation_ != Position(0, 0, 0)) {
      write_dataset(group, "translation", translation_);
    }
    if (!rotation_.empty()) {
      if (rotation_.size() == 12) {
        std::array<double, 3> rot {rotation_[9], rotation_[10], rotation_[11]};
        write_dataset(group, "rotation", rot);
      } else {
        write_dataset(group, "rotation", rotation_);
      }
    }

  } else if (type_ == Fill::LATTICE) {
    write_dataset(group, "fill_type", "lattice");
    write_dataset(group, "lattice", model::lattices[fill_]->id_);
  }

  close_group(group);
}

BoundingBox CSGCell::bounding_box_simple() const {
  BoundingBox bbox;
  for (int32_t token : region_) {
    bbox &= model::surfaces[abs(token)-1].bounding_box(token > 0);
  }
  return bbox;
}

BoundingBox CSGCell::bounding_box() const {
  return bounding_box_simple();
}

//==============================================================================
// UniversePartitioner implementation
//==============================================================================

UniversePartitioner::UniversePartitioner(const Universe& univ)
{
  // Define an ordered set of surface indices that point to z-planes.  Use a
  // functor to to order the set by the z0_ values of the corresponding planes.
  struct compare_surfs {
    bool operator()(const int32_t& i_surf, const int32_t& j_surf) const
    {
      const double zi = model::surfaces[i_surf].storage_.zp.z0_;
      const double zj = model::surfaces[j_surf].storage_.zp.z0_;
      return zi < zj;
    }
  };
  std::set<int32_t, compare_surfs> surf_set;

  // Find all of the z-planes in this universe.  A set is used here for the
  // O(log(n)) insertions that will ensure entries are not repeated.
  for (auto i_cell : univ.cells_) {
    for (auto token : model::cells[i_cell]->region_) {
      auto i_surf = std::abs(token) - 1;
      if (model::surfaces[i_surf].type_ == Surface::SurfaceType::zplane)
        surf_set.insert(i_surf);
    }
  }

  // Populate the surfs_ vector from the ordered set.
  surfs_.insert(surfs_.begin(), surf_set.begin(), surf_set.end());

  // Populate the partition lists.
  partitions_.resize(surfs_.size() + 1);
  for (auto i_cell : univ.cells_) {

    // Find the tokens for bounding z-planes.
    int32_t lower_token = 0, upper_token = 0;
    double min_z, max_z;
    for (auto token : model::cells[i_cell]->region_) {
      if (model::surfaces[std::abs(token)-1].type_ == Surface::SurfaceType::zplane) {
        const auto* zplane = &model::surfaces[std::abs(token)-1].storage_.zp;
        if (lower_token == 0 || zplane->z0_ < min_z) {
          lower_token = token;
          min_z = zplane->z0_;
        }
        if (upper_token == 0 || zplane->z0_ > max_z) {
          upper_token = token;
          max_z = zplane->z0_;
        }
      }
    }

    // If there are no bounding z-planes, add this cell to all partitions.
    if (lower_token == 0) {
      for (auto& p : partitions_) p.push_back(i_cell);
      continue;
    }

    // Find the first partition this cell lies in.  If the lower_token indicates
    // a negative halfspace, then the cell is unbounded in the lower direction
    // and it lies in the first partition onward.  Otherwise, it is bounded by
    // the positive halfspace given by the lower_token.
    int first_partition = 0;
    if (lower_token > 0) {
      for (int i = 0; i < surfs_.size(); ++i) {
        if (lower_token == surfs_[i] + 1) {
          first_partition = i + 1;
          break;
        }
      }
    }

    // Find the last partition this cell lies in.  The logic is analogous to the
    // logic for first_partition.
    int last_partition = surfs_.size();
    if (upper_token < 0) {
      for (int i = first_partition; i < surfs_.size(); ++i) {
        if (upper_token == -(surfs_[i] + 1)) {
          last_partition = i;
          break;
        }
      }
    }

    // Add the cell to all relevant partitions.
    for (int i = first_partition; i <= last_partition; ++i) {
      partitions_[i].push_back(i_cell);
    }
  }
}

HD const vector<int32_t>& UniversePartitioner::get_cells(
  Position r, Direction u) const
{
#ifdef __CUDA_ARCH__
  using gpu::surfaces;
#else
  using model::surfaces;
#endif
  // Perform a binary search for the partition containing the given coordinates.
  int left = 0;
  int middle = (surfs_.size() - 1) / 2;
  int right = surfs_.size() - 1;
  while (true) {
    // Check the sense of the coordinates for the current surface.
    const auto& surf = surfaces[surfs_[middle]];
    if (surf.sense(r, u)) {
      // The coordinates lie in the positive halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the positive
      // side of this surface.
      int right_leaf = right - (right - middle) / 2;
      if (right_leaf != middle) {
        left = middle + 1;
        middle = right_leaf;
      } else {
        return partitions_[middle+1];
      }

    } else {
      // The coordinates lie in the negative halfspace.  Recurse if there are
      // more surfaces to check.  Otherwise, return the cells on the negative
      // side of this surface.
      int left_leaf = left + (middle - left) / 2;
      if (left_leaf != middle) {
        right = middle-1;
        middle = left_leaf;
      } else {
        return partitions_[middle];
      }
    }
  }
}

//==============================================================================
// Non-method functions
//==============================================================================

void read_cells(pugi::xml_node node)
{
  // Count the number of cells.
  int n_cells = 0;
  for (pugi::xml_node cell_node: node.children("cell")) {n_cells++;}
  if (n_cells == 0) {
    fatal_error("No cells found in geometry.xml!");
  }

  // Loop over XML cell elements and populate the array.
  model::cells.reserve(n_cells);
  for (pugi::xml_node cell_node : node.children("cell")) {
    model::cells.push_back(make_unique<CSGCell>(cell_node));
  }

  // Fill the cell map.
  for (int i = 0; i < model::cells.size(); i++) {
    int32_t id = model::cells[i]->id_;
    auto search = model::cell_map.find(id);
    if (search == model::cell_map.end()) {
      model::cell_map[id] = i;
    } else {
      fatal_error(fmt::format("Two or more cells use the same unique ID: {}", id));
    }
  }

  // Populate the Universe vector and map.
  for (int i = 0; i < model::cells.size(); i++) {
    int32_t uid = model::cells[i]->universe_;
    auto it = model::universe_map.find(uid);
    if (it == model::universe_map.end()) {
      model::universes.push_back(make_unique<Universe>());
      model::universes.back()->id_ = uid;
      model::universes.back()->cells_.push_back(i);
      model::universe_map[uid] = model::universes.size() - 1;
    } else {
      model::universes[it->second]->cells_.push_back(i);
    }
  }
  model::universes.shrink_to_fit();

  // Allocate the cell overlap count if necessary.
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }

#ifdef __CUDACC__
  // Put pointers to start of universe and cell pointer arrays into GPU constant
  // memory
  auto first_cell = model::cells.data();
  auto first_universe = model::universes.data();
  cudaMemcpyToSymbol(gpu::cells, &first_cell, sizeof(unique_ptr<CSGCell>*));
  cudaMemcpyToSymbol(
    gpu::universes, &first_universe, sizeof(unique_ptr<Universe>*));
#endif
}

//==============================================================================
// C-API functions
//==============================================================================

extern "C" int
openmc_cell_get_fill(int32_t index, int* type, int32_t** indices, int32_t* n)
{
  if (index >= 0 && index < model::cells.size()) {
    Cell& c {*model::cells[index]};
    *type = static_cast<int>(c.type_);
    if (c.type_ == Fill::MATERIAL) {
      *indices = c.material_.data();
      *n = c.material_.size();
    } else {
      *indices = &c.fill_;
      *n = 1;
    }
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

extern "C" int
openmc_cell_set_fill(int32_t index, int type, int32_t n,
                     const int32_t* indices)
{
  Fill filltype = static_cast<Fill>(type);
  if (index >= 0 && index < model::cells.size()) {
    Cell& c {*model::cells[index]};
    if (filltype == Fill::MATERIAL) {
      c.type_ = Fill::MATERIAL;
      c.material_.clear();
      for (int i = 0; i < n; i++) {
        int i_mat = indices[i];
        if (i_mat == MATERIAL_VOID) {
          c.material_.push_back(MATERIAL_VOID);
        } else if (i_mat >= 0 && i_mat < model::materials.size()) {
          c.material_.push_back(i_mat);
        } else {
          set_errmsg("Index in materials array is out of bounds.");
          return OPENMC_E_OUT_OF_BOUNDS;
        }
      }
      c.material_.shrink_to_fit();
    } else if (filltype == Fill::UNIVERSE) {
      c.type_ = Fill::UNIVERSE;
    } else {
      c.type_ = Fill::LATTICE;
    }
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
  return 0;
}

extern "C" int
openmc_cell_set_temperature(int32_t index, double T, const int32_t* instance, bool set_contained)
{
  if (index < 0 || index >= model::cells.size()) {
    strcpy(openmc_err_msg, "Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  int32_t instance_index = instance ? *instance : -1;
  try {
    model::cells[index]->set_temperature(T, instance_index, set_contained);
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_UNASSIGNED;
  }
  return 0;
}

extern "C" int
openmc_cell_get_temperature(int32_t index, const int32_t* instance, double* T)
{
  if (index < 0 || index >= model::cells.size()) {
    strcpy(openmc_err_msg, "Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  int32_t instance_index = instance ? *instance : -1;
  try {
    *T = model::cells[index]->temperature(instance_index);
  } catch (const std::exception& e) {
    set_errmsg(e.what());
    return OPENMC_E_UNASSIGNED;
  }
  return 0;
}

//! Get the bounding box of a cell
extern "C" int
openmc_cell_bounding_box(const int32_t index, double* llc, double* urc) {

  BoundingBox bbox;

  const auto& c = model::cells[index];
  bbox = c->bounding_box();

  // set lower left corner values
  llc[0] = bbox.xmin;
  llc[1] = bbox.ymin;
  llc[2] = bbox.zmin;

  // set upper right corner values
  urc[0] = bbox.xmax;
  urc[1] = bbox.ymax;
  urc[2] = bbox.zmax;

  return 0;
}

//! Get the name of a cell
extern "C" int
openmc_cell_get_name(int32_t index, const char** name) {
  if (index < 0 || index >= model::cells.size()) {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  *name = model::cells[index]->name().data();

  return 0;
}

//! Set the name of a cell
extern "C" int
openmc_cell_set_name(int32_t index, const char* name) {
  if (index < 0 || index >= model::cells.size()) {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }

  model::cells[index]->set_name(name);

  return 0;
}

std::unordered_map<int32_t, vector<int32_t>> Cell::get_contained_cells() const
{
  std::unordered_map<int32_t, vector<int32_t>> contained_cells;
  vector<ParentCell> parent_cells;

  // if this cell is filled w/ a material, it contains no other cells
  if (type_ != Fill::MATERIAL) {
    this->get_contained_cells_inner(contained_cells, parent_cells);
  }

  return contained_cells;
}

//! Get all cells within this cell
void Cell::get_contained_cells_inner(
  std::unordered_map<int32_t, vector<int32_t>>& contained_cells,
  vector<ParentCell>& parent_cells) const
{

  // filled by material, determine instance based on parent cells
  if (type_ == Fill::MATERIAL) {
    int instance = 0;
    if (this->distribcell_index_ >= 0) {
      for (auto& parent_cell : parent_cells) {
        auto& cell = openmc::model::cells[parent_cell.cell_index];
        if (cell->type_ == Fill::UNIVERSE) {
          instance += cell->offset_[distribcell_index_];
        } else if (cell->type_ == Fill::LATTICE) {
          auto& lattice = model::lattices[cell->fill_];
          instance += lattice->offset(this->distribcell_index_, parent_cell.lattice_index);
        }
      }
    }
    // add entry to contained cells
    contained_cells[model::cell_map[id_]].push_back(instance);
  // filled with universe, add the containing cell to the parent cells
  // and recurse
  } else if (type_ == Fill::UNIVERSE) {
    parent_cells.push_back({model::cell_map[id_], -1});
    auto& univ = model::universes[fill_];
    for(auto cell_index : univ->cells_) {
      auto& cell = model::cells[cell_index];
      cell->get_contained_cells_inner(contained_cells, parent_cells);
    }
    parent_cells.pop_back();
  // filled with a lattice, visit each universe in the lattice
  // with a recursive call to collect the cell instances
  } else if (type_ == Fill::LATTICE) {
    auto& lattice = model::lattices[fill_];
    for (auto i = lattice->begin(); i != lattice->end(); ++i) {
      auto& univ = model::universes[*i];
      parent_cells.push_back({model::cell_map[id_], i.indx_});
      for (auto cell_index : univ->cells_) {
        auto& cell = model::cells[cell_index];
        cell->get_contained_cells_inner(contained_cells, parent_cells);
      }
      parent_cells.pop_back();
    }
  }
}

//! Return the index in the cells array of a cell with a given ID
extern "C" int
openmc_get_cell_index(int32_t id, int32_t* index)
{
  auto it = model::cell_map.find(id);
  if (it != model::cell_map.end()) {
    *index = it->second;
    return 0;
  } else {
    set_errmsg("No cell exists with ID=" + std::to_string(id) + ".");
    return OPENMC_E_INVALID_ID;
  }
}

//! Return the ID of a cell
extern "C" int
openmc_cell_get_id(int32_t index, int32_t* id)
{
  if (index >= 0 && index < model::cells.size()) {
    *id = model::cells[index]->id_;
    return 0;
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//! Set the ID of a cell
extern "C" int
openmc_cell_set_id(int32_t index, int32_t id)
{
  if (index >= 0 && index < model::cells.size()) {
    model::cells[index]->id_ = id;
    model::cell_map[id] = index;
    return 0;
  } else {
    set_errmsg("Index in cells array is out of bounds.");
    return OPENMC_E_OUT_OF_BOUNDS;
  }
}

//! Extend the cells array by n elements
extern "C" int
openmc_extend_cells(int32_t n, int32_t* index_start, int32_t* index_end)
{
  if (index_start) *index_start = model::cells.size();
  if (index_end) *index_end = model::cells.size() + n - 1;
  for (int32_t i = 0; i < n; i++) {
    model::cells.push_back(make_unique<CSGCell>());
  }
  return 0;
}

#ifdef DAGMC
int32_t next_cell(DAGCell* cur_cell, DAGSurface* surf_xed)
{
  moab::EntityHandle surf =
    surf_xed->dagmc_ptr_->entity_by_index(2, surf_xed->dag_index_);
  moab::EntityHandle vol =
    cur_cell->dagmc_ptr_->entity_by_index(3, cur_cell->dag_index_);

  moab::EntityHandle new_vol;
  cur_cell->dagmc_ptr_->next_vol(surf, vol, new_vol);

  return cur_cell->dagmc_ptr_->index_by_handle(new_vol);
}
#endif

extern "C" int cells_size() { return model::cells.size(); }

} // namespace openmc
