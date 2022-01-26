#include "openmc/surface.h"

#include <cmath>
#include <utility>
#include <set>

#include <fmt/core.h>
#include <gsl/gsl>

#include "openmc/array.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/hdf5_interface.h"
#include "openmc/math_functions.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/string_utils.h"
#include "openmc/xml_interface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
  std::unordered_map<int, int> surface_map;
  vector<Surface> surfaces;
} // namespace model

#ifdef __CUDACC__
namespace gpu {
// Pointer to start of vector of surface pointers on device
__constant__ Surface* surfaces;
} // namespace gpu
#endif

//==============================================================================
// Helper functions for reading the "coeffs" node of an XML surface element
//==============================================================================

void read_coeffs(pugi::xml_node surf_node, double &c1)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 1) {
    fatal_error(fmt::format("Surface expects 1 coeff but was given {}",
      n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf", &c1);
  if (stat != 1) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, double &c1, double &c2,
                 double &c3)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 3) {
    fatal_error(fmt::format("Surface expects 3 coeffs but was given {}",
      n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf", &c1, &c2, &c3);
  if (stat != 3) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, double &c1, double &c2,
                 double &c3, double &c4)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 4) {
    fatal_error(fmt::format("Surface expects 4 coeffs but was given ",
      n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf", &c1, &c2, &c3, &c4);
  if (stat != 4) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

void read_coeffs(pugi::xml_node surf_node, double &c1, double &c2,
                 double &c3, double &c4, double &c5, double &c6, double &c7,
                 double &c8, double &c9, double &c10)
{
  // Check the given number of coefficients.
  std::string coeffs = get_node_value(surf_node, "coeffs");
  int n_words = word_count(coeffs);
  if (n_words != 10) {
    fatal_error(fmt::format("Surface expects 10 coeffs but was given {}",
      n_words));
  }

  // Parse the coefficients.
  int stat = sscanf(coeffs.c_str(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                    &c1, &c2, &c3, &c4, &c5, &c6, &c7, &c8, &c9, &c10);
  if (stat != 10) {
    fatal_error("Something went wrong reading surface coeffs");
  }
}

//==============================================================================
// Surface implementation
//==============================================================================

Surface::Surface(pugi::xml_node surf_node, SurfaceType type)
{

  if (check_for_node(surf_node, "id")) {
    id_ = std::stoi(get_node_value(surf_node, "id"));
    if (contains(settings::source_write_surf_id, id_)) {
      surf_source_ = true;
    }
  } else {
    fatal_error("Must specify id of surface in geometry XML file.");
  }

  if (check_for_node(surf_node, "name")) {
    name_ = get_node_value(surf_node, "name", false);
  }

  if (check_for_node(surf_node, "boundary")) {
    std::string surf_bc = get_node_value(surf_node, "boundary", true, true);

    if (surf_bc == "transmission" || surf_bc == "transmit" ||surf_bc.empty()) {
      // Leave the bc_ a nullptr
    } else if (surf_bc == "vacuum") {
      bc_ = make_unique<VacuumBC>();
    } else if (surf_bc == "reflective" || surf_bc == "reflect"
               || surf_bc == "reflecting") {
      bc_ = make_unique<ReflectiveBC>();
    } else if (surf_bc == "white") {
      bc_ = make_unique<WhiteBC>();
    } else if (surf_bc == "periodic") {
      // periodic BC's are handled separately
    } else {
      fatal_error(fmt::format("Unknown boundary condition \"{}\" specified "
        "on surface {}", surf_bc, id_));
    }
  }

  type_ = type;
  switch (type) {
    case SurfaceType::xplane:
      new (&storage_.xp) SurfaceXPlane(surf_node);
      break;
    case SurfaceType::yplane:
      new (&storage_.yp) SurfaceYPlane(surf_node);
      break;
    case SurfaceType::zplane:
      new (&storage_.zp) SurfaceZPlane(surf_node);
      break;
    case SurfaceType::plane:
      new (&storage_.p) SurfacePlane(surf_node);
      break;
    case SurfaceType::xcylinder:
      new (&storage_.xc) SurfaceXCylinder(surf_node);
      break;
    case SurfaceType::ycylinder:
      new (&storage_.yc) SurfaceYCylinder(surf_node);
      break;
    case SurfaceType::zcylinder:
      new (&storage_.zc) SurfaceZCylinder(surf_node);
      break;
    case SurfaceType::sphere:
      new (&storage_.sph) SurfaceSphere(surf_node);
      break;
    case SurfaceType::xcone:
      new (&storage_.xco) SurfaceXCone(surf_node);
      break;
    case SurfaceType::ycone:
      new (&storage_.yco) SurfaceYCone(surf_node);
      break;
    case SurfaceType::zcone:
      new (&storage_.zco) SurfaceZCone(surf_node);
      break;
    case SurfaceType::quadric:
      new (&storage_.q) SurfaceQuadric(surf_node);
      break;
  }
}

void
Surface::to_hdf5(hid_t group_id) const
{
  std::string group_name {"surface "};
  group_name += std::to_string(id_);

  hid_t surf_group = create_group(group_id, group_name);

  if (bc_) {
    write_string(surf_group, "boundary_type", bc_->type(), false);
  } else {
    write_string(surf_group, "boundary_type", "transmission", false);
  }

  if (!name_.empty()) {
    write_string(surf_group, "name", name_, false);
  }

  to_hdf5_inner(surf_group);

  close_group(surf_group);
}

SurfaceXPlane::SurfaceXPlane(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0_);
}

void SurfaceXPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-plane", false);
  std::array<double, 1> coeffs {{x0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox
SurfaceXPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {x0_, INFTY, -INFTY, INFTY, -INFTY, INFTY};
  } else {
    return {-INFTY, x0_, -INFTY, INFTY, -INFTY, INFTY};
  }
}

//==============================================================================
// SurfaceYPlane implementation
//==============================================================================

SurfaceYPlane::SurfaceYPlane(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, y0_);
}

void SurfaceYPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-plane", false);
  std::array<double, 1> coeffs {{y0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox
SurfaceYPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {-INFTY, INFTY, y0_, INFTY, -INFTY, INFTY};
  } else {
    return {-INFTY, INFTY, -INFTY, y0_, -INFTY, INFTY};
  }
}

//==============================================================================
// SurfaceZPlane implementation
//==============================================================================

SurfaceZPlane::SurfaceZPlane(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, z0_);
}

void SurfaceZPlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-plane", false);
  std::array<double, 1> coeffs {{z0_}};
  write_dataset(group_id, "coefficients", coeffs);
}

BoundingBox
SurfaceZPlane::bounding_box(bool pos_side) const
{
  if (pos_side) {
    return {-INFTY, INFTY, -INFTY, INFTY, z0_, INFTY};
  } else {
    return {-INFTY, INFTY, -INFTY, INFTY, -INFTY, z0_};
  }
}

//==============================================================================
// SurfacePlane implementation
//==============================================================================

SurfacePlane::SurfacePlane(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, A_, B_, C_, D_);
}


void SurfacePlane::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "plane", false);
  std::array<double, 4> coeffs {{A_, B_, C_, D_}};
  write_dataset(group_id, "coefficients", coeffs);
}


//==============================================================================
// SurfaceXCylinder implementation
//==============================================================================

SurfaceXCylinder::SurfaceXCylinder(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, y0_, z0_, radius_);
}

void SurfaceXCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cylinder", false);
  std::array<double, 3> coeffs {{y0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

HD BoundingBox SurfaceXCylinder::bounding_box(bool pos_side) const
{
  if (!pos_side) {
    return {-INFTY, INFTY, y0_ - radius_, y0_ + radius_, z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}
//==============================================================================
// SurfaceYCylinder implementation
//==============================================================================

SurfaceYCylinder::SurfaceYCylinder(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0_, z0_, radius_);
}

void SurfaceYCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cylinder", false);
  std::array<double, 3> coeffs {{x0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

HD BoundingBox SurfaceYCylinder::bounding_box(bool pos_side) const
{
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, -INFTY, INFTY, z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}

//==============================================================================
// SurfaceZCylinder implementation
//==============================================================================

SurfaceZCylinder::SurfaceZCylinder(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0_, y0_, radius_);
}

void SurfaceZCylinder::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cylinder", false);
  std::array<double, 3> coeffs {{x0_, y0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

HD BoundingBox SurfaceZCylinder::bounding_box(bool pos_side) const
{
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_, y0_ - radius_, y0_ + radius_, -INFTY, INFTY};
  } else {
    return {};
  }
}

//==============================================================================
// SurfaceSphere implementation
//==============================================================================

SurfaceSphere::SurfaceSphere(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0_, y0_, z0_, radius_);
}


HD BoundingBox SurfaceSphere::bounding_box(bool pos_side) const
{
  if (!pos_side) {
    return {x0_ - radius_, x0_ + radius_,
            y0_ - radius_, y0_ + radius_,
            z0_ - radius_, z0_ + radius_};
  } else {
    return {};
  }
}


void SurfaceSphere::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "sphere", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceXCone implementation
//==============================================================================

SurfaceXCone::SurfaceXCone(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0_, y0_, z0_, radius_sq_);
}

void SurfaceXCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "x-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceYCone implementation
//==============================================================================

SurfaceYCone::SurfaceYCone(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0_, y0_, z0_, radius_sq_);
}

void SurfaceYCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "y-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceZCone implementation
//==============================================================================

SurfaceZCone::SurfaceZCone(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, x0_, y0_, z0_, radius_sq_);
}

void SurfaceZCone::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "z-cone", false);
  std::array<double, 4> coeffs {{x0_, y0_, z0_, radius_sq_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================
// SurfaceQuadric implementation
//==============================================================================

SurfaceQuadric::SurfaceQuadric(pugi::xml_node surf_node)
{
  read_coeffs(surf_node, A_, B_, C_, D_, E_, F_, G_, H_, J_, K_);
}

void SurfaceQuadric::to_hdf5_inner(hid_t group_id) const
{
  write_string(group_id, "type", "quadric", false);
  std::array<double, 10> coeffs {{A_, B_, C_, D_, E_, F_, G_, H_, J_, K_}};
  write_dataset(group_id, "coefficients", coeffs);
}

//==============================================================================

void read_surfaces(pugi::xml_node node)
{
  // Count the number of surfaces
  int n_surfaces = 0;
  for (pugi::xml_node surf_node : node.children("surface")) {n_surfaces++;}
  if (n_surfaces == 0) {
    fatal_error("No surfaces found in geometry.xml!");
  }

  // Loop over XML surface elements and populate the array.  Keep track of
  // periodic surfaces.
  model::surfaces.reserve(n_surfaces);
  std::set<std::pair<int, int>> periodic_pairs;
  {
    pugi::xml_node surf_node;
    int i_surf;
    for (surf_node = node.child("surface"), i_surf = 0; surf_node;
         surf_node = surf_node.next_sibling("surface"), i_surf++) {
      std::string surf_type = get_node_value(surf_node, "type", true, true);

      // Allocate and initialize the new surface

      Surface::SurfaceType type;
      if (surf_type == "x-plane") {
        type = Surface::SurfaceType::xplane;
      } else if (surf_type == "y-plane") {
        type = Surface::SurfaceType::yplane;
      } else if (surf_type == "z-plane") {
        type = Surface::SurfaceType::zplane;
      } else if (surf_type == "plane") {
        type = Surface::SurfaceType::plane;
      } else if (surf_type == "x-cylinder") {
        type = Surface::SurfaceType::xcylinder;
      } else if (surf_type == "y-cylinder") {
        type = Surface::SurfaceType::ycylinder;
      } else if (surf_type == "z-cylinder") {
        type = Surface::SurfaceType::zcylinder;
      } else if (surf_type == "sphere") {
        type = Surface::SurfaceType::sphere;
      } else if (surf_type == "x-cone") {
        type = Surface::SurfaceType::xcone;
      } else if (surf_type == "y-cone") {
        type = Surface::SurfaceType::ycone;
      } else if (surf_type == "z-cone") {
        type = Surface::SurfaceType::zcone;
      } else if (surf_type == "quadric") {
        type = Surface::SurfaceType::quadric;
      } else {
        fatal_error(fmt::format("Invalid surface type, \"{}\"", surf_type));
      }
      model::surfaces.emplace_back(surf_node, type);

      // Check for a periodic surface
      if (check_for_node(surf_node, "boundary")) {
        std::string surf_bc = get_node_value(surf_node, "boundary", true, true);
        if (surf_bc == "periodic") {
          if (check_for_node(surf_node, "periodic_surface_id")) {
            int i_periodic = std::stoi(get_node_value(surf_node,
                                                      "periodic_surface_id"));
            int lo_id = std::min(model::surfaces.back().id_, i_periodic);
            int hi_id = std::max(model::surfaces.back().id_, i_periodic);
            periodic_pairs.insert({lo_id, hi_id});
          } else {
            periodic_pairs.insert({model::surfaces.back().id_, -1});
          }
        }
      }
    }
  }

  // Fill the surface map
  for (int i_surf = 0; i_surf < model::surfaces.size(); i_surf++) {
    int id = model::surfaces[i_surf].id_;
    auto in_map = model::surface_map.find(id);
    if (in_map == model::surface_map.end()) {
      model::surface_map[id] = i_surf;
    } else {
      fatal_error(fmt::format(
        "Two or more surfaces use the same unique ID: {}", id));
    }
  }

  // Resolve unpaired periodic surfaces.  A lambda function is used with
  // std::find_if to identify the unpaired surfaces.
  auto is_unresolved_pair =
    [](const std::pair<int, int> p){return p.second == -1;};
  auto first_unresolved = std::find_if(periodic_pairs.begin(),
    periodic_pairs.end(), is_unresolved_pair);
  if (first_unresolved != periodic_pairs.end()) {
    // Found one unpaired surface; search for a second one
    auto next_elem = first_unresolved;
    next_elem++;
    auto second_unresolved = std::find_if(next_elem, periodic_pairs.end(),
      is_unresolved_pair);
    if (second_unresolved == periodic_pairs.end()) {
      fatal_error("Found only one periodic surface without a specified partner."
        " Please specify the partner for each periodic surface.");
    }

    // Make sure there isn't a third unpaired surface
    next_elem = second_unresolved;
    next_elem++;
    auto third_unresolved = std::find_if(next_elem,
      periodic_pairs.end(), is_unresolved_pair);
    if (third_unresolved != periodic_pairs.end()) {
      fatal_error("Found at least three periodic surfaces without a specified "
        "partner. Please specify the partner for each periodic surface.");
    }

    // Add the completed pair and remove the old, unpaired entries
    int lo_id = std::min(first_unresolved->first, second_unresolved->first);
    int hi_id = std::max(first_unresolved->first, second_unresolved->first);
    periodic_pairs.insert({lo_id, hi_id});
    periodic_pairs.erase(first_unresolved);
    periodic_pairs.erase(second_unresolved);
  }

  // Assign the periodic boundary conditions
  for (auto periodic_pair : periodic_pairs) {
    int i_surf = model::surface_map[periodic_pair.first];
    int j_surf = model::surface_map[periodic_pair.second];
    Surface& surf1 {model::surfaces[i_surf]};
    Surface& surf2 {model::surfaces[j_surf]};

    // Compute the dot product of the surface normals
    Direction norm1 = surf1.normal({0, 0, 0});
    Direction norm2 = surf2.normal({0, 0, 0});
    norm1 /= norm1.norm();
    norm2 /= norm2.norm();
    double dot_prod = norm1.dot(norm2);

    // If the dot product is 1 (to within floating point precision) then the
    // planes are parallel which indicates a translational periodic boundary
    // condition.  Otherwise, it is a rotational periodic BC.
    if (std::abs(1.0 - dot_prod) < FP_PRECISION) {
      surf1.bc_ = make_unique<TranslationalPeriodicBC>(i_surf, j_surf);
      surf2.bc_ = make_unique<TranslationalPeriodicBC>(i_surf, j_surf);
    } else {
      surf1.bc_ = make_unique<RotationalPeriodicBC>(i_surf, j_surf);
      surf2.bc_ = make_unique<RotationalPeriodicBC>(i_surf, j_surf);
    }
  }

  // Check to make sure a boundary condition was applied to at least one
  // surface
  bool boundary_exists = false;
  for (const auto& surf : model::surfaces) {
    if (surf.bc_) {
      boundary_exists = true;
      break;
    }
  }
  if (settings::run_mode != RunMode::PLOTTING && !boundary_exists) {
    fatal_error("No boundary conditions were applied to any surfaces!");
  }

#ifdef __CUDACC__
  // Save pointer to vector of surface pointers on GPU, since global variables
  // on device are kept separately
  Surface* first_surface_ptr = model::surfaces.data();
  cudaMemcpyToSymbol(
    gpu::surfaces, &first_surface_ptr, sizeof(Surface*));
#endif
}

void free_memory_surfaces()
{
  model::surfaces.clear();
  model::surface_map.clear();
}

} // namespace openmc
