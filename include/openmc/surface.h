#ifndef OPENMC_SURFACE_H
#define OPENMC_SURFACE_H

#include "openmc/string.h"
#include <limits> // For numeric_limits
#include <unordered_map>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/boundary_condition.h"
#include "openmc/constants.h"
#include "openmc/math_functions.h" // for rotate_angle
#include "openmc/memory.h" // for unique_ptr
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/string.h"
#include "openmc/vector.h"
#include "openmc/string.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

class Surface;

namespace model {
  extern std::unordered_map<int, int> surface_map;
  extern vector<Surface> surfaces;
} // namespace model

#ifdef __CUDACC__
namespace gpu {
// Pointer to start of vector of surface pointers on device
extern __constant__ Surface* surfaces;
} // namespace gpu
#endif

//==============================================================================
//! Coordinates for an axis-aligned cuboid bounding a geometric object.
//==============================================================================

struct BoundingBox
{
  double xmin = -INFTY;
  double xmax = INFTY;
  double ymin = -INFTY;
  double ymax = INFTY;
  double zmin = -INFTY;
  double zmax = INFTY;

  BoundingBox() = default;
  BoundingBox(BoundingBox&&) = default;
  BoundingBox(BoundingBox const&) = default;
  BoundingBox& operator=(BoundingBox const&) = default;

  HD inline BoundingBox operator&(const BoundingBox& other)
  {
    BoundingBox result = *this;
    return result &= other;
  }

  HD inline BoundingBox operator|(const BoundingBox& other)
  {
    BoundingBox result = *this;
    return result |= other;
  }

  // intersect operator
  HD inline BoundingBox& operator&=(const BoundingBox& other)
  {
    xmin = std::max(xmin, other.xmin);
    xmax = std::min(xmax, other.xmax);
    ymin = std::max(ymin, other.ymin);
    ymax = std::min(ymax, other.ymax);
    zmin = std::max(zmin, other.zmin);
    zmax = std::min(zmax, other.zmax);
    return *this;
  }

  // union operator
  HD inline BoundingBox& operator|=(const BoundingBox& other)
  {
    xmin = std::min(xmin, other.xmin);
    xmax = std::max(xmax, other.xmax);
    ymin = std::min(ymin, other.ymin);
    ymax = std::max(ymax, other.ymax);
    zmin = std::min(zmin, other.zmin);
    zmax = std::max(zmax, other.zmax);
    return *this;
  }
};

//==============================================================================
//! A plane perpendicular to the x-axis.
//
//! The plane is described by the equation \f$x - x_0 = 0\f$
//==============================================================================

class SurfaceXPlane
{
public:
  explicit SurfaceXPlane(pugi::xml_node surf_node);
  SurfaceXPlane(SurfaceXPlane&& other) = default;

  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;
  HD BoundingBox bounding_box(bool pos_side) const;

  double x0_;
};

//==============================================================================
//! A plane perpendicular to the y-axis.
//
//! The plane is described by the equation \f$y - y_0 = 0\f$
//==============================================================================

class SurfaceYPlane
{
public:
  explicit SurfaceYPlane(pugi::xml_node surf_node);
  SurfaceYPlane(SurfaceYPlane&& other) = default;

  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;
  HD BoundingBox bounding_box(bool pos_side) const;

  double y0_;
};

//==============================================================================
//! A plane perpendicular to the z-axis.
//
//! The plane is described by the equation \f$z - z_0 = 0\f$
//==============================================================================

class SurfaceZPlane
{
public:
  explicit SurfaceZPlane(pugi::xml_node surf_node);
  SurfaceZPlane(SurfaceZPlane&& other) = default;
  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;
  HD BoundingBox bounding_box(bool pos_side) const;

  double z0_;
};

//==============================================================================
//! A general plane.
//
//! The plane is described by the equation \f$A x + B y + C z - D = 0\f$
//==============================================================================

class SurfacePlane
{
public:
  explicit SurfacePlane(pugi::xml_node surf_node);

  SurfacePlane(SurfacePlane&& other) = default;

  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double A_, B_, C_, D_;
};

//==============================================================================
//! A cylinder aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceXCylinder
{
public:
  explicit SurfaceXCylinder(pugi::xml_node surf_node);
  SurfaceXCylinder(SurfaceXCylinder&& other) = default;
  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;
  HD BoundingBox bounding_box(bool pos_side) const;

  double y0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceYCylinder
{
public:
  explicit SurfaceYCylinder(pugi::xml_node surf_node);
  SurfaceYCylinder(SurfaceYCylinder&& other) = default;
  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;
  HD BoundingBox bounding_box(bool pos_side) const;

  double x0_, z0_, radius_;
};

//==============================================================================
//! A cylinder aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceZCylinder
{
public:
  explicit SurfaceZCylinder(pugi::xml_node surf_node);
  SurfaceZCylinder(SurfaceZCylinder&& other) = default;

  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;
  HD BoundingBox bounding_box(bool pos_side) const;

  double x0_, y0_, radius_;
};

//==============================================================================
//! A sphere.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 + (z - z_0)^2 - R^2 = 0\f$
//==============================================================================

class SurfaceSphere
{
public:
  explicit SurfaceSphere(pugi::xml_node surf_node);
  SurfaceSphere(SurfaceSphere&& other) = default;
  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;
  HD BoundingBox bounding_box(bool pos_side) const;

  double x0_, y0_, z0_, radius_;
};

//==============================================================================
//! A cone aligned along the x-axis.
//
//! The cylinder is described by the equation
//! \f$(y - y_0)^2 + (z - z_0)^2 - R^2 (x - x_0)^2 = 0\f$
//==============================================================================

class SurfaceXCone
{
public:
  explicit SurfaceXCone(pugi::xml_node surf_node);
  SurfaceXCone(SurfaceXCone&& other) = default;
  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the y-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (z - z_0)^2 - R^2 (y - y_0)^2 = 0\f$
//==============================================================================

class SurfaceYCone
{
public:
  explicit SurfaceYCone(pugi::xml_node surf_node);
  SurfaceYCone(SurfaceYCone&& other) = default;
  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A cone aligned along the z-axis.
//
//! The cylinder is described by the equation
//! \f$(x - x_0)^2 + (y - y_0)^2 - R^2 (z - z_0)^2 = 0\f$
//==============================================================================

class SurfaceZCone
{
public:
  explicit SurfaceZCone(pugi::xml_node surf_node);
  SurfaceZCone(SurfaceZCone&& other) = default;
  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;

  double x0_, y0_, z0_, radius_sq_;
};

//==============================================================================
//! A general surface described by a quadratic equation.
//
//! \f$A x^2 + B y^2 + C z^2 + D x y + E y z + F x z + G x + H y + J z + K = 0\f$
//==============================================================================

class SurfaceQuadric
{
public:
  explicit SurfaceQuadric(pugi::xml_node surf_node);
  SurfaceQuadric(SurfaceQuadric&& other) = default;

  HD double evaluate(Position const& r) const;
  HD double distance(Position const& r, Direction const& u, bool coincident) const;
  HD Direction normal(Position const& r) const;
  void to_hdf5_inner(hid_t group_id) const;

  // Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Jz + K = 0
  double A_, B_, C_, D_, E_, F_, G_, H_, J_, K_;
};

//==============================================================================
//! A geometry primitive used to define regions of 3D space.
//==============================================================================

/*
 * Tagged union defining any of the quadratic CSG surface types.
 */
class Surface
{
public:

  int id_; //!< Unique ID
  string name_; //!< User-defined name
  unique_ptr<BoundaryCondition> bc_ {nullptr}; //!< Boundary condition
  bool surf_source_ {false}; //!< Activate source banking for the surface?

  union SurfaceStorage {
    SurfaceXPlane xp;
    SurfaceYPlane yp;
    SurfaceZPlane zp;
    SurfacePlane   p;
    SurfaceXCylinder xc;
    SurfaceYCylinder yc;
    SurfaceZCylinder zc;
    SurfaceSphere sph;
    SurfaceXCone xco;
    SurfaceYCone yco;
    SurfaceZCone zco;
    SurfaceQuadric q;
    SurfaceStorage() {memset(this, 0, sizeof(SurfaceStorage));}
  } storage_;

  enum class SurfaceType {
    xplane,
    yplane,
    zplane,
    plane,
    xcylinder,
    ycylinder,
    zcylinder,
    sphere,
    xcone,
    ycone,
    zcone,
    quadric
  } type_;

  explicit Surface(pugi::xml_node surf_node, SurfaceType type);
  Surface() = default;

  Surface(Surface&& other) = default;

  //! Determine the direction of a ray reflected from the surface.
  //! \param[in] r The point at which the ray is incident.
  //! \param[in] u Incident direction of the ray
  //! \param[inout] p Pointer to the particle
  //! \return Outgoing direction of the ray
  HD Direction reflect(Position r, Direction u, Particle* p) const {
    Direction n;
    switch (type_) {
      case SurfaceType::xplane:
        n = storage_.xp.normal(r);
        break;
      case SurfaceType::yplane:
        n = storage_.yp.normal(r);
        break;
      case SurfaceType::zplane:
        n = storage_.zp.normal(r);
        break;
      case SurfaceType::plane:
        n = storage_.p.normal(r);
        break;
      case SurfaceType::xcylinder:
        n = storage_.xc.normal(r);
        break;
      case SurfaceType::ycylinder:
        n = storage_.yc.normal(r);
        break;
      case SurfaceType::zcylinder:
        n = storage_.zc.normal(r);
        break;
      case SurfaceType::sphere:
        n = storage_.sph.normal(r);
        break;
      case SurfaceType::xcone:
        n = storage_.xco.normal(r);
        break;
      case SurfaceType::ycone:
        n = storage_.yco.normal(r);
        break;
      case SurfaceType::zcone:
        n = storage_.zco.normal(r);
        break;
      case SurfaceType::quadric:
        n = storage_.q.normal(r);
    }
    return u.reflect(n);
  }

  HD Direction diffuse_reflect(
    Position r, Direction u, uint64_t* seed) const {
    Direction n;
    switch (type_) {
      case SurfaceType::xplane:
        n = storage_.xp.normal(r);
        break;
      case SurfaceType::yplane:
        n = storage_.yp.normal(r);
        break;
      case SurfaceType::zplane:
        n = storage_.zp.normal(r);
        break;
      case SurfaceType::plane:
        n = storage_.p.normal(r);
        break;
      case SurfaceType::xcylinder:
        n = storage_.xc.normal(r);
        break;
      case SurfaceType::ycylinder:
        n = storage_.yc.normal(r);
        break;
      case SurfaceType::zcylinder:
        n = storage_.zc.normal(r);
        break;
      case SurfaceType::sphere:
        n = storage_.sph.normal(r);
        break;
      case SurfaceType::xcone:
        n = storage_.xco.normal(r);
        break;
      case SurfaceType::ycone:
        n = storage_.yco.normal(r);
        break;
      case SurfaceType::zcone:
        n = storage_.zco.normal(r);
        break;
      case SurfaceType::quadric:
        n = storage_.q.normal(r);
    }

		n /= n.norm();
		const double projection = n.dot(u);

		// sample from inverse function, u=sqrt(rand) since p(u)=2u, so F(u)=u^2
		const double mu = (projection>=0.0) ?
										-std::sqrt(prn(seed)) : std::sqrt(prn(seed));

		// sample azimuthal distribution uniformly
		u = rotate_angle(n, mu, nullptr, seed);

		// normalize the direction
		return u/u.norm();
  }

  //! Evaluate the equation describing the surface.
  //!
  //! Surfaces can be described by some function f(x, y, z) = 0.  This member
  //! function evaluates that mathematical function.
  //! \param r A 3D Cartesian coordinate.
  HD double evaluate(Position const& r) const {
    switch (type_) {
      case SurfaceType::xplane:
        return storage_.xp.evaluate(r);
      case SurfaceType::yplane:
        return storage_.yp.evaluate(r);
      case SurfaceType::zplane:
        return storage_.zp.evaluate(r);
      case SurfaceType::plane:
        return storage_.p.evaluate(r);
      case SurfaceType::xcylinder:
        return storage_.xc.evaluate(r);
      case SurfaceType::ycylinder:
        return storage_.yc.evaluate(r);
      case SurfaceType::zcylinder:
        return storage_.zc.evaluate(r);
      case SurfaceType::sphere:
        return storage_.sph.evaluate(r);
      case SurfaceType::xcone:
        return storage_.xco.evaluate(r);
      case SurfaceType::ycone:
        return storage_.yco.evaluate(r);
      case SurfaceType::zcone:
        return storage_.zco.evaluate(r);
      case SurfaceType::quadric:
        return storage_.q.evaluate(r);
    }
    return 0.0;
  }

  //! Determine which side of a surface a point lies on.
  //! \param r The 3D Cartesian coordinate of a point.
  //! \param u A direction used to "break ties" and pick a sense when the
  //!   point is very close to the surface.
  //! \return true if the point is on the "positive" side of the surface and
  //!   false otherwise.
  HD bool sense(Position const& r, Direction const& u) const {
    double f;// = evaluate(r);
    double udotnormalr; //=u.dot(normal(r)) > 0.0
    // Note: avoids the use of two jump tables!
    switch (type_) {
      case SurfaceType::xplane:
        f = storage_.xp.evaluate(r);
        udotnormalr = u.dot(storage_.xp.normal(r));
        break;
      case SurfaceType::yplane:
        f = storage_.yp.evaluate(r);
        udotnormalr = u.dot(storage_.yp.normal(r));
        break;
      case SurfaceType::zplane:
        f = storage_.zp.evaluate(r);
        udotnormalr = u.dot(storage_.zp.normal(r));
        break;
      case SurfaceType::plane:
        f = storage_.p.evaluate(r);
        udotnormalr = u.dot(storage_.p.normal(r));
        break;
      case SurfaceType::xcylinder:
        f = storage_.xc.evaluate(r);
        udotnormalr = u.dot(storage_.xc.normal(r));
        break;
      case SurfaceType::ycylinder:
        f = storage_.yc.evaluate(r);
        udotnormalr = u.dot(storage_.yc.normal(r));
        break;
      case SurfaceType::zcylinder:
        f = storage_.zc.evaluate(r);
        udotnormalr = u.dot(storage_.zc.normal(r));
        break;
      case SurfaceType::sphere:
        f = storage_.sph.evaluate(r);
        udotnormalr = u.dot(storage_.sph.normal(r));
        break;
      case SurfaceType::xcone:
        f = storage_.xco.evaluate(r);
        udotnormalr = u.dot(storage_.xco.normal(r));
        break;
      case SurfaceType::ycone:
        f = storage_.yco.evaluate(r);
        udotnormalr = u.dot(storage_.yco.normal(r));
        break;
      case SurfaceType::zcone:
        f = storage_.zco.evaluate(r);
        udotnormalr = u.dot(storage_.zco.normal(r));
        break;
      case SurfaceType::quadric:
        f = storage_.q.evaluate(r);
        udotnormalr = u.dot(storage_.q.normal(r));
    }
    const bool c1 = std::abs(f) < FP_COINCIDENT;
    const bool c2 = udotnormalr > 0.0;
    const bool c3 = f > 0.0;
    return (c1 && c2) || (!c1 && c3); // le epic branch-free logic
  }

  //! Compute the distance between a point and the surface along a ray.
  //! \param r A 3D Cartesian coordinate.
  //! \param u The direction of the ray.
  //! \param coincident A hint to the code that the given point should lie
  //!   exactly on the surface.
  HD double distance(
    Position const& r, Direction const& u, bool coincident) {
    switch (type_) {
      case SurfaceType::xplane:
        return storage_.xp.distance(r, u, coincident);
      case SurfaceType::yplane:
        return storage_.yp.distance(r, u, coincident);
      case SurfaceType::zplane:
        return storage_.zp.distance(r, u, coincident);
      case SurfaceType::plane:
        return storage_.p.distance(r, u, coincident);
      case SurfaceType::xcylinder:
        return storage_.xc.distance(r, u, coincident);
      case SurfaceType::ycylinder:
        return storage_.yc.distance(r, u, coincident);
      case SurfaceType::zcylinder:
        return storage_.zc.distance(r, u, coincident);
      case SurfaceType::sphere:
        return storage_.sph.distance(r, u, coincident);
      case SurfaceType::xcone:
        return storage_.xco.distance(r, u, coincident);
      case SurfaceType::ycone:
        return storage_.yco.distance(r, u, coincident);
      case SurfaceType::zcone:
        return storage_.zco.distance(r, u, coincident);
      case SurfaceType::quadric:
        return storage_.q.distance(r, u, coincident);
    }
    return 0.0;
  }

  //! Compute the local outward normal direction of the surface.
  //! \param r A 3D Cartesian coordinate.
  //! \return Normal direction
  HD Direction normal(Position const& r) const {
    switch (type_) {
      case SurfaceType::xplane:
        return storage_.xp.normal(r);
      case SurfaceType::yplane:
        return storage_.yp.normal(r);
      case SurfaceType::zplane:
        return storage_.zp.normal(r);
      case SurfaceType::plane:
        return storage_.p.normal(r);
      case SurfaceType::xcylinder:
        return storage_.xc.normal(r);
      case SurfaceType::ycylinder:
        return storage_.yc.normal(r);
      case SurfaceType::zcylinder:
        return storage_.zc.normal(r);
      case SurfaceType::sphere:
        return storage_.sph.normal(r);
      case SurfaceType::xcone:
        return storage_.xco.normal(r);
      case SurfaceType::ycone:
        return storage_.yco.normal(r);
      case SurfaceType::zcone:
        return storage_.zco.normal(r);
      case SurfaceType::quadric:
        return storage_.q.normal(r);
    }
    return {};
  }

  //! Write all information needed to reconstruct the surface to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

  void to_hdf5_inner(hid_t group_id) const {
    switch (type_) {
      case SurfaceType::xplane:
        storage_.xp.to_hdf5_inner(group_id);
        break;
      case SurfaceType::yplane:
        storage_.yp.to_hdf5_inner(group_id);
        break;
      case SurfaceType::zplane:
        storage_.zp.to_hdf5_inner(group_id);
        break;
      case SurfaceType::plane:
        storage_.p.to_hdf5_inner(group_id);
        break;
      case SurfaceType::xcylinder:
        storage_.xc.to_hdf5_inner(group_id);
        break;
      case SurfaceType::ycylinder:
        storage_.yc.to_hdf5_inner(group_id);
        break;
      case SurfaceType::zcylinder:
        storage_.zc.to_hdf5_inner(group_id);
        break;
      case SurfaceType::sphere:
        storage_.sph.to_hdf5_inner(group_id);
        break;
      case SurfaceType::xcone:
        storage_.xco.to_hdf5_inner(group_id);
        break;
      case SurfaceType::ycone:
        storage_.yco.to_hdf5_inner(group_id);
        break;
      case SurfaceType::zcone:
        storage_.zco.to_hdf5_inner(group_id);
        break;
      case SurfaceType::quadric:
        storage_.q.to_hdf5_inner(group_id);
    }
  }

  //! Get the BoundingBox for this surface.
  HD BoundingBox bounding_box(bool pos_side) const {
    switch (type_) {
      case SurfaceType::xplane:
        return storage_.xp.bounding_box(pos_side);
      case SurfaceType::yplane:
        return storage_.yp.bounding_box(pos_side);
      case SurfaceType::zplane:
        return storage_.zp.bounding_box(pos_side);
      case SurfaceType::plane:
        return {};
      case SurfaceType::xcylinder:
        return storage_.xc.bounding_box(pos_side);
      case SurfaceType::ycylinder:
        return storage_.yc.bounding_box(pos_side);
      case SurfaceType::zcylinder:
        return storage_.zc.bounding_box(pos_side);
      case SurfaceType::sphere:
        return storage_.sph.bounding_box(pos_side);
      case SurfaceType::xcone:
        return {};
      case SurfaceType::ycone:
        return {};
      case SurfaceType::zcone:
        return {};
      case SurfaceType::quadric:
        return {};
    }
    return {};
  }
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_surfaces(pugi::xml_node node);

void free_memory_surfaces();

} // namespace openmc
#endif // OPENMC_SURFACE_H
