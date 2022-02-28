#ifndef OPENMC_SURFACE_H
#define OPENMC_SURFACE_H

#include "openmc/string.h"
#include <limits> // For numeric_limits
#include <unordered_map>
#include <cstring> // memset

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

  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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

  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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
  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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

  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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
  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;

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
  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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

  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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
  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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
  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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
  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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
  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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

  inline HD double evaluate(Position const& r) const;
  inline HD double distance(Position const& r, Direction const& u, bool coincident) const;
  inline HD Direction normal(Position const& r) const;
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

// The template parameter indicates the axis normal to the plane.
template<int i>
double HD axis_aligned_plane_distance(
  Position const& r, Direction const& u, bool coincident, double const& offset)
{
  const double f = offset - r[i];
  if (coincident || std::abs(f) < FP_COINCIDENT || u[i] == 0.0) return INFTY;
  const double d = f / u[i];
  if (d < 0.0) return INFTY;
  return d;
}

HD double SurfaceXPlane::evaluate(Position const& r) const
{
  return r.x - x0_;
}

HD double SurfaceXPlane::distance(
  Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_plane_distance<0>(r, u, coincident, x0_);
}

HD Direction SurfaceXPlane::normal(Position const& r) const
{
  return {1., 0., 0.};
}

HD double SurfaceYPlane::evaluate(Position const& r) const
{
  return r.y - y0_;
}

HD double SurfaceYPlane::distance(
  Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_plane_distance<1>(r, u, coincident, y0_);
}

HD Direction SurfaceYPlane::normal(Position const& r) const
{
  return {0., 1., 0.};
}


HD double SurfaceZPlane::evaluate(Position const& r) const
{
  return r.z - z0_;
}

HD double SurfaceZPlane::distance(
  Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_plane_distance<2>(r, u, coincident, z0_);
}

HD Direction SurfaceZPlane::normal(Position const& r) const
{
  return {0., 0., 1.};
}

double HD SurfacePlane::evaluate(Position const& r) const
{
  return A_*r.x + B_*r.y + C_*r.z - D_;
}

double HD SurfacePlane::distance(Position const& r, Direction const& u, bool coincident) const
{
  const double f = A_*r.x + B_*r.y + C_*r.z - D_;
  const double projection = A_*u.x + B_*u.y + C_*u.z;
  if (coincident || std::abs(f) < FP_COINCIDENT || projection == 0.0) {
    return INFTY;
  } else {
    const double d = -f / projection;
    if (d < 0.0) return INFTY;
    return d;
  }
}

Direction HD SurfacePlane::normal(Position const& r) const
{
  return {A_, B_, C_};
}

// The template parameters indicate the axes perpendicular to the axis of the
// cylinder.  offset1 and offset2 should correspond with i1 and i2,
// respectively.
template<int i1, int i2>
double HD axis_aligned_cylinder_evaluate(
  Position const& r, double const& offset1, double const& offset2, double const& radius)
{
  const double r1 = r.get<i1>() - offset1;
  const double r2 = r.get<i2>() - offset2;
  return r1*r1 + r2*r2 - radius*radius;
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3>
double HD axis_aligned_cylinder_distance(Position const& r, Direction const& u,
  bool coincident, double const& offset1, double const& offset2, double const& radius)
{
  const double a = 1.0 - u.get<i1>() * u.get<i1>(); // u^2 + v^2
  if (a == 0.0) return INFTY;

  const double r2 = r.get<i2>() - offset1;
  const double r3 = r.get<i3>() - offset2;
  const double k = r2 * u.get<i2>() + r3 * u.get<i3>();
  const double c = r2*r2 + r3*r3 - radius*radius;
  const double quad = k*k - a*c;

  if (quad < 0.0) {
    // No intersection with cylinder.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the cylinder, thus one distance is positive/negative
    // and the other is zero. The sign of k determines if we are facing in or
    // out.
    if (k >= 0.0) {
      return INFTY;
    } else {
      return (-k + sqrt(quad)) / a;
    }

  } else if (c < 0.0) {
    // Particle is inside the cylinder, thus one distance must be negative
    // and one must be positive. The positive distance will be the one with
    // negative sign on sqrt(quad).
    return (-k + sqrt(quad)) / a;

  } else {
    // Particle is outside the cylinder, thus both distances are either
    // positive or negative. If positive, the smaller distance is the one
    // with positive sign on sqrt(quad).
    const double d = (-k - sqrt(quad)) / a;
    if (d < 0.0) return INFTY;
    return d;
  }
}

// The first template parameter indicates which axis the cylinder is aligned to.
// The other two parameters indicate the other two axes.  offset1 and offset2
// should correspond with i2 and i3, respectively.
template<int i1, int i2, int i3>
Direction HD axis_aligned_cylinder_normal(
  Position const& r, double const& offset1, double const& offset2)
{
  Direction u;
  u.get<i2>() = 2.0 * (r.get<i2>() - offset1);
  u.get<i3>() = 2.0 * (r.get<i3>() - offset2);
  u.get<i1>() = 0.0;
  return u;
}

HD double SurfaceYCylinder::evaluate(Position const& r) const
{
  return axis_aligned_cylinder_evaluate<0, 2>(r, x0_, z0_, radius_);
}

HD double SurfaceYCylinder::distance(
  Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_cylinder_distance<1, 0, 2>(r, u, coincident, x0_, z0_,
                                                 radius_);
}

HD Direction SurfaceYCylinder::normal(Position const& r) const
{
  return axis_aligned_cylinder_normal<1, 0, 2>(r, x0_, z0_);
}

HD double SurfaceXCylinder::evaluate(Position const& r) const
{
  return axis_aligned_cylinder_evaluate<1, 2>(r, y0_, z0_, radius_);
}

HD double SurfaceXCylinder::distance(
  Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_cylinder_distance<0, 1, 2>(r, u, coincident, y0_, z0_,
                                                 radius_);
}

HD Direction SurfaceXCylinder::normal(Position const& r) const
{
  return axis_aligned_cylinder_normal<0, 1, 2>(r, y0_, z0_);
}

HD double SurfaceZCylinder::evaluate(Position const& r) const
{
  return axis_aligned_cylinder_evaluate<0, 1>(r, x0_, y0_, radius_);
}

HD double SurfaceZCylinder::distance(
  Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_cylinder_distance<2, 0, 1>(r, u, coincident, x0_, y0_,
                                                 radius_);
}

HD Direction SurfaceZCylinder::normal(Position const& r) const
{
  return axis_aligned_cylinder_normal<2, 0, 1>(r, x0_, y0_);
}

HD double SurfaceSphere::evaluate(Position const& r) const
{
  const double x = r.x - x0_;
  const double y = r.y - y0_;
  const double z = r.z - z0_;
  return x*x + y*y + z*z - radius_*radius_;
}

HD double SurfaceSphere::distance(
  Position const& r, Direction const& u, bool coincident) const
{
  const double x = r.x - x0_;
  const double y = r.y - y0_;
  const double z = r.z - z0_;
  const double k = x*u.x + y*u.y + z*u.z;
  const double c = x*x + y*y + z*z - radius_*radius_;
  const double quad = k*k - c;

  if (quad < 0.0) {
    // No intersection with sphere.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the sphere, thus one distance is positive/negative and
    // the other is zero. The sign of k determines if we are facing in or out.
    if (k >= 0.0) {
      return INFTY;
    } else {
      return -k + sqrt(quad);
    }

  } else if (c < 0.0) {
    // Particle is inside the sphere, thus one distance must be negative and
    // one must be positive. The positive distance will be the one with
    // negative sign on sqrt(quad)
    return -k + sqrt(quad);

  } else {
    // Particle is outside the sphere, thus both distances are either positive
    // or negative. If positive, the smaller distance is the one with positive
    // sign on sqrt(quad).
    const double d = -k - sqrt(quad);
    if (d < 0.0) return INFTY;
    return d;
  }
}

HD Direction SurfaceSphere::normal(Position const& r) const
{
  return {2.0*(r.x - x0_), 2.0*(r.y - y0_), 2.0*(r.z - z0_)};
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3>
double HD axis_aligned_cone_evaluate(
  Position const& r, double const& offset1, double const& offset2, double const& offset3, double const& radius_sq)
{
  const double r1 = r.get<i1>() - offset1;
  const double r2 = r.get<i2>() - offset2;
  const double r3 = r.get<i3>() - offset3;
  return r2*r2 + r3*r3 - radius_sq*r1*r1;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3>
double HD axis_aligned_cone_distance(Position const& r, Direction const& u, bool coincident,
  double const& offset1, double const& offset2, double const& offset3, double const& radius_sq)
{
  const double r1 = r.get<i1>() - offset1;
  const double r2 = r.get<i2>() - offset2;
  const double r3 = r.get<i3>() - offset3;
  const double a = u.get<i2>() * u.get<i2>() + u.get<i3>() * u.get<i3>() -
                   radius_sq * u.get<i1>() * u.get<i1>();
  const double k =
    r2 * u.get<i2>() + r3 * u.get<i3>() - radius_sq * r1 * u.get<i1>();
  const double c = r2*r2 + r3*r3 - radius_sq*r1*r1;
  double quad = k*k - a*c;

  double d;

  if (quad < 0.0) {
    // No intersection with cone.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the cone, thus one distance is positive/negative
    // and the other is zero. The sign of k determines if we are facing in or
    // out.
    if (k >= 0.0) {
      d = (-k - sqrt(quad)) / a;
    } else {
      d = (-k + sqrt(quad)) / a;
    }

  } else {
    // Calculate both solutions to the quadratic.
    quad = sqrt(quad);
    d = (-k - quad) / a;
    const double b = (-k + quad) / a;

    // Determine the smallest positive solution.
    if (d < 0.0) {
      if (b > 0.0) d = b;
    } else {
      if (b > 0.0) {
        if (b < d) d = b;
      }
    }
  }

  // If the distance was negative, set boundary distance to infinity.
  if (d <= 0.0) return INFTY;
  return d;
}

// The first template parameter indicates which axis the cone is aligned to.
// The other two parameters indicate the other two axes.  offset1, offset2,
// and offset3 should correspond with i1, i2, and i3, respectively.
template<int i1, int i2, int i3>
Direction HD axis_aligned_cone_normal(
  Position const& r, double const& offset1, double const& offset2, double const& offset3, double const& radius_sq)
{
  Direction u;
  u.get<i1>() = -2.0 * radius_sq * (r.get<i1>() - offset1);
  u.get<i2>() = 2.0 * (r.get<i2>() - offset2);
  u.get<i3>() = 2.0 * (r.get<i3>() - offset3);
  return u;
}

HD double SurfaceXCone::evaluate(Position const& r) const
{
  return axis_aligned_cone_evaluate<0, 1, 2>(r, x0_, y0_, z0_, radius_sq_);
}

HD double SurfaceXCone::distance(Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_cone_distance<0, 1, 2>(r, u, coincident, x0_, y0_, z0_,
                                             radius_sq_);
}

HD Direction SurfaceXCone::normal(Position const& r) const
{
  return axis_aligned_cone_normal<0, 1, 2>(r, x0_, y0_, z0_, radius_sq_);
}

HD double SurfaceYCone::evaluate(Position const& r) const
{
  return axis_aligned_cone_evaluate<1, 0, 2>(r, y0_, x0_, z0_, radius_sq_);
}

HD double SurfaceYCone::distance(Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_cone_distance<1, 0, 2>(r, u, coincident, y0_, x0_, z0_,
                                             radius_sq_);
}

HD Direction SurfaceYCone::normal(Position const& r) const
{
  return axis_aligned_cone_normal<1, 0, 2>(r, y0_, x0_, z0_, radius_sq_);
}

HD double SurfaceZCone::evaluate(Position const& r) const
{
  return axis_aligned_cone_evaluate<2, 0, 1>(r, z0_, x0_, y0_, radius_sq_);
}

HD double SurfaceZCone::distance(Position const& r, Direction const& u, bool coincident) const
{
  return axis_aligned_cone_distance<2, 0, 1>(r, u, coincident, z0_, x0_, y0_,
                                             radius_sq_);
}

HD Direction SurfaceZCone::normal(Position const& r) const
{
  return axis_aligned_cone_normal<2, 0, 1>(r, z0_, x0_, y0_, radius_sq_);
}

double HD SurfaceQuadric::evaluate(Position const& r) const
{
  const double x = r.x;
  const double y = r.y;
  const double z = r.z;
  return x*(A_*x + D_*y + G_) +
         y*(B_*y + E_*z + H_) +
         z*(C_*z + F_*x + J_) + K_;
}

double HD SurfaceQuadric::distance(
  Position const& r, Direction const& ang, bool coincident) const
{
  const double &x = r.x;
  const double &y = r.y;
  const double &z = r.z;
  const double &u = ang.x;
  const double &v = ang.y;
  const double &w = ang.z;

  const double a = A_*u*u + B_*v*v + C_*w*w + D_*u*v + E_*v*w + F_*u*w;
  const double k = A_*u*x + B_*v*y + C_*w*z + 0.5*(D_*(u*y + v*x)
                   + E_*(v*z + w*y) + F_*(w*x + u*z) + G_*u + H_*v + J_*w);
  const double c = A_*x*x + B_*y*y + C_*z*z + D_*x*y + E_*y*z +  F_*x*z + G_*x
                   + H_*y + J_*z + K_;
  double quad = k*k - a*c;

  double d;

  if (quad < 0.0) {
    // No intersection with surface.
    return INFTY;

  } else if (coincident || std::abs(c) < FP_COINCIDENT) {
    // Particle is on the surface, thus one distance is positive/negative and
    // the other is zero. The sign of k determines which distance is zero and
    // which is not. Additionally, if a is zero, it means the particle is on
    // a plane-like surface.
    if (a == 0.0) {
      d = INFTY; // see the below explanation
    } else if (k >= 0.0) {
      d = (-k - sqrt(quad)) / a;
    } else {
      d = (-k + sqrt(quad)) / a;
    }

  } else if (a == 0.0) {
    // Given the orientation of the particle, the quadric looks like a plane in
    // this case, and thus we have only one solution despite potentially having
    // quad > 0.0. While the term under the square root may be real, in one
    // case of the +/- of the quadratic formula, 0/0 results, and in another, a
    // finite value over 0 results. Applying L'Hopital's to the 0/0 case gives
    // the below. Alternatively this can be found by simply putting a=0 in the
    // equation ax^2 + bx + c = 0.
    d = -0.5 * c / k;
  } else {
    // Calculate both solutions to the quadratic.
    quad = sqrt(quad);
    d = (-k - quad) / a;
    double b = (-k + quad) / a;

    // Determine the smallest positive solution.
    if (d < 0.0) {
      if (b > 0.0) d = b;
    } else {
      if (b > 0.0) {
        if (b < d) d = b;
      }
    }
  }

  // If the distance was negative, set boundary distance to infinity.
  if (d <= 0.0) return INFTY;
  return d;
}

Direction HD SurfaceQuadric::normal(Position const& r) const
{
  const double &x = r.x;
  const double &y = r.y;
  const double &z = r.z;
  return {2.0*A_*x + D_*y + F_*z + G_,
          2.0*B_*y + D_*x + E_*z + H_,
          2.0*C_*z + E_*y + F_*x + J_};
}

//==============================================================================
// Non-member functions
//==============================================================================

void read_surfaces(pugi::xml_node node);

void free_memory_surfaces();

} // namespace openmc
#endif // OPENMC_SURFACE_H
