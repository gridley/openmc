#ifndef OPENMC_CELL_H
#define OPENMC_CELL_H

#include <cstdint>
#include <functional> // for hash
#include <limits>
#include <unordered_map>

#include <gsl/gsl>
#include "hdf5.h"
#include "pugixml.hpp"
#include "dagmc.h"

#include "openmc/constants.h"
#include "openmc/memory.h" // for unique_ptr
#include "openmc/string.h"
#include "openmc/vector.h"
#include "openmc/neighbor_list.h"
#include "openmc/position.h"
#include "openmc/surface.h"
#include "openmc/vector.h"

namespace openmc {

//==============================================================================
// Constants
//==============================================================================

enum class Fill {
  MATERIAL,
  UNIVERSE,
  LATTICE
};

// TODO: Convert to enum
constexpr int32_t OP_LEFT_PAREN   {std::numeric_limits<int32_t>::max()};
constexpr int32_t OP_RIGHT_PAREN  {std::numeric_limits<int32_t>::max() - 1};
constexpr int32_t OP_COMPLEMENT   {std::numeric_limits<int32_t>::max() - 2};
constexpr int32_t OP_INTERSECTION {std::numeric_limits<int32_t>::max() - 3};
constexpr int32_t OP_UNION        {std::numeric_limits<int32_t>::max() - 4};

//==============================================================================
// Global variables
//==============================================================================

class CSGCell;
class ParentCell;
class CellInstance;
class Universe;
class UniversePartitioner;

namespace model {
  extern std::unordered_map<int32_t, int32_t> cell_map;
  extern vector<unique_ptr<CSGCell>> cells;

  extern std::unordered_map<int32_t, int32_t> universe_map;
  extern vector<unique_ptr<Universe>> universes;
} // namespace model

#ifdef __CUDACC__
namespace gpu {
extern __constant__ unique_ptr<CSGCell>* cells;
extern __constant__ unique_ptr<Universe>* universes;
} // namespace gpu
#endif

//==============================================================================
//! A geometry primitive that fills all space and contains cells.
//==============================================================================

class Universe
{
public:
  Universe() = default;

  int32_t id_;                  //!< Unique ID
  vector<int32_t> cells_;       //!< Cells within this universe

  //! \brief Write universe information to an HDF5 group.
  //! \param group_id An HDF5 group id.
  void to_hdf5(hid_t group_id) const;

  BoundingBox bounding_box() const;

  unique_ptr<UniversePartitioner> partitioner_;
};

//==============================================================================
//==============================================================================

class Cell {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors, factory functions

  explicit Cell(pugi::xml_node cell_node);
  Cell() = default;

  //----------------------------------------------------------------------------
  // Accessors

  //! Get the temperature of a cell instance
  //! \param[in] instance Instance index. If -1 is given, the temperature for
  //!   the first instance is returned.
  //! \return Temperature in [K]
  HD double temperature(int32_t instance = -1) const;

  //! Set the temperature of a cell instance
  //! \param[in] T Temperature in [K]
  //! \param[in] instance Instance index. If -1 is given, the temperature for
  //!   all instances is set.
  //! \param[in] set_contained If this cell is not filled with a material,
  //!   collect all contained cells with material fills and set their
  //!   temperatures.
  void set_temperature(double T, int32_t instance = -1, bool set_contained = false);

  //! Get the name of a cell
  //! \return Cell name
  const string& name() const { return name_; };

  //! Set the temperature of a cell instance
  //! \param[in] name Cell name
  void set_name(const std::string& name) { name_ = name; };

  //! Get all cell instances contained by this cell
  //! \return Map with cell indexes as keys and instances as values
  std::unordered_map<int32_t, vector<int32_t>> get_contained_cells() const;

protected:
  void get_contained_cells_inner(
    std::unordered_map<int32_t, vector<int32_t>>& contained_cells,
    vector<ParentCell>& parent_cells) const;

public:
  //----------------------------------------------------------------------------
  // Data members

  int32_t id_;                //!< Unique ID
  string name_;               //!< User-defined name
  Fill type_;                 //!< Material, universe, or lattice
  int32_t universe_;          //!< Universe # this cell is in
  int32_t fill_;              //!< Universe # filling this cell
  int32_t n_instances_{0};    //!< Number of instances of this cell

  //! \brief Index corresponding to this cell in distribcell arrays
  int distribcell_index_{C_NONE};

  //! \brief Material(s) within this cell.
  //!
  //! May be multiple materials for distribcell.
  vector<int32_t> material_;

  //! \brief Temperature(s) within this cell.
  //!
  //! The stored values are actually sqrt(k_Boltzmann * T) for each temperature
  //! T. The units are sqrt(eV).
  vector<double> sqrtkT_;

  //! Definition of spatial region as Boolean expression of half-spaces
  vector<int32_t> region_;

  //! \brief Neighboring cells in the same universe.
  NeighborList neighbors_;

  Position translation_ {0, 0, 0}; //!< Translation vector for filled universe

  //! \brief Rotational tranfsormation of the filled universe.
  //
  //! The vector is empty if there is no rotation. Otherwise, the first 9 values
  //! give the rotation matrix in row-major order. When the user specifies
  //! rotation angles about the x-, y- and z- axes in degrees, these values are
  //! also present at the end of the vector, making it of length 12.
  vector<double> rotation_;

  vector<int32_t> offset_; //!< Distribcell offset table
};

struct CellInstanceItem {
  int32_t index {-1};        //! Index into global cells array
  int     lattice_indx{-1};  //! Flat index value of the lattice cell
};

//==============================================================================

class CSGCell : public Cell
{
public:
  explicit CSGCell(pugi::xml_node cell_node);
  CSGCell() = default;

  // This is pretty performance critical so it's inlined
  bool HD contains(Position const& r, Direction const& u, int32_t const& on_surface) const
  {
    for (int32_t token : region_) {
      // Assume that no tokens are operators. Evaluate the sense of particle with
      // respect to the surface and see if the token matches the sense. If the
      // particle's surface attribute is set and matches the token, that
      // overrides the determination based on sense().
      if (token == on_surface) {
      } else if (-token == on_surface) {
        return false;
      } else {
        // Note the off-by-one indexing
#ifdef __CUDA_ARCH__
        bool sense = gpu::surfaces[abs(token) - 1].sense(r, u);
#else
        bool sense = model::surfaces[abs(token)-1].sense(r, u);
#endif
        if (sense != (token > 0)) {return false;}
      }
    }
    return true;
  }

  std::pair<double, int32_t> HD distance(
    Position r, Direction u, int32_t on_surface, Particle* p) const
  {
      #ifdef __CUDA_ARCH__
      using gpu::surfaces;
      #else
      using model::surfaces;
      #endif

      double min_dist {INFTY};
      constexpr int32_t default_cell_index = std::numeric_limits<int32_t>::max();
      int32_t i_surf {default_cell_index};

      for (int32_t token : region_) {

        // Calculate the distance to this surface.
        // Note the off-by-one indexing
        bool coincident {std::abs(token) == std::abs(on_surface)};
        double d {surfaces[abs(token) - 1].distance(r, u, coincident)};

        // Check if this distance is the new minimum.
        if (d < min_dist) {
          if (min_dist - d >= FP_PRECISION*min_dist) {
            min_dist = d;
            i_surf = -token;
          }
        }
      }

    return {min_dist, i_surf};
  }

  void to_hdf5(hid_t group_id) const;

  BoundingBox bounding_box() const;

protected:
  BoundingBox bounding_box_simple() const;
};

//==============================================================================
//! Speeds up geometry searches by grouping cells in a search tree.
//
//! Currently this object only works with universes that are divided up by a
//! bunch of z-planes.  It could be generalized to other planes, cylinders,
//! and spheres.
//==============================================================================

class UniversePartitioner
{
public:
  explicit UniversePartitioner(const Universe& univ);

  //! Return the list of cells that could contain the given coordinates.
  HD const vector<int32_t>& get_cells(Position r, Direction u) const;

private:
  //! A sorted vector of indices to surfaces that partition the universe
  vector<int32_t> surfs_;

  //! Vectors listing the indices of the cells that lie within each partition
  //
  //! There are n+1 partitions with n surfaces.  `partitions_.front()` gives the
  //! cells that lie on the negative side of `surfs_.front()`.
  //! `partitions_.back()` gives the cells that lie on the positive side of
  //! `surfs_.back()`.  Otherwise, `partitions_[i]` gives cells sandwiched
  //! between `surfs_[i-1]` and `surfs_[i]`.
  vector<vector<int32_t>> partitions_;
};


//==============================================================================
//! Define a containing (parent) cell
//==============================================================================

struct ParentCell {
  gsl::index cell_index;
  gsl::index lattice_index;
};

//==============================================================================
//! Define an instance of a particular cell
//==============================================================================

struct CellInstance {
  //! Check for equality
  bool operator==(const CellInstance& other) const
  { return index_cell == other.index_cell && instance == other.instance; }

  gsl::index index_cell;
  gsl::index instance;
};

struct CellInstanceHash {
  std::size_t operator()(const CellInstance& k) const
  {
    return 4096*k.index_cell + k.instance;
  }
};

//==============================================================================
// Non-member functions
//==============================================================================

void read_cells(pugi::xml_node node);

#ifdef DAGMC
int32_t next_cell(DAGCell* cur_cell, DAGSurface* surf_xed);
#endif

} // namespace openmc
#endif // OPENMC_CELL_H
