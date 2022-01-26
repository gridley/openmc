#include "openmc/array.h"

namespace openmc {

// OpenMC major, minor, and release numbers
constexpr int VERSION_MAJOR {1000};
constexpr int VERSION_MINOR {123};
constexpr int VERSION_RELEASE {1234};
constexpr bool VERSION_DEV {true};
constexpr std::array<int, 3> VERSION {VERSION_MAJOR, VERSION_MINOR, VERSION_RELEASE};

}
