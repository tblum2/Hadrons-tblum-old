#include <Hadrons/Modules/MAction/NaiveStaggered.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MAction;

template class Grid::Hadrons::MAction::TNaiveStaggered<STAGIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MAction::TNaiveStaggered<STAGIMPLF>;
#endif
