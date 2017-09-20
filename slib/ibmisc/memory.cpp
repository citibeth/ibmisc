#include <ibmisc/memory.hpp>

namespace ibmisc {

void TmpAlloc::free() {
    for (auto ii(deleters.begin()); ii != deleters.end(); ++ii) {
        (*ii)();
    }
    deleters.clear();
}

}   // namespace ibmisc
