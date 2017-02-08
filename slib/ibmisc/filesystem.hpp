#include <vector>

namespace ibmisc {

/** Locate existing files in a path given by an environnment variable. */
class EnvSearchPath {

    std::string env_var;
    std::vector<std::string> path;

public:
    EnvSearchPath(std::string const &_env_var);

    std::string locate(std::string const &file_name);
};

}
