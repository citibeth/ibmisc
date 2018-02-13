#ifndef IBMISC_FILESYSTEM_HPP
#define IBMISC_FILESYSTEM_HPP

#include <vector>
#include <string>

namespace ibmisc {


struct FileLocator {
    virtual std::string locate(std::string const &file_name) const = 0;
    ~FileLocator() {}
};

/** Locate existing files in a path given by an environnment variable. */
class EnvSearchPath : public FileLocator {

    std::string env_var;
    std::vector<std::string> path;

public:
    EnvSearchPath(std::string const &_env_var);

    std::string locate(std::string const &file_name) const;
};

}    // namespace ibmisc
#endif    // IBMISC_FILESYSTEM_HPP
