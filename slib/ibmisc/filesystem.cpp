#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <ibmisc/ibmisc.hpp>
#include <ibmisc/filesystem.hpp>

namespace ibmisc {

EnvSearchPath::EnvSearchPath(std::string const &_env_var)
    : env_var(_env_var)
{
    std::string PATH(getenv(env_var.c_str()));
    boost::algorithm::split(path, PATH,
        boost::is_any_of(":;"));
}

std::string EnvSearchPath::locate(std::string const &file_name) const
{
    for (auto &dir : path) {
        auto fname(boost::filesystem::path(dir) / file_name);
        if (boost::filesystem::exists(fname))
            return fname.string();
    }
    (*ibmisc_error)(-1,
        "Cannot locate file %s in path $%s", file_name.c_str(), env_var.c_str());
}

}
