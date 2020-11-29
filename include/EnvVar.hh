// Simple object for querying (at construction time) some environment
// variable (substituting a default fallback value if not found) and
// storing it for later access.

#ifndef MYRICUBE_ENVVAR_HH_
#define MYRICUBE_ENVVAR_HH_

#include "myricube.hh"

#include <stdint.h>
#include <stdlib.h>
#include <string>

namespace myricube {

// Environment variable parsed as 64 bit integer.
class EnvVar64
{
    int64_t var;

  public:
    EnvVar64(const char* name, int64_t fallback)
    {
        const char* env = getenv(name);
        if (env == nullptr) {
            var = fallback;
        }
        else {
            // Want to allow hexadecimal but I hate leading 0 octal.
            bool leading_zeros = false;
            while (*env == '0') {
                env++;
                leading_zeros = true;
            }
            if (*env == '\0' and leading_zeros) {
                var = 0;
            }
            else {
                char* endptr;
                #ifdef MYRICUBE_WINDOWS
                    var = strtoll(env, &endptr, 0);
                #else
                    var = strtol(env, &endptr, 0);
                #endif
                if (*endptr != '\0') {
                    throw std::runtime_error(
                        "Could not parse environment variable "
                        + std::string(env) + " as integer: '" + env + "'");
                }
            }
        }
    }

    operator int64_t() const {
        return var;
    }
};

// String environment variable.
class EnvVarString
{
    std::string var;

  public:
    EnvVarString(const char* name, const char* fallback)
    {
        const char* env = getenv(name);
        var = env ? env : fallback;
    }

    operator const std::string& () const {
        return var;
    }
};

// File name environment variable (only different on Windows).
#if !defined(MYRICUBE_WINDOWS)

using EnvVarFilename = EnvVarString;

#else
class EnvVarFilename
{
    filename_string var;

  public:
    EnvVarFilename(const char* name, const char* fallback)
    {
        filename_string empty;
        const wchar_t* env =
            _wgetenv(filename_concat_c_str(empty, name).c_str());
        if (env == nullptr) {
            var = filename_concat_c_str(empty, fallback);
        }
        else {
            var = env;
        }
    }

    operator const filename_string& () const {
        return var;
    }
};
#endif

} // end namespace myricube

#endif /* !MYRICUBE_ENVVAR_HH_ */
