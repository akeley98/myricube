// Implement the global name-to-app-factory map.

#include "app.hh"

#include <unordered_map>
#include <stdio.h>

namespace myricube {

static std::unordered_map<std::string, App* (*) ()>& get_map()
{
    static std::unordered_map<std::string, App* (*) ()> map;
    return map;
}

AppNamer::AppNamer(std::string name, const char* filename, App* (*factory)())
{
    auto& map = get_map();

    auto it = map.find(name);
    if (it == map.end()) {
        map.emplace(std::move(name), factory);
    }
    else {
        fprintf(stderr, "%s:\nSkipped %s due to name collision.\n",
            filename, name.c_str());
    }
}

App* new_named_app(const std::string& name)
{
    auto& map = get_map();
    auto it = map.find(name);
    if (it != map.end()) {
        return it->second();
    }
    return nullptr;
}

void stderr_dump_app_names()
{
    for (auto& pair : get_map()) {
        fprintf(stderr, "%s\n", pair.first.c_str());
    }
}

} // end namespace
