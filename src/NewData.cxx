#include <vector>

#include <NewData.hpp>
#include <util/Log.hpp>

namespace sisi4s {
namespace data {

static size_t SYMBOL_COUNTER = 0;

std::string gensym(std::string const &prefix) {
  std::stringstream s;
  s << prefix << "data-field-" << SYMBOL_COUNTER++;
  return s.str();
}

bool null(StorePair const &p) { return p.first == "" && p.second == nullptr; }

StorePair find_by_name(std::string const &name) {
  auto it = std::find_if(
      DATA_STORE.begin(),
      DATA_STORE.end(),
      [&name](StorePair &p) -> bool { return p.second->name == name; });
  return (it == DATA_STORE.end()) ? StorePair() : *it;
}

StorePair make_data() {
  std::string db_index = gensym("");
  DATA_STORE[db_index] = new Data();
  return {db_index, DATA_STORE[db_index]};
}

Data *getraw(std::string const &db_index) {
  auto ptr = DATA_STORE[db_index];
  if (ptr) return ptr;
  else {
    std::stringstream s;
    s << "Data with index '" << db_index << "' not found";
    throw std::domain_error(s.str());
  }
}

void setraw(std::string const &db_index, Data *ptr) {
  // auto _ptr = DATA_STORE[db_index];
  // if (_ptr)
  //   throw std::domain_error(
  //       _FORMAT("db: Data with index '%s' already exists",
  //       db_index.c_str()));
  DATA_STORE[db_index] = ptr;
}

#define INSTANTIATE(type, user_name)                                           \
  std::string Namer<type>::name() { return user_name; }                        \
  std::string Namer<std::vector<type>>::name() {                               \
    return "Vector of " + std::string(user_name);                              \
  }
INSTANTIATE(int, "Integer")
INSTANTIATE(int64_t, "64 bit integer")
INSTANTIATE(std::string, "String")
INSTANTIATE(bool, "Boolean value")
INSTANTIATE(double, "Double precission number")
INSTANTIATE(float, "Single precission number")
#undef INSTANTIATE

bool exists(std::string const &db_index) { return DATA_STORE[db_index]; }

} // namespace data
} // namespace sisi4s
