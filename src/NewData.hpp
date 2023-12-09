#ifndef ___NEWDATA_HPP_
#define ___NEWDATA_HPP_

#include <typeinfo>
#include <typeindex>
#include <sstream>
#include <map>
#include <algorithm>
#include <iostream>
#include <vector>

namespace sisi4s {
namespace data {

class Data {
public:
  std::type_index idx;
  void *data;
  size_t size;
  size_t state;
  std::string name;

  enum State {
    ALLOCATED = 1 << 1,
    TYPED = 1 << 2,
    FREE = 1 << 3,
  };

  bool is_state(size_t query) { return query & state; }

  Data()
      : idx(std::type_index(typeid(void)))
      , data(nullptr)
      , size(0)
      , state(FREE) {}
};

using DataStore = std::map<std::string, Data *>;
using StorePair = DataStore::value_type;
static DataStore DATA_STORE;

std::string gensym(std::string const &prefix);

StorePair make_data();

StorePair find_by_name(std::string const &name);

bool null(StorePair const &p);

Data *getraw(std::string const &db_index);

void setraw(std::string const &db_index, Data *ptr);

bool exists(std::string const &db_index);

template <typename F>
void declare_type(std::string const &db_index) {
  auto ptr = getraw(db_index);
  ptr->state |= Data::State::TYPED;
  if (!ptr) return;
  ptr->idx = std::type_index(typeid(F));
}

template <typename F>
class Namer {
public:
  static std::string name() { return std::type_index(typeid(F)).name(); }
};

template <typename F>
class Namer<std::vector<F>> {
public:
  static std::string name() { return "Vector of " + Namer<F>::name(); }
};

template <typename F>
std::string type_name() {
  // return std::type_index(typeid(F)).name();
  return Namer<F>::name();
}

#define INSTANTIATE(type)                                                      \
  template <>                                                                  \
  class Namer<type> {                                                          \
  public:                                                                      \
    static std::string name();                                                 \
  };                                                                           \
  template <>                                                                  \
  class Namer<std::vector<type>> {                                             \
  public:                                                                      \
    static std::string name();                                                 \
  };
INSTANTIATE(int)
INSTANTIATE(int64_t)
INSTANTIATE(std::string)
INSTANTIATE(bool)
INSTANTIATE(double)
INSTANTIATE(float)
#undef INSTANTIATE

template <typename F>
bool istype(std::string const &db_index) {
  auto ptr = getraw(db_index);
  const std::string name1 = ptr->idx.name(),
                    name2 = std::type_index(typeid(F)).name();
  return name1 == name2;
}

template <typename F>
F *get(std::string const &db_index) {
  auto ptr = getraw(db_index);
  const std::string name1 = ptr->idx.name(),
                    name2 = std::type_index(typeid(F)).name();
  if (!istype<F>(db_index)) {
    std::cout << "db: WARNING: "
              << "getting variable \"" << db_index << "\" with type \n\t"
              << name2 << " \nalthough in the database it has the type \n\t"
              << name1 << std::endl;
  }
  return static_cast<F *>(ptr->data);
}

template <typename F>
StorePair put(std::string const &name, F *const &value) {
  auto pair = make_data();
  const std::string db_index = pair.first;
  auto ptr = pair.second;
  ptr->state |= Data::State::ALLOCATED;
  ptr->state ^= Data::State::FREE;
  ptr->name = name;
  ptr->data = (void *)value;
  ptr->size = sizeof(F);
  setraw(db_index, ptr);
  declare_type<F>(db_index);
  return pair;
}

template <typename F>
StorePair put(std::string const &name, F const &value) {
  F *copy = new F(value);
  return put<F>(name, copy);
}

template <typename F>
StorePair put(F const &value) {
  const std::string name = gensym("var-");
  return put<F>(name, value);
}

template <typename F>
void free(std::string name) {
  auto ptr = get<F>(name);
  delete ptr;
}

} // namespace data

} // namespace sisi4s

#endif
