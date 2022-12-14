
#include <Data.hpp>
#include <sstream>

using namespace sisi4s;

Data::Data(std::string const &name_)
    : name(name_)
    , stage(MENTIONED) {
  std::stringstream sStream;
  sStream << name_ << " of yet unknown type";
  typeName = sStream.str();
  dataMap[name_] = this;
}

std::map<std::string, Data *> Data::dataMap;

int TypedData::nextId;
