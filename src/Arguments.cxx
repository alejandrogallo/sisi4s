#include <regex>
#include <string>

#include <Arguments.hpp>

namespace sisi4s {

std::string Arguments::normalize_name(std::string const &name) {
  const std::regex symbols("[-_ ]");
  std::string new_name = std::regex_replace(name, symbols, "");
  std::transform(new_name.begin(),
                 new_name.end(),
                 new_name.begin(),
                 [](unsigned char const &c) { return std::tolower(c); });
  return new_name;
}

void Arguments::push(std::string const &name, std::string const &db_index) {
  arguments[normalize_name(name)] = db_index;
}

bool Arguments::present(std::string const &name) {
  return arguments.find(normalize_name(name)) != arguments.end();
}

data::StorePair Arguments::get_data(std::string const &name) {
  const std::string nname = normalize_name(name);
  auto dataIterator(arguments.find(nname));
  if (dataIterator == arguments.end()) {
    std::stringstream sStream;
    sStream << "Missing argument: " << name;
    //    throw new EXCEPTION(std::stringstream() << "Missing argument: " <<
    //    name);
    throw new EXCEPTION(sStream.str());
  }
  std::string db_index = dataIterator->second;
  // Data *data = Data::get(db_index);
  data::Data *data = data::getraw(db_index);
  if (!data) {
    throw new EXCEPTION(_FORMAT("Missing data in db '%s' for key '%s' ",
                                db_index.c_str(),
                                name.c_str()));
  }
  return {db_index, data};
}

std::string Arguments::get_var(std::string const &name) {
  return get_data(normalize_name(name)).second->name;
}

} // namespace sisi4s
