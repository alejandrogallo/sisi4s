#include <vector>

#include <util/Yaml.hpp>
#include <Sisi4s.hpp>
#include <algorithms/Read.hpp>
#include <algorithms/TensorWriter.hpp>
#include <util/Tensor.hpp>
#include <vendor/filesystem.hpp>

namespace sisi4s {

namespace fs = ghc::filesystem;

template <typename F>
static Tensor<F> *new_tensor_from_dimensions(cc4s::Dimensions const &dims) {
  std::vector<TensorIndex> lens(dims.size()), syms(dims.size(), NS);
  std::transform(dims.begin(),
                 dims.end(),
                 lens.begin(),
                 [](cc4s::Dimension const &d) { return d.length; });

  return new Tensor<F>(dims.size(), syms.data(), lens.data(), *Sisi4s::world);
}

IMPLEMENT_EMPTY_DRYRUN(Read) {}

IMPLEMENT_ALGORITHM(Read) {

  const std::string fileName = getTextArgument("fileName");

  const auto dataPath =
      fs::path(fileName).replace_extension(fs::path("elements"));

  const auto node = YAML::LoadFile(fileName);
  cc4s::ReadHeader header = node.as<cc4s::ReadHeader>();
  LOG(0, "Read") << "version: " << header.version << "\n";
  for (auto const &d : header.dimensions) {
    LOG(0, "Read:Header") << "#dimen length: " << d.length << "\n";
    LOG(0, "Read:Header") << "#dimen type  : " << d.type << "\n";
  }

  size_t count = std::accumulate(
      header.dimensions.begin(),
      header.dimensions.end(),
      1UL,
      [](size_t a, cc4s::Dimension const &b) { return a * b.length; });

  LOG(0, "Read:Header") << "#elems : " << count << "\n";

  switch (header.scalarType) {
  case cc4s::ScalarType::Real64: {
    LOG(0, "Read") << "Real64 tensor"
                   << "\n";
    using F = typename cc4s::ScalarTypeTraits<cc4s::ScalarType::Real64>::type;
    auto t = new_tensor_from_dimensions<F>(header.dimensions);

    switch (header.elementsType) {
    case cc4s::ElementFileType::TextFile: {
      LOG(0, "Read") << "Text tensor"
                     << "\n";
      if (Sisi4s::world->rank == 0) {
        std::vector<F> values(count);
        size_t idx = 0;
        {
          std::ifstream s(dataPath.c_str());
          for (std::string line; std::getline(s, line);)
            values[idx++] = std::atof(line.c_str());
        }

        {
          // write to tensor
          std::vector<TensorIndex> indices(count);
          std::iota(indices.begin(), indices.end(), 0);
          t->write(count, indices.data(), values.data());
        }
      }
      break;
    }
    case cc4s::ElementFileType::IeeeBinaryFile:
      LOG(0, "Read") << "Binary tensor"
                     << "\n";
      t->read_dense_from_file(dataPath.c_str());
      break;
    }

    allocatedTensorArgument<F>("destination", t);
    break;
  }
  case cc4s::ScalarType::Complex64: {
    LOG(0, "Read") << "Complex64 tensor"
                   << "\n";
    using F =
        typename cc4s::ScalarTypeTraits<cc4s::ScalarType::Complex64>::type;
    auto t = new_tensor_from_dimensions<F>(header.dimensions);
    t->read_dense_from_file(dataPath.c_str());
    allocatedTensorArgument<F>("destination", t);
    break;
  }
  }
}

} // namespace sisi4s

// This is the uninteresting YAML details

namespace YAML {

using namespace sisi4s::cc4s;

template <>
struct convert<AxisType> {

  static Node encode(AxisType const &i) {
    Node n;
    switch (i) {
    case AxisType::AuxiliaryField: n = "AuxiliaryField"; break;
    case AxisType::State: n = "State"; break;
    }
    return n;
  }

  static bool decode(Node const &n, AxisType &r) {
    std::string _n = n.as<std::string>();
    if (!n.IsScalar()) return false;
    if (_n == "AuxiliaryField") {
      r = AxisType::AuxiliaryField;
      return true;
    } else if (_n == "State") {
      r = AxisType::State;
      return true;
    }
    return false;
  }
};

template <>
struct convert<ElementFileType> {
  using Native = ElementFileType;

  static Node encode(Native const &i) {
    Node n;
    switch (i) {
    case Native::TextFile: n = "TextFile"; break;
    case Native::IeeeBinaryFile: n = "IeeeBinaryFile"; break;
    }
    return n;
  }

  static bool decode(Node const &n, Native &r) {
    std::string _n = n.as<std::string>();
    if (!n.IsScalar()) return false;
    if (_n == "TextFile") {
      r = Native::TextFile;
      return true;
    } else if (_n == "IeeeBinaryFile") {
      r = Native::IeeeBinaryFile;
      return true;
    }
    return false;
  }
};

template <>
struct convert<ScalarType> {
  using Native = ScalarType;

  static Node encode(Native const &i) {
    Node n;
    switch (i) {
    case Native::Real64: n = "Real64"; break;
    case Native::Complex64: n = "Complex64"; break;
    }
    return n;
  }

  static bool decode(Node const &n, Native &r) {
    std::string _n = n.as<std::string>();
    if (!n.IsScalar()) return false;
    if (_n == "Real64") {
      r = Native::Real64;
      return true;
    } else if (_n == "Complex64") {
      r = Native::Complex64;
      return true;
    }
    return false;
  }
};

template <>
struct convert<ReadHeader::Version> {
  using Native = ReadHeader::Version;

  static Node encode(Native const &i) {
    Node n;
    switch (i) {
    case Native::ONE: n = "100"; break;
    }
    return n;
  }

  static bool decode(Node const &n, Native &r) {
    if (!n.IsScalar()) return false;
    std::string _n = n.as<std::string>();
    if (_n == "100") {
      r = Native::ONE;
      return true;
    }
    return false;
  }
};

template <>
struct convert<ReadableType> {
  using Native = ReadableType;

  static Node encode(Native const &i) {
    Node n;
    switch (i) {
    case Native::Tensor: n = "Tensor"; break;
    }
    return n;
  }

  static bool decode(Node const &n, Native &r) {
    if (!n.IsScalar()) return false;
    std::string _n = n.as<std::string>();
    if (_n == "Tensor") {
      r = Native::Tensor;
      return true;
    }
    return false;
  }
};

template <>
struct convert<Dimension> {
  using Native = Dimension;

  static Node encode(Native const &i) {
    Node n;
    n["length"] = i.length;
    n["type"] = i.type;
    return n;
  }

  static bool decode(Node const &n, Native &r) {
    if (!n.IsMap()) return false;
    YAML_ASSERT_KEY(n, "length");
    YAML_ASSERT_KEY(n, "type");
    r.length = n["length"].as<size_t>();
    r.type = n["type"].as<AxisType>();
    return true;
  }
};

template <>
struct convert<ReadHeader> {
  using Native = ReadHeader;

  static Node encode(Native const &i) {
    Node n;
    n["version"] = i.version;
    n["type"] = i.type;
    n["dimensions"] = i.dimensions;
    n["scalarType"] = i.scalarType;
    n["elements"]["type"] = i.elementsType;
    n["unit"] = i.unit;
    return n;
  }

  static bool decode(Node const &n, Native &r) {
    if (!n.IsMap()) return false;
    YAML_ASSERT_KEY(n, "dimensions");
    YAML_ASSERT_KEY(n, "elements");
    YAML_ASSERT_KEY(n, "scalarType");
    YAML_ASSERT_KEY(n, "type");
    YAML_ASSERT_KEY(n, "unit");
    YAML_ASSERT_KEY(n, "version");
    r.version = n["version"].as<ReadHeader::Version>();
    r.type = n["type"].as<ReadableType>();
    r.dimensions = n["dimensions"].as<Dimensions>();
    r.scalarType = n["scalarType"].as<ScalarType>();
    r.elementsType = n["elements"]["type"].as<ElementFileType>();
    r.unit = n["unit"].as<double>();
    return true;
  };
};

} // namespace YAML

namespace sisi4s {

IMPLEMENT_EMPTY_DRYRUN(Write) {}

IMPLEMENT_ALGORITHM(Write) {

  using namespace sisi4s::cc4s;
  const bool binary_p = getBooleanArgument("binary", true);

  ReadHeader header;
  header.version = ReadHeader::Version::ONE;
  header.type = ReadableType::Tensor;
  header.elementsType =
      binary_p ? ElementFileType::IeeeBinaryFile : ElementFileType::TextFile;
  header.unit = 1.0;

  std::string

      source_name = "source",
      dataName = getArgumentData(source_name)->getName(),
      fileNameInput = getTextArgument("fileName", dataName);

  const auto dataPath = fs::path(fileNameInput)
                            .replace_extension(fs::path("elements")),
             yamlPath =
                 fs::path(fileNameInput).replace_extension(fs::path("yaml"));

  if (typeid(double) == TensorWriter::check_type(getArgumentData("source"))) {
    LOG(1, "TensorWriter") << "Writing real tensor" << std::endl;
    auto t = getTensorArgument<Float64>(source_name);
    for (size_t i = 0; i < t->order; i++)
      header.dimensions.push_back({t->lens[i], AxisType::State});
    TensorWriter::write<Float64>(dataName, dataPath, t, binary_p, "", "", "");
    header.scalarType = ScalarType::Real64;
  } else {
    LOG(1, "TensorWriter") << "Writing complex tensor" << std::endl;
    auto t = getTensorArgument<sisi4s::Complex64>(source_name);
    for (size_t i = 0; i < t->order; i++)
      header.dimensions.push_back({t->lens[i], AxisType::State});
    TensorWriter::write<sisi4s::complex>(dataName,
                                         dataPath,
                                         t,
                                         binary_p,
                                         "",
                                         "",
                                         "");
    header.scalarType = ScalarType::Complex64;
  }

  YAML::Emitter yout;
  yout << YAML::convert<ReadHeader>::encode(header);
  // YAML::Node node = header.scalarType;
  // yout << YAML::convert<ReadHeader>(header);
  // yout << header;
  std::ofstream fout(yamlPath);
  fout << yout.c_str();
}
} // namespace sisi4s
