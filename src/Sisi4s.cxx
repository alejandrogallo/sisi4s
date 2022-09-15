/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Sisi4s.hpp>
#include <Parser.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/Timer.hpp>
#include <DryTensor.hpp>
#include <util/FlopsCounter.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/Emitter.hpp>
#include <util/Exception.hpp>
#include <fstream>
#include <string>
#include <sstream>

// TODO: to be removed from the main class
#include <math/MathFunctions.hpp>

using namespace sisi4s;
using namespace CTF;

Sisi4s::Sisi4s() {
}

Sisi4s::~Sisi4s() {
}

void Sisi4s::run() {
  EMIT() << YAML::BeginMap;
  printBanner();
  listHosts();
  std::vector<Algorithm *> algorithms;
  if (options->cc4s) {
    LOG(0, "root")
      << "WARNING: "
      << "You are using the old cc4s input language, which is deprecated\n"
      << "Consider converting your input file to the yaml DSL using the"
      << " cc4s-to-yaml.py tool found in your package source."
      << "\n";
    InputFileParser<InputFileFormat::HUMMEL> parser(options->inFile);
    algorithms = parser.parse();
  } else {
    InputFileParser<InputFileFormat::YAML> parser(options->inFile);
    algorithms = parser.parse();
  }
  LOG(0, "root") <<
    "execution plan read, steps=" << algorithms.size() << std::endl;
  EMIT() <<
    YAML::Key << "execution-plan-size" << YAML::Value << algorithms.size();

  EMIT() << YAML::Key << "steps" << YAML::Value << YAML::BeginSeq;

  int64_t rootFlops, totalFlops;
  Time totalTime;
  {
    FlopsCounter rootCounter(&rootFlops);
    FlopsCounter totalCounter(&totalFlops, world->comm);
    Timer totalTimer(&totalTime);

    for (unsigned int i(0); i < algorithms.size(); ++i) {
      EMIT() << YAML::BeginMap;
      LOG(0, "root") << "step=" << (i+1) << ", " << algorithms[i]->getName() << std::endl;
      EMIT() << YAML::Key << "step" << YAML::Value << (i+1)
        << YAML::Key << "name" << YAML::Value << algorithms[i]->getName();

      int64_t flops;
      Time time;
      {
        FlopsCounter flopsCounter(&flops);
        Timer timer(&time);
        algorithms[i]->run();
        delete algorithms[i];
      }

      std::stringstream realtime;
      realtime << time;
      LOG(1, "root") << "step=" << (i+1) << ", realtime=" << realtime.str() << " s"
        << ", operations=" << flops / 1e9 << " GFLOPS/core"
        << ", speed=" << flops / 1e9 / time.getFractionalSeconds() << " GFLOPS/s/core" << std::endl;
      EMIT() << YAML::Key << "realtime" << YAML::Value << realtime.str()
        << YAML::Comment(" seconds")
        << YAML::Key << "floating-point-operations" << YAML::Value << flops
        << YAML::Comment("on root process")
        << YAML::Key << "flops" << YAML::Value << flops / time.getFractionalSeconds();
      printStatistics();
      EMIT() << YAML::EndMap;
    }
  }

  EMIT() << YAML::EndSeq;

  OUT() << std::endl;
  std::stringstream totalRealtime;
  totalRealtime << totalTime;
  LOG(0, "root") << "total realtime=" << totalRealtime.str() << " s" << std::endl;
  LOG(0, "root") << "total operations=" << rootFlops / 1e9 << " GFLOPS/core"
    << " speed=" << rootFlops/1e9 / totalTime.getFractionalSeconds() << " GFLOPS/s/core" << std::endl;
  LOG(0, "root") << "overall operations=" << totalFlops / 1.e9 << " GFLOPS"
    << std::endl;
  EMIT() << YAML::Key << "realtime" << YAML::Value << totalRealtime.str()
    << YAML::Key << "floating-point-operations" << YAML::Value << rootFlops
    << YAML::Comment("on root process")
    << YAML::Key << "flops" << YAML::Value << rootFlops / totalTime.getFractionalSeconds()
    << YAML::Key << "total-floating-point-operations" << totalFlops
    << YAML::Comment("of all processes");

  EMIT() << YAML::EndMap;
}

void Sisi4s::dryRun() {
  EMIT() << YAML::BeginMap;
  printBanner();
  LOG(0, "root") <<
    "DRY RUN - nothing will be calculated" << std::endl;
  OUT() << std::endl;
  InputFileParser<InputFileFormat::YAML> parser(options->inFile);
  std::vector<Algorithm *> algorithms(parser.parse());
  LOG(0, "root") <<
    "execution plan read, steps=" << algorithms.size() << std::endl;
  EMIT() <<
    YAML::Key << "execution-plan-size" << YAML::Value << algorithms.size();

  EMIT() << YAML::Key << "steps" << YAML::Value << YAML::BeginSeq;

  for (unsigned int i(0); i < algorithms.size(); ++i) {
    EMIT() << YAML::BeginMap;
    LOG(0, "root") << "step=" << (i+1) << ", " << algorithms[i]->getName() << std::endl;
    EMIT() << YAML::Key << "step" << YAML::Value << (i+1)
      << YAML::Key << "name" << YAML::Value << algorithms[i]->getName();
    algorithms[i]->dryRun();
    LOG(0, "root")
      << "estimated memory=" << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
      << " GB" << std::endl;
    EMIT()
      << YAML::Key << "estimated-total-memory"
      << YAML::Value << DryMemory::maxTotalSize / (1024.0*1024.0*1024.0)
      << YAML::Comment("GB");
    EMIT() << YAML::EndMap;
  }
  EMIT() << YAML::EndSeq;
  EMIT() << YAML::EndMap;
}


void Sisi4s::printBanner() {
  std::stringstream buildDate;
  buildDate << __DATE__ << " " << __TIME__;

  OUT() <<
    (" ____    __ __    ____     ____      " "\n"
     "/\\  _`\\ /\\ \\\\ \\  /\\  _`\\  /\\  _`\\    " "\n"
     "\\ \\,\\L\\_\\ \\ \\\\ \\ \\ \\ \\/\\_\\\\ \\ \\/\\_\\  " "\n"
     " \\/_\\__ \\\\ \\ \\\\ \\_\\ \\ \\/_/_\\ \\ \\/_/_ " "\n"
     "   /\\ \\L\\ \\ \\__ ,__\\ \\ \\L\\ \\\\ \\ \\L\\ \\" "\n"
     "   \\ `\\____\\/_/\\_\\_/\\ \\____/ \\ \\____/" "\n"
     "    \\/_____/  \\/_/   \\/___/   \\/___/ " "\n"
     "\n");

  LOG(0, "root") << "version=" << SISI4S_VERSION <<
    ", date=" << SISI4S_DATE << std::endl;
  LOG(0, "root") << "build date=" << buildDate.str() << std::endl;
  LOG(0, "root") << "compiler=" << COMPILER_VERSION << std::endl;
  LOG(0, "root") << "total processes=" << Sisi4s::world->np << std::endl;
  OUT() << std::endl;

  EMIT()
    << YAML::Key << "version" << YAML::Value << SISI4S_VERSION
    << YAML::Key << "build-date" << YAML::Value << buildDate.str()
    << YAML::Key << "compiler" << YAML::Value << COMPILER_VERSION
    << YAML::Key << "total-processes" << Sisi4s::world->np;
}

void Sisi4s::printStatistics() {
  std::string fieldName;
  int64_t peakVirtualSize, peakPhysicalSize;
  // assuming LINUX
  std::ifstream statusStream("/proc/self/status", std::ios_base::in);
  std::string line;
  while (std::getline(statusStream, line)) {
    std::istringstream lineStream(line);
    lineStream >> fieldName;
    if (fieldName == "VmPeak:") {
      lineStream >> peakVirtualSize;
    } else if (fieldName == "VmHWM:") {
      lineStream >> peakPhysicalSize;
    }
    // TODO: check memory unit, currently assumed to be kB
  }
  statusStream.close();
  real unitsPerGB(1024.0*1024.0);
  LOG(0, "root") << "peak physical memory=" << peakPhysicalSize / unitsPerGB << " GB/core"
    << ", peak virtual memory: " << peakVirtualSize / unitsPerGB << " GB/core" << std::endl;
  EMIT()
    << YAML::Key << "peak-physical-memory" << YAML::Value << peakPhysicalSize / unitsPerGB
    << YAML::Comment("GB/core")
    << YAML::Key << "peak-virtual-memory" << YAML::Value <<  peakVirtualSize / unitsPerGB
    << YAML::Comment("GB/core");

  int64_t globalPeakVirtualSize, globalPeakPhysicalSize;
  MpiCommunicator communicator(world->rank, world->np, world->comm);
  communicator.reduce(peakPhysicalSize, globalPeakPhysicalSize);
  communicator.reduce(peakVirtualSize, globalPeakVirtualSize);
  LOG(0, "root") << "overall peak physical memory="
    << globalPeakPhysicalSize / unitsPerGB << " GB"
    << ", overall virtual memory=" << globalPeakVirtualSize / unitsPerGB << " GB" << std::endl;
  EMIT()
    << YAML::Key << "total-peak-physical-memory" << YAML::Value << globalPeakPhysicalSize / unitsPerGB
    << YAML::Comment("GB")
    << YAML::Key << "peak-virtual-memory" << YAML::Value <<  globalPeakVirtualSize / unitsPerGB
    << YAML::Comment("GB");
}

bool Sisi4s::isDebugged() {
  // assuming LINUX
  std::ifstream statusStream("/proc/self/status", std::ios_base::in);
  std::string line;
  while (std::getline(statusStream, line)) {
    std::string pidField("TracerPid:");
    size_t position(line.find(pidField));
    if (position != std::string::npos) {
      std::stringstream pidStream(line.substr(position + pidField.length()));
      size_t pid; pidStream >> pid;
      if (pid > 0) {
        LOG(0, "root") << "Debugger present" << std::endl;
      }
      return pid > 0;
    }
  }
  return false;
}

void Sisi4s::listHosts() {
  char ownName[MPI_MAX_PROCESSOR_NAME];
  int nameLength;
  MPI_Get_processor_name(ownName, &nameLength);
  ownName[nameLength] = 0;

  if (world->rank == 0) {
    std::map<std::string,std::vector<int>> ranksOfHosts;
    // enter hostname/rank
    ranksOfHosts[ownName].push_back(0);
    // receive all names but own name from remote ranks
    for (int remoteRank(1); remoteRank < world->np; ++remoteRank) {
      char remoteName[MPI_MAX_PROCESSOR_NAME];
      MPI_Recv(
        remoteName, MPI_MAX_PROCESSOR_NAME, MPI_BYTE,
        remoteRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE
      );
      // and enter remote name/rank into map
      ranksOfHosts[remoteName].push_back(remoteRank);
    }
    EMIT() << YAML::Key << "hosts" << YAML::Value;
    EMIT() << YAML::BeginSeq;
    for (auto &ranksOfHost: ranksOfHosts) {
      EMIT() << YAML::BeginMap
        << YAML::Key << "host" << YAML::Value << ranksOfHost.first
        << YAML::Key << "ranks" << YAML::Value;
      EMIT() << YAML::Flow << YAML::BeginSeq;
      for (auto &rank: ranksOfHost.second) {
        EMIT() << rank;
      }
      EMIT() << YAML::EndSeq;
      EMIT() << YAML::EndMap;
    }
    EMIT() << YAML::EndSeq;
  } else {
    // send own name 
    MPI_Send(
      ownName, MPI_MAX_PROCESSOR_NAME, MPI_BYTE,
      0, 0, MPI_COMM_WORLD
    );
  }
}


World *Sisi4s::world;
Options *Sisi4s::options;


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  Sisi4s::world = new World(argumentCount, arguments);
  Sisi4s::options = new Options(argumentCount, arguments);
  Sisi4s::options->parse();
  Log::setRank(Sisi4s::world->rank);
  Log::setFileName(Sisi4s::options->logFile);
  Log::setLogLevel(Sisi4s::options->logLevel);
  Emitter::setFileName(Sisi4s::options->yamlOutFile);
  Emitter::setRank(Sisi4s::world->rank);

  Sisi4s sisi4s;
  if (Sisi4s::isDebugged()) {
    // run without try-catch in debugger to allow tracing throwing code
    if (Sisi4s::options->dryRun) sisi4s.dryRun();
    else sisi4s.run();
  } else {
    // without debugger: catch and write exception cause
    try {
      if (Sisi4s::options->dryRun) sisi4s.dryRun();
      else sisi4s.run();
    } catch (const std::string &msg) {
      LOG(0) << std::endl << msg << std::endl;
      throw msg;
    } catch (const char *msg) {
      LOG(0) << std::endl << msg << std::endl;
      throw msg;
    } catch (DetailedException *cause) {
      LOG(0) << std::endl << cause->getMessage() << std::endl;
    }
  }

  EMIT_FLUSH();
  MPI_Finalize();
  return 0;
}

