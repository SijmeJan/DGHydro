#ifndef DG_CONFIGFILE_HPP
#define DG_CONFIGFILE_HPP

#include <map>
#include <string>
#include <iostream>

namespace DGHydro {

  class ConfigFile {
  public:
    ConfigFile(char *fileName, int rank=0);
    ~ConfigFile(void);

    void List();

    template<typename T>
    T GetParameter(const std::string p) {
      if (parameterMap.count(p) != 1) {
        throw std::runtime_error("Parameter not found in map: " + p);
      }
      return (T) parameterMap[p];
    }

  private:
    std::map<std::string, double> parameterMap;

  };

} // namespace DGHydro

#endif  // DG_CONFIGFILE_HPP
