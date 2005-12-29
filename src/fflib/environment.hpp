#include "config-wrapper.h"
#include <string>
#include <list>
#include <map>
typedef std::list<std::string>  OneEnvironmentData;
typedef std::map<std::string,OneEnvironmentData > EnvironmentData;

extern EnvironmentData  environment;
extern long verbosity;

bool EnvironmentInsert(std::string key,std::string item,std::string before);

void GetEnvironment();