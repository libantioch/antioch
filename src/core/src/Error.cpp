#include "antioch/Error.hpp"

namespace Antioch
{

Error::Error(const std::string &type,const std::string &method, const std::string &message)
{
  printError(type,method,message);
  if(isError(type))exit(1);
}

void Error::printError(std::string const & type, std::string const &method, std::string const &message) const
{
  std::cerr << std::endl;
  for(int i = 0; i < 50; i++){std::cerr << "*";}
  std::cerr << std::endl
       << type << " in method " << method << ":" << std::endl
       << message << std::endl;
  for(int i = 0; i < 50; i++){std::cerr << "*";}
  std::cerr << std::endl << std::endl;
}


bool Error::isError(const std::string &type) const
{
  return (type.find("error") != std::string::npos ||
          type.find("Error") != std::string::npos ||
          type.find("ERROR") != std::string::npos);
}

void antiochWarning(const std::string &method, const std::string &message)
{
  Error war("Warning",method,message);
}

void antiochError(const std::string &method, const std::string &message)
{
  Error err("Error",method,message);
}

}
