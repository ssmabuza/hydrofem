// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_IOBase.hpp"

namespace hydrofem
{

void IOBase::createPVUOutputFile() const
{
  
  std::ofstream outw(getPVUFilename());
  if (!outw.good())
  {
    std::stringstream ss; 
    ss << "Error in \"createPVUOutputFile\"." <<
          std::endl << "Can not open output file '" << getPVUFilename() << "'.";
    throw std::runtime_error(ss.str());
  }
  
  outw << "<?xml version=\"1.0\"?>" << std::endl;
  outw << R"(<VTKFile type="Collection" version="0.1">)" << std::endl;
  outw << "  <Collection>" << std::endl;
  
}

void IOBase::updatePVUOutputFile(const std::string file_name, const double time) const
{
  
  std::ofstream outw(getPVUFilename(),std::ios::app);
  if (!outw.good())
  {
    std::stringstream ss; 
    ss << "Error in \"updatePVUOutputFile\"." <<
          std::endl << "Can not open output file '" << getPVUFilename() << "'.";
    throw std::runtime_error(ss.str());
  }
  outw << "    <DataSet timestep=\"" << std::setprecision(16) << time
  << "\"" << R"( part="0" file=")" << file_name << "\" />" << std::endl;
  
}

void IOBase::finalizePVUOutputFile() const
{
  
  std::ofstream outw(getPVUFilename(),std::ios::app);
  if (!outw.good())
  {
    std::stringstream ss; 
    ss << "Error in \"finalizePVUOutputFile\"." <<
          std::endl << "Can not open output file '" << getPVUFilename() << "'.";
    throw std::runtime_error(ss.str());
  }
  outw << "  </Collection>" << std::endl;
  outw << "</VTKFile>" << std::endl;
  
}


}
// end namespace hydrofem

