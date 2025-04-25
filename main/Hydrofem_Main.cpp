// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_DriverFactory.hpp"

using namespace Hydrofem;

int main(int argc, char** argv)
{
  Eigen::initParallel();  
  // program options handler object
  try {
    
    std::string configfilename;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v", "print version string")
      ("help", "produce help message")
      ("config", po::value<std::string>(&configfilename)->default_value("default.config"), "Configuration or run script");
  
    // Declare a group of options that will be 
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()("input-file", po::value<std::vector<std::string>>(), "input file");
    {
      po::variables_map vm_;
      po::options_description cmdline_options_;
      cmdline_options_.add(generic).add(config).add(hidden);
      po::options_description config_file_options_;
      config_file_options_.add(config).add(hidden);
      po::options_description visible_("Allowed options");
      visible_.add(generic).add(config);
      po::parsed_options parsed_ = po::command_line_parser(argc, argv).options(cmdline_options_).allow_unregistered().run();
      po::store(parsed_,vm_);
      po::notify(vm_);
    }
    {
      po::variables_map vm_;
      po::options_description cmdline_options_;
      cmdline_options_.add(generic).add(config).add(hidden);
      po::options_description config_file_options_;
      config_file_options_.add(config).add(hidden);
      po::options_description visible_("Allowed options");
      visible_.add(generic).add(config);
      po::parsed_options parsed_ =po::command_line_parser(argc, argv).options(cmdline_options_).allow_unregistered().run();
      po::store(parsed_,vm_);
      std::ifstream ifs_(configfilename.c_str());
      if (ifs_.fail())
      {
        std::cout << "Could not open configfile " << configfilename << "\nAborting.\n";
        exit(1);
      }
      po::store(po::parse_config_file(ifs_, config_file_options_,true), vm_);
      po::notify(vm_);
    }
     
    // build option parser 
    auto option_handler = std::make_shared<OptionHandler>(argc,argv,generic,hidden,configfilename);
    // build driver factory 
    auto driver_factory = std::make_shared<DriverFactory>(option_handler);
    // build driver
    std::shared_ptr<Driver> driver = driver_factory->buildDriver();
    // setup the driver
    driver->setup();
    // solve the problem and write solution to file
    driver->solve();
    
  } catch (std::exception &e) {
    
    std::cout << "Exception: " << e.what() << std::endl;
    return 1;
    
  } catch (...) {
    
    std::cout << "Other exceptions..." << std::endl;
    
  }
  
  std::cout << "All simulations completed." << std::endl;
  return 0;
}
