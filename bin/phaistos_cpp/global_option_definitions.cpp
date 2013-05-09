namespace module_camshift {

// Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {                    
     }

     // Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {

          // Import namespace for make_vector
          using namespace boost::fusion;

// Define defaults for the different modes
                     ModeDefinitions mode_definitions(target, chain);


          // CamShift energy term
          for (int counter = occurrences[prefix+"-camshift"]; counter > 0; counter--) {

               // typedef typename TermCamshift::Settings Settings;
               // boost::shared_ptr<Settings> settings(new Settings());
               typedef TermCamshift EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Camshift chemical shift energy (" + prefix + ")",
                         prefix+"-camshift", settings,
                         make_vector(
                              make_vector(std::string("star-filename"),
                                          std::string("Chemical shifts file from RefDB"),
                                          &settings->star_filename),
                              make_vector(std::string("energy-type"),
                                          std::string("Type of energy expression: 1=Robustelli 2=Jeffrey's Prior 3=Naive Gaussian 4=ProCS"),
                                          &settings->energy_type),
                              make_vector(std::string("flat-bottom-tolerance"),
                                          std::string("The tolerance parameter for flat bottom potentials"),
                                          &settings->flat_bottom_tolerance)
                              )),
                    super_group, counter==1);
          }

          // CamShift energy term
          for (int counter = occurrences[prefix+"-camshift-cached"]; counter > 0; counter--) {

               //typedef typename TermCamshiftCached::Settings Settings;
               // booost::shared_ptr<Settings> settings(new Settings());
               typedef TermCamshiftCached EnergyTerm;
                typedef EnergyTerm::Settings Settings;
                 boost::shared_ptr<Settings> settings(
                               SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));


               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Camshift chemical shift energy - cached version (" + prefix + ")",
                         prefix+"-camshift-cached", settings,
                         make_vector(
                              make_vector(std::string("star-filename"),
                                          std::string("Chemical shifts file from RefDB"),
                                          &settings->star_filename),
                              make_vector(std::string("energy-type"),
                                          std::string("Type of energy expression: 1=Robustelli 2=Jeffrey's Prior 3=Naive Gaussian 4=ProCS"),
                                          &settings->energy_type),
                              make_vector(std::string("cutoff-distance"),
                                          std::string("Distance beyond which contributions are set to 0. Use only for debugging!"),
                                          reinterpret_cast<ProgramOptionParser::WrappedDoublePointer*>(&settings->cutoff_distance)),
                              make_vector(std::string("flat-bottom-tolerance"),
                                          std::string("The tolerance parameter for flat bottom potentials"),
                                          &settings->flat_bottom_tolerance)
                              )),
                    super_group, counter==1);
          }
    }
};

}
