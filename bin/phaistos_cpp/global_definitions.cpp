namespace module_camshift {

// Module: energy term initialization
struct EnergyInitialization {


     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyInitialization(const Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
                                       Energy<CHAIN_TYPE> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                                       std::string prefix="") {
     }

     // Constructor - template specific case
     template <typename DBN_TYPE>
     EnergyInitialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                                       Energy<ChainFB> *energy, std::vector<RandomNumberEngine *> *random_number_generators, 
                                       std::string prefix="") {

          Options::OptionValue option;

          // CamShift energy term
          option = options[prefix+"-camshift"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermCamshift::Settings Settings;

               // Add energy term
               energy->add_term(new TermCamshift(chain,
                                                 options.get_settings<Settings>(option,i)));
          }


          // CamShift energy term - cached version
          option = options[prefix+"-camshift-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef typename TermCamshiftCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermCamshiftCached(chain,
                                                       options.get_settings<Settings>(option,i)));
          }

     }

};


}
