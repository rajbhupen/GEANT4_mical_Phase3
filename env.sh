 


  #Pythia 6                                                                                              
  export PYTHIA6=/products/genie/v3_0_6_sl7/PYTHIA6/pythia6428
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PYTHIA6

  #GSL-2.6                                                                                               
  export GSLHOME=/products/genie/v3_0_6_sl7/GSL26/gsl26
  export GSL_ROOT_DIR=/products/genie/v3_0_6_sl7/GSL26/gsl26
  export GSL_INCLUDE_DIR=$GSLHOME/include
  export GSL_LIBRARY=$GSLHOME/lib
  export PATH=$GSLHOME/bin:$PATH
  export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GSLHOME/lib

  #ROOT 6.20.04                                                                                          
  source /products/genie/v3_0_6_sl7/ROOT6/root-6.20.04/bin/thisroot.sh

  #GENIE                                                                                                 
  export GENIEBASE=/products/genie/v3_0_6_sl7
  export GENIE=$GENIEBASE/Generator-R-3_00_06
  export LOG4CPP_INC=$GENIEBASE/LOG4CPP/log4cpp/include
  export LOG4CPP_LIB=$GENIEBASE/LOG4CPP/log4cpp/lib
  export LHAPDF=$GENIEBASE/LHAPDF6/lhapdf630
  export LHAPATH=$GENIEBASE/LHAPDF6/lhapdf630/share/LHAPDF
  export LHAPDF_INC=$GENIEBASE/LHAPDF6/lhapdf630/include
  export LHAPDF_LIB=$GENIEBASE/LHAPDF6/lhapdf630/lib
  export LD_LIBRARY_PATH=$LHAPDF_LIB:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$LOG4CPP_LIB:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$PYTHIA6:$LD_LIBRARY_PATH
  export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
  export PATH=$GENIE/bin:$ROOTSYS/bin:$GSLHOME/bin:$GENIE/src/stdapp:$PATH
  export LD_LIBRARY_PATH=$GENIE/lib:$GSL_LIB/lib:$LD_LIBRARY_PATH

  #CLHEP                                                                                                 
  export CLHEP_BASE_DIR=/products/GEANT4/sl7/CLHEP/clhep2404
  export PATH=$CLHEP_BASE_DIR/bin:$PATH
  export LD_LIBRARY_PATH=$CLHEP_BASE_DIR/lib:$LD_LIBRARY_PATH
  #GEANT                                                                                                 
  export G4INSTALL=/products/GEANT4/sl7/geant4.10.04.p03-install/share/Geant4-10.4.3/geant4make
  # export G4WORKDIR=/home/user/G4WORK                                                                   
  source $G4INSTALL/geant4make.sh
  source /products/GEANT4/sl7/geant4.10.04.p03-install/bin/geant4.sh


