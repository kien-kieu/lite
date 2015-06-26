/* Line Tessellation (LiTe) library
   |||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
   |||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

/******************************************************************************/
/*   PSEUDO-LIKELIHOOD INFERENCE BASED ON STOCHASTIC APPROXIMATION            */
/*   SIMULATION ON CRTT MODEL                                                 */
/******************************************************************************/

#include "ttessel.h"

int main(int argc, char** argv) {

  int                      opt, seed=5;
  extern char*             optarg;
  Rectangle                w;
  int                      width = 1;        // Width of the square window
  int                      n_burnin=5000;    // Number of iterations for burnin
  double                   tau=2;            // Model parameter tau
  double                   e0=5e-2;          // Initial value for epsilon
  int                      n_sa=20;          // Number of loops for stochastic 
                                             // approximation
 
  TTessel                  obsTes;
  Energy                   trueMod, fitMod;

  while((opt=getopt(argc,argv,"t:s:w:b:c:n:"))!=EOF) {
    switch(opt) {
    case 't':
      tau = atof(optarg);
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 'w':
      width = atoi(optarg);
      break;
    case 'b':
      n_burnin = atoi(optarg);
      break;
    case 'c':
      e0 = atof(optarg);
      break;
    case 'n':
      n_sa = atoi(optarg);
      break;
    default:
      std::cout << "Invalid argument" << std::endl;
      break;
    }
  }
  
  rnd = new CGAL::Random(seed);
  
  // Start generating the "observed tessellation"
  obsTes.insert_window(Rectangle(Point2(0,0),Point2(width,width)));

  SMFChain smf = SMFChain(&trueMod,0.33,0.33); 
  
  trueMod.add_theta_segs(log(tau));
  trueMod.add_features_segs(minus_is_segment_internal);
  
  trueMod.set_ttessel(&obsTes);
  
  std::cout << "True parameter value " << trueMod.get_theta().segs[0] << std::endl;  
  smf.step(n_burnin);
  std::cout << std::endl;
  std::cout << "Some features of the \"observed\" tessellation:" << std::endl;
  std::cout << "\tnumber of internal segments ";
  std::cout << obsTes.number_of_segments()-4 << std::endl;
  std::cout << "\tnumber of internal non-blocking segments ";
  std::cout << obsTes.number_of_non_blocking_segments() << std::endl;
  std::cout << "\tnumber of internal blocking segments ";
  std::cout << obsTes.number_of_blocking_segments() << std::endl;
  std::cout << "\ttotal perimeter " << obsTes.get_total_internal_length();
  std::cout << std::endl;
  std::cout << "Exact pseudo-likelihood estimate ";
  std::cout << log(obsTes.number_of_non_blocking_segments()*CGAL_PI/obsTes.get_total_internal_length()) << std::endl << std::endl;
  // Generation of observed tessellation: done
  // Create fitting object
  fitMod.add_theta_segs(0.0);
  fitMod.add_features_segs(minus_is_segment_internal);
  fitMod.set_ttessel(&obsTes);
  PseudoLikStochApprox plsa(&fitMod,50.0,e0,0.67);
  for (int i=0;i!=n_sa;i++) {
    plsa.step();
    std::cout << "parameter estimate: ";
    std::cout << plsa.GetEnergy()->get_theta().segs[0] << std::endl;
  }
  delete rnd;
  return 0;
}
