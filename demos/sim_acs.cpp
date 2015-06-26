/* Line Tessellation (LiTe) library
   |||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
   |||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

/******************************************************************************/
/*             LINE BASED TESSELLATION CLASS  - MAIN PROGRAM                   */
/******************************************************************************/

#include "ttessel.h"

int main(int argc, char** argv) {

  int                      opt, seed=5;
  extern char*             optarg;
  Rectangle                w;
  int                      width = 1;        // Width of the square window
  int                      n_i=10;           // Number of loops
  int                      nipl=250;          // Number of iterations per loop
  double                   tau=12;           // Default value for ACS
  double                   fi1=0;            // Default: ACS
  double                   fi2=0;            // Default: ACS
 
  TTessel                  tesl;
  Energy                   eng;

  while((opt=getopt(argc,argv,"t:f:g:s:w:i:j:"))!=EOF) {
    switch(opt) {
    case 't':
      tau = atof(optarg);
      break;
    case 'f':
      fi1 = atof(optarg);
      break; 
    case 'g':
      fi2 = atof(optarg);
      break; 
    case 's':
      seed = atoi(optarg);
      break;
    case 'w':
      width = atoi(optarg);
      break;
    case 'i':
      n_i = atoi(optarg);
      break;
    case 'j':
      nipl = atoi(optarg);
      break;
    default:
      std::cout << "Invalid argument" << std::endl;
      break;
    }
  }
  
  rnd = new CGAL::Random(seed);
  
  tesl.insert_window(Rectangle(Point2(0,0),Point2(width,width)));
  
  SMFChain smf = SMFChain(&eng,0.33,0.33); 
  
  eng.add_theta_segs(-log(tau));
  eng.add_features_segs(is_segment_internal);
  
  eng.add_theta_edges((tau-1)/CGAL_PI);
  eng.add_features_edges(edge_length);
  
  eng.add_features_faces(face_area_2); // sum of squared face areas
  eng.add_theta_faces(tau*tau*tau*tau*fi1);
  
  eng.add_features_faces(face_sum_of_angles); // sum of complementary vertex angles
  eng.add_theta_faces(fi2);
  
  eng.set_ttessel(&tesl);
  
  
  std::cerr << "step" << " " << "energy" << " ";
  std::cerr << "nb_of_int_segments" << " " << "int_length" << std::endl;
  for (int i=0;i!=n_i;i++) {
    smf.step(nipl);
    std::cerr << i << " " << eng.get_value() << " ";
    std::cerr << tesl.number_of_segments()-4 << " ";
    std::cerr << (tesl.get_total_internal_length()-4*width)/2 << std::endl;
  }
 
  std::cerr << std::endl;
  std::cerr << "Tessellation features:" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Number of segments: " << tesl.number_of_segments()-4 << std::endl;
  std::cerr << "Total internal length: " << tesl.get_total_internal_length()-4*width << std::endl; 
  std::cerr << "Sum of squared area of faces: "  << sum_of_faces_squared_areas(&tesl) << std::endl;
  std::cerr << "Angle statistic: "  << sum_of_min_angles(&tesl) << std::endl;
  std::cerr << std::endl;
  std::cout << "Coordinates of segments: " << std::endl;
  tesl.print_all_segments();
  
  
  delete rnd;
  return 0;

}


