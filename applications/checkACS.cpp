/* Line Tessellation (LiTe) library
   |||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu. 
   |||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */


# include "ttessel.h"
# include <fstream>
# include <sstream>
# include <string>
 
  



using std::ofstream;
using std::string;
using std::ostringstream;

 int main(int argc, char** argv) {
  
 int                      opt, seed=5;
 extern char*             optarg;
 int                      width = 1;       // Width of the square window
 int                      n_i =100;        // Number of loops
 int                      nipl=1;          // Number of iterations per loop
 double                   tau=6;           // Scale parameter in ACS model
 int                      prec=16;          // Precision of length output
 
 TTessel                  tesl;
 Energy                   eng;
 ModifCounts              m_stat;

 std::vector<Segment>     inter_segs;      
 Segment                  s;
 CGAL::Object             inter;
 Point2                   p; 
 
 ofstream                 inter_out; 
 ofstream                 stat_out; 
 ofstream                 summary_out;
 ofstream                 tessel_out;

 char inittessel[]="t00.txt";
 ostringstream case_number;
 ostringstream statfile;
 ostringstream interfile;
 ostringstream tesselfile;

 while((opt=getopt(argc,argv,"t:s:w:i:j:"))!=EOF) {
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

 


 
  /* Preparing output files */

  case_number << "_" << tau << "_" << seed << "_" << width << "_" <<  n_i << "_" << nipl << "_" << inittessel;
  statfile   << "stat"   << case_number.str() ;
  interfile  << "inter"  << case_number.str() ;
  tesselfile << "tessel" << case_number.str() ;


  stat_out.open(statfile.str().data(),std::ios::out); 
  inter_out.open(interfile.str().data(),std::ios::out); 
  tessel_out.open(tesselfile.str().data(),std::ios::out); 
  summary_out.open("summary.txt",std::ios::app);

  /* Writing sample informations */

  
  summary_out   << "init="        << inittessel       << " "
                << "tau="         << tau              << " "
                << "seed="        << seed             << " "
                << "n_sim="       << n_i              << " "
	        << "step="        << nipl             << " "
                << "stat="        << statfile.str()   << " "
                << "inter="       << interfile.str()  << " "
               << "last_tessel=" << tesselfile.str()
               << std::endl;
 
  
  rnd = new CGAL::Random(seed);
  tesl.insert_window(Rectangle(Point2(0,0),Point2(width,width)));
  SMFChain smf = SMFChain(&eng,0.33,0.33); 
 

  /*  Energy of ACS model  */

  eng.add_features_segs(is_segment_internal);
  eng.add_features_edges(edge_length);

  eng.add_theta_segs(-log(tau));
  eng.add_theta_edges((tau-1)/CGAL_PI);
  
  eng.set_ttessel(&tesl);


  /* Initializing segments intersecting tessellation */
  
  
  for (TTessel::Edge_iterator e=tesl.edges_begin();e!=tesl.edges_end();e++){
  inter_segs.push_back(Segment(e->source()->point(),e->target()->point()));
  }

  inter_segs.push_back(Segment(Point2(0,width/2),Point2(width,width/2)));
  inter_segs.push_back(Segment(Point2(width/2,0),Point2(width/2,width)));
 
 

   /* Header of statfile */

   stat_out  << "i"        << " "  // Sample id
            << "l(T)"     << " "   // Total internal length
            << "n(T)"     << " "   // Number of internal segments
            << "n_bl(T)"  << " "   // Number of blocking segments
            << "n_nbl(T)" << " "   // Number of non-blocking segments 
            << "m(T)"     << " "   // Number of internal vertices
            << "prop_S"   << " "   // Number of proposed splits
            << "acc_S"    << " "   // Number of accepted splits
	    << "prop_M"   << " "   // Number of proposed merges 
            << "acc_M"    << " "   // Number of accepted merges
	    << "prop_F"   << " "   // Number of proposed flips 
            << "acc_F"             // Number of accepted flips
            << std::endl;
   stat_out << std::setprecision(prec);

  /* Header of interfile */

   inter_out <<  "i" << " " << "seg"<< " " << "x" << " " << "y"  << std::endl;
  
 /* Sampling */
 
  for (int i=0;i!=n_i;i++) { // loop over a number of samples
	m_stat=smf.step(nipl); 
        /* Writing sample statistics */
          stat_out << i << " "
                  <<((tesl.get_total_internal_length()-4*width)/2) << " "
                  << (tesl.number_of_segments()-4) << " "
                  << tesl. number_of_blocking_segments() << " "
                  << tesl. number_of_non_blocking_segments() << " "
                  << number_of_internal_vertices(tesl) << " "
                  << m_stat.proposed_S << " "
                  << m_stat.accepted_S << " "
                  << m_stat.proposed_M << " "
                  << m_stat.accepted_M << " "
                  << m_stat.proposed_F << " "
                  << m_stat.accepted_F 
                  << std::endl;

	   /* Writing intersections */
	 
	 
          for (TTessel::Edge_iterator e=tesl.edges_begin();e!=tesl.edges_end();e++){
  	      if (!e->face()->is_unbounded() && !e->twin()->face()->is_unbounded()) {	
  		 for (unsigned int k=0;k!=inter_segs.size();k++){
  		      s=inter_segs[k];
  		      inter = CGAL::intersection(s,Segment(e->source()->point(), e->target()->point()));
  		      if (CGAL::assign(p,inter)) 
			inter_out <<  i << " " << k << " " << p << std::endl;
  		 }
  	      }
  	   }
	  
  }   
 
  // Save final tessellation into a data file
  tesl.write(tessel_out);
  tessel_out.close();
  
return 0;
}


