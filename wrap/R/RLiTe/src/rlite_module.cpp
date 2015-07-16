/* Line Tessellation (LiTe) library
   Release 1.1
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
   Copyright INRA 2006-2015.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

//include <RcppCommon.h>
#include <iostream>
#include <fstream>
#define STRICT_R_HEADERS /* otherwise conflict between R_ext/Constants.h and 
			    CGAL/Nef_polyhedron_2.h about PI. */
#include <Rcpp.h>
#include "ttessel.h"
namespace Rcpp {
  template<> SEXP wrap(const CatVector&);
  template<> CatVector as(SEXP);
  template<> SEXP wrap(const Polygons&);
}
RCPP_EXPOSED_CLASS(LineTes)
RCPP_EXPOSED_CLASS(TTessel)
RCPP_EXPOSED_CLASS(Energy)
RCPP_EXPOSED_CLASS(SMFChain)
RCPP_EXPOSED_CLASS(PseudoLikDiscrete)
RCPP_EXPOSED_CLASS(PLInferenceNOIS)

void setDomain(LineTes *tes, double width, double height) {
  tes->clear();
  tes->insert_window(Rectangle(Point2(0,0),Point2(width,height)));
}
void setPolygonalDomain(LineTes *tes, Rcpp::NumericMatrix coords) {
  tes->clear();
  Polygon domain;
  for (int i=0;i<coords.nrow();i++) {
    NT x(coords(i,0));
    NT y(coords(i,1));
    Point2 corner(x,y);
    domain.push_back(corner);
  }
  tes->insert_window(domain);
}
Rcpp::NumericMatrix getDomain(LineTes *tes) {
  Polygon w = tes->get_window();
  Rcpp::NumericMatrix rw(w.size(),2);
  for (int i=0;i!=w.size();i++) {
    rw(i,0) = CGAL::to_double(w[i][0]);
    rw(i,1) = CGAL::to_double(w[i][1]);
  }
  return rw;
}
Rcpp::NumericMatrix getSegmentCoords(LineTes *tes) {
  int n = (int) tes->number_of_segments();
  Rcpp::NumericMatrix coords(n,4); 
  int i = 0;
  for (LineTes::Seg_list_iterator si=tes->segments_begin();
       si!=tes->segments_end();si++) {
    coords(i,0) = CGAL::to_double((*si)->halfedges_start()->source()->point().x().exact());
    coords(i,1) = CGAL::to_double((*si)->halfedges_start()->source()->point().y().exact());
    coords(i,2) = CGAL::to_double((*si)->halfedges_end()->target()->point().x().exact());
    coords(i,3) = CGAL::to_double((*si)->halfedges_end()->target()->point().y().exact());
    i++;
  }
  Rcpp::List cnames = Rcpp::List::create(R_NilValue,
                                          Rcpp::CharacterVector::create("x0",
                                                                        "y0",
                                                                        "x1",
                                                                        "y1"));
  coords.attr("dimnames") = cnames;
  return coords;
}
void insert_segment(LineTes *tes,double x0,double y0,double x1,double y1) {
  Segment seg(Point2(x0,y0),Point2(x1,y1));
  tes->insert_segment(seg);
}
void remove_ivertices(LineTes *tes, int i=0, bool v=false) {
  tes->remove_ivertices(i,v);
}
void remove_lvertices(LineTes *tes, int i=0, bool v=false) {
  tes->remove_lvertices(i,v);
}
bool is_a_T_tessellation(LineTes *tes, bool verbose) {
  return tes->is_a_T_tessellation(verbose,Rcpp::Rcout);
}
void read(LineTes *tes,std::string filename) {
  std::ifstream input_file;
  input_file.open(filename.c_str(),std::ios::in);
  if (!input_file.is_open()) {
    Rcpp::Rcout << "Could not open file " << filename << "." << std::endl;
    Rcpp::Rcout << "Reading a line tessellation from file aborted." ;
    Rcpp::Rcout << std::endl;
    return;
  }
  tes->read(input_file);
  return;
}
void write(LineTes *tes,std::string filename) {
  std::ofstream output_file;
  output_file.open(filename.c_str(),std::ios::out);
  if (!output_file.is_open()) {
    Rcpp::Rcout << "Could not open file " << filename << "." << std::endl;
    Rcpp::Rcout << "Writing a line tessellation from file aborted." ;
    Rcpp::Rcout << std::endl;
    return;
  }
  tes->write(output_file);
  return;
}

int getSegmentNumber(TTessel *tes) {
  return (int) tes->number_of_segments();
}
int getInternalSegmentNumber(TTessel *tes) {
  int n;
  n = (int) tes->number_of_segments();
  n -= (int) tes->number_of_window_edges();
  return n;
}
double getInternalLength(TTessel *tes) {
  double l;
  l = tes->get_total_internal_length();
  l -= tes->get_window_perimeter();
  return l;
}
double getSumAngles(TTessel *tes) {
  double a = sum_of_min_angles(tes);
  return a;
  }  
void makeASplit(TTessel *tes) {    
  tes->update(tes->propose_split());
}
Rcpp::NumericVector getAcuteAngles(TTessel *tes, bool onlyInternal) {
  Rcpp::NumericVector res(2*tes->number_of_segments(),0.0);
  Rcpp::NumericVector::iterator res_it = res.begin();
  for (TTessel::Face_iterator f = tes->faces_begin();f!=tes->faces_end();f++) {
    if(f->is_unbounded()) continue;
    Polygon poly = face2poly(f);
    for(int i = 0;i != poly.size();i++) { 
      int j = (i-1)%((int)poly.size());
      if (j <0 ) j = j+poly.size();
      Point2 p = poly[i];
      double tokeep ;
      if(onlyInternal) {
	tokeep = is_point_inside_window(p,tes);
      } else {
	tokeep = 1.0;
      }
      if (tokeep>.5) {
	Vector v1 = Vector(poly[i],poly[j]);
	Vector v2 = Vector(poly[i],poly[(i+1)%((int)poly.size())]);
	double angle = angle_between_vectors(v1,v2);
	if (angle<CGAL_PI/2) {
	  *res_it += angle;
	  res_it++;
	}
      }
    }
  }
  Rcpp::NumericVector res_red(res.begin(),res_it);
  return res_red;
}
// To be changed when the domain will be a polygon
Rcpp::NumericVector getCellAreas(TTessel *tes) {
  Rcpp::NumericVector res(tes->number_of_segments()-4+1);
  Rcpp::NumericVector::iterator res_it = res.begin();
  for (TTessel::Face_iterator f = tes->faces_begin();f!=tes->faces_end();f++) {
    if (!f->is_unbounded()) { 
      Polygon poly = face2poly(f);
      *res_it = to_double(poly.area());
      res_it++;
    }
  }
  return res;
}
double getSumSquaredFaceAreas(TTessel *tes) {
  double a2;
  a2 = sum_of_faces_squared_areas(tes);
  return a2;
}

void ttessel_print_rcali(TTessel *tes) {
  tes->printRCALI(Rcpp::Rcout);
}

Rcpp::NumericMatrix getAdjacencyLengths(TTessel *tes) {
  int n = getInternalSegmentNumber(tes)+1;
  TTessel::Unbounded_face_iterator fstart = tes->unbounded_faces_begin();
  TTessel::Unbounded_face_iterator fend = tes->unbounded_faces_end();
  Rcpp::NumericMatrix res(n,n);
  for (TTessel::Edge_iterator e=tes->edges_begin();e!=tes->edges_end();e++) {
    double len = e->get_length();
    TTessel::Face_handle f1 = e->face();
    TTessel::Face_handle f2 = e->twin()->face();
    if (!f1->is_unbounded() && !f2->is_unbounded()) {
      // find face indices. Rather tedious!
      int i1 = n;
      int i2 = n;
      int i = 0;
      for (TTessel::Face_iterator fi=fstart;fi!=fend;fi++) {
	if(fi->is_unbounded())
	  continue;
	TTessel::Face_handle f = fi;
	if (f==f1) 
	  i1 = i;
	if (f==f2)
	  i2 = i;
	if (i1<n && i2<n)
	  break;
	i++;
      }
      if (i1>=n || i2>=n) {
	Rcpp::Rcout << "error: face index out of range" << std::endl;
      } else {
	res(i1,i2) = len;
	res(i2,i1) = len;
      }
    }
  }
  return res;
}

Energy ModelCRTT(double tau) {
  Energy mod;
  mod.add_features_segs(minus_is_segment_internal);
  mod.add_theta_segs(log(tau));
  return mod;
}
Energy ModelACS(double tau) {
  Energy mod;
  mod.add_features_segs(minus_is_segment_internal);
  mod.add_theta_segs(log(tau));
  mod.add_features_edges(edge_length);
  mod.add_theta_edges(tau/CGAL_PI);
  mod.add_features_vertices(is_point_inside_window);
  mod.add_theta_vertices(log(2.0));
  return mod;
}
Energy ModelSquaredAreas(double tau, double alpha) {
  Energy mod;
  mod.add_features_segs(minus_is_segment_internal);
  mod.add_theta_segs(log(tau));
  mod.add_features_faces(face_area_2);
  mod.add_theta_faces(alpha);
  return mod;
}
Energy ModelAcuteAngles(double tau, double beta) {
  Energy mod;
  mod.add_features_segs(minus_is_segment_internal);
  mod.add_theta_segs(log(tau));
  mod.add_features_faces(face_sum_of_angles);
  mod.add_theta_faces(beta);
  return mod;
}
Energy ModelAreasAngles(double tau, double alpha, double beta) {
  Energy mod;
  mod.add_features_segs(minus_is_segment_internal);
  mod.add_theta_segs(log(tau));
  mod.add_features_faces(face_area_2);
  mod.add_theta_faces(alpha);
  mod.add_features_faces(face_sum_of_angles);
  mod.add_theta_faces(beta);
  return mod;
}

namespace Rcpp {
  template<> SEXP wrap(const CatVector& par) {
    std::map<std::string,std::vector<double> > list;
    std::vector<double> vertices;
    std::vector<double> edges;
    std::vector<double> faces;
    std::vector<double> segs;
    for (unsigned int i=0; i!=par.vertices.size();i++)
      vertices.push_back(par.vertices[i]);
    list["vertices"] = vertices;
    for (unsigned int i=0; i!=par.edges.size();i++)
      vertices.push_back(par.edges[i]);
    list["edges"] = edges;
    for (unsigned int i=0; i!=par.faces.size();i++)
      faces.push_back(par.faces[i]);
    list["faces"] = faces;
    for (unsigned int i=0; i!=par.segs.size();i++)
      segs.push_back(par.segs[i]);
    list["segs"] = segs;
    return Rcpp::wrap(list);
  }
}

namespace Rcpp {
  template <> CatVector as(SEXP R_list) {
    typedef std::vector<double> nv;
    CatVector par;
    SEXP names = Rf_getAttrib(R_list,R_NamesSymbol);
    for(int i=0;i<Rf_length(R_list);i++) {
      if (!strcmp(CHAR(STRING_ELT(names, i)), "vertices")) {
	Rcpp::NumericVector x(VECTOR_ELT(R_list,i));
	par.vertices.assign(x.begin(),x.end());
      } else if (!strcmp(CHAR(STRING_ELT(names, i)), "edges")) {
	Rcpp::NumericVector x(VECTOR_ELT(R_list,i));
	par.edges.assign(x.begin(),x.end());
      } else if (!strcmp(CHAR(STRING_ELT(names, i)), "faces")) {
	Rcpp::NumericVector x(VECTOR_ELT(R_list,i));
	par.faces.assign(x.begin(),x.end());
      } else if (!strcmp(CHAR(STRING_ELT(names, i)), "segs")) {
	Rcpp::NumericVector x(VECTOR_ELT(R_list,i));
	par.segs.assign(x.begin(),x.end());
      }
    }
    return par;
  }
}

namespace Rcpp {
  template<> SEXP wrap(const Polygons& polygons) {
    std::vector<NumericMatrix> list;
    for (unsigned int i=0; i!=polygons.size(); i++) {
      Polygon p = polygons[i];
      NumericMatrix coord(p.size(),2);
      for (unsigned int j=0; j!=p.size(); j++) {
	coord(j,0) = CGAL::to_double(p[j][0]);
	coord(j,1) = CGAL::to_double(p[j][1]);
      }
      list.push_back(coord);
    }
    return Rcpp::wrap(list);
  }
}
// Rcpp::NumericVector R_get_theta(Energy *e) {
//   Parameters par = e->get_theta();
//   Rcpp::NumericVector theta = Parameters2NumericVector(par);
//   return theta;
// }

// void R_set_theta(Energy *e, Rcpp::NumericVector theta) {
//   Parameters par = NumericVector2Parameters(theta,e);
//   e->set_theta(par);
// }

void setLiTeSeed(int seed) {
  if(rnd) 
    delete rnd;
  rnd = new CGAL::Random(seed);
}

Rcpp::NumericMatrix SMFStep(SMFChain *smf, int n) {
  ModifCounts track;
  Rcpp::NumericMatrix res(2,3);
  track = smf->step(n);
  res(0,0) = track.proposed_S;
  res(1,0) = track.accepted_S;
  res(0,1) = track.proposed_M;
  res(1,1) = track.accepted_M;
  res(0,2) = track.proposed_F;
  res(1,2) = track.accepted_F;
  Rcpp::List dnames = Rcpp::List::create(Rcpp::CharacterVector::create("proposed",
								       "accepted"),
					 Rcpp::CharacterVector::create("split",
								       "merge",
								       "flip"));
  res.attr("dimnames") = dnames;
  return res;
}

Rcpp::NumericMatrix PLGetHessian(PseudoLikDiscrete* pld,CatVector theta) {
  CatMatrix h0 = pld->GetHessian(theta);
  FMatrix h1 = asFMatrix(h0);
  Rcpp::NumericMatrix h2(h1.row_dimension(),h1.column_dimension());
  for (int i=0;i!=h1.row_dimension();i++) {
    for (int j=0;j!=h1.column_dimension();j++) {
      h2(i,j) = h1(i,j);
    }
  }
  return h2;
}
// Rcpp::NumericVector PseudoLikDiscreteGetGradient(PseudoLikDiscrete *pld, 
// 				  Rcpp::NumericVector RTheta) {
//   Parameters theta = NumericVector2Parameters(RTheta,
// 					      pld->GetEnergy());
//   Parameters g = pld->GetGradient(theta);
//   Rcpp::NumericVector Rg = Parameters2NumericVector(g);
//   return Rg;
// }

RCPP_MODULE(lite){
  using namespace Rcpp;
  class_<LineTes>("LineTes")
    .constructor("generate an empty LineTes object")
    /* Something weird here. It is not possible to define the setDomain method
       for LineTes and let TTessel inherit it. Furthermore, it is not possible
       to define two methods for LineTes and TTessel with the same name but
       "different" implementations.*/
    .method("is_valid",&LineTes::is_valid,
	    "test validity of the line tessellation")
    .method("setDomain",&setDomain,
	    "define the rectangular domain to be tessellated")
    .method("setDomain",&setPolygonalDomain,
	    "define the polygonal domain to be tessellated")
    .method("getDomain",&getDomain,"get the domain to be tessellated")
    .method("getSegmentCoords",&getSegmentCoords,
	    "return the segment coordinates")
    .method("insert_segment",&insert_segment,"insert a line segment in a tessellation")
    .method("remove_ivertices",&remove_ivertices,"remove I-vertices from a tessellation")
    .method("remove_lvertices",&remove_lvertices,"remove internal L-vertices from a tessellation")
    .method("is_a_T_tessellation",&is_a_T_tessellation,"test whether the tessellation is of T-type")
    .method("read",&read,"read a line tessellation from file")
    .method("write",&write,"write a line tessellation to a file")
    ;
  class_<TTessel>("TTessel")
    .derives<LineTes>("LineTes")
    .constructor("generate an empty TTessel object")
    .constructor<LineTes&>("generate a TTessel object from a LineTes object")
    .method("is_valid",&TTessel::is_valid,"test validity of a T-tessellation")
    .method("clear",&TTessel::clear,"empty the tessellation")
    .method("getSegmentNumber",&getSegmentNumber,
	    "return the total number of segments")
    .method("getInternalSegmentNumber",&getInternalSegmentNumber,
	    "return the number of internal segments")
    .method("number_of_blocking_segments",&TTessel::number_of_blocking_segments,
	    "return the number of internal blocking segments")
    .method("number_of_non_blocking_segments",
	    &TTessel::number_of_non_blocking_segments,
	    "return the number of internal non-blocking segments")
    .method("totalPerimeter",&TTessel::get_total_internal_length,
	    "return the sum of all cell perimeters")
    .method("getInternalLength",&getInternalLength,
	    "return the total length of internal segments")
    .method("all_faces",&TTessel::all_faces,"return all tessellation faces")
    .method("getSumAngles",&getSumAngles,
	    "return the sum of vertex acute angles")
    .method("makeASplit",&makeASplit,
	    "split the T-tessellation at random")
    .method("getAcuteAngles",&getAcuteAngles,"return vertex acute angles")
    .method("getCellAreas",&getCellAreas,"return cell areas")
    .method("getSumSquaredFaceAreas",&getSumSquaredFaceAreas,
	    "return the sum of squared face areas")
    .method("printRCALI",&ttessel_print_rcali,
	    "Print the tessellation in a format supported by califlopp")
    .method("getAdjacencyLengths",&getAdjacencyLengths,
	    "return length wieighted adjacency matrix for cells")
    ;
  /* Pour tester l'estimation, il faut que je puisse spécifier des sous-modèles
   du modèle le plus complet (celui qui est implémenté actuellement). */
  class_<Energy>("Energy")
    .constructor()
    .method("set_ttessel",&Energy::set_ttessel,"set the TTessel data member")
    .method("get_ttessel",&Energy::get_ttessel,"get the TTessel data member")
    .method("getValue",&Energy::get_value,"return the energy value")
    .method("set_theta",&Energy::set_theta,"set the parameter vector")
    .property("theta",&Energy::get_theta,&Energy::set_theta,
	      "access to the parameter vector")
  ;
  function("ModelCRTT",&ModelCRTT,List::create(_["tau"]=1.0),
	   "create CRTT model as an Energy object");
  function("ModelACS",&ModelACS,List::create(_["tau"]),
	   "create ACS model as an Energy object");
  function("ModelSquaredAreas",&ModelSquaredAreas,
	   List::create(_["tau"],_["alpha"]),
	   "create model with penalty on squared areas as an Energy object");
  function("ModelAcuteAngles",&ModelAcuteAngles,
	   List::create(_["tau"],_["beta"]),
	   "create model with penalty on acute angles as an Energy object");
  function("ModelAreasAngles",&ModelAreasAngles,
	   List::create(_["tau"],_["alpha"],_["beta"]),
	   "create model with penalty on squared areas and acute angles");
  class_<SMFChain>("SMFChain")
  .constructor()
  .constructor<Energy*,double,double>()
    .method("set_energy",&SMFChain::set_energy,"set the model to be simulated")
    .method("get_energy",&SMFChain::get_energy,"get the model to be simulated")
    .method("setSMFProb",&SMFChain::set_smf_prob,
	    "set split and merge probabilities")
    .method("step",&SMFStep,"proceed to simulation iteration")
  ;
  class_<PseudoLikDiscrete>("PseudoLikDiscrete")
    .constructor("generate an empty PseudoLikDiscrete object")
    .constructor<Energy*>("generate a PseudoLikDiscrete object associated to an Energy object to be fitted")
    .method("GetEnergy",&PseudoLikDiscrete::GetEnergy,
	    "Access to current fitted Model")
    .method("SetEnergy",&PseudoLikDiscrete::SetEnergy,
	    "Modify the current fitted Model")
    .method("AddSplits",&PseudoLikDiscrete::AddSplits,
	    "Increment the sample of splits involved in the integral component of the pseudo-likelihood")
    .method("ClearSplits",&PseudoLikDiscrete::ClearSplits,
	    "Clear the sample of splits")
    .method("GetValue",&PseudoLikDiscrete::GetValue,
	    "Access to the current value of the discretized pseudo-likelihood")
    .method("GetGradient",&PseudoLikDiscrete::GetGradient,
    "Access to the gradient of the discretized pseudo-likelihood")
    .method("GetHessian",&PLGetHessian,
	    "Access to the Hessian of the discretized pseudo-likelihood")
    ;
  class_<PLInferenceNOIS>("PLInferenceNOIS")
    .derives<PseudoLikDiscrete>("PseudoLikDiscrete")
    .constructor<Energy*,bool,double>("PLInferenceNOIS constructor")
    .method("Step",&PLInferenceNOIS::Step,
	    "Iterate the estimation algorithm")
    .method("GetEstimate",&PLInferenceNOIS::GetEstimate,"Access to the"
	    " current parameter estimate")
    .method("Run",&PLInferenceNOIS::Run,
	    "Iterate until convergence")
    .method("GetEstimates",&PLInferenceNOIS::GetEstimates,
	    "Access to all intermediate estimates")
    .method("GetValues",&PLInferenceNOIS::GetValues,
	    "Access to all intermediate log-pseudolikelihood values")
    ;
  function("setLiTeSeed",&setLiTeSeed,"Initialize the random generator");
}
