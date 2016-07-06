/* Line Tessellation (LiTe) library
   |||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
   |||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

/* Modifs de Kasia pour retirer les sommets en X
getVerticesCoords
clipLineByPolygonR
remove_xvertices
*/

#include "ttessel.h"

// To be moved to rlite_module.cpp
/* 
Rcpp::NumericMatrix getVerticesCoords(LineTes *tes) {  
  int n = (int) tes->number_of_vertices();
  Rcpp::NumericMatrix coords(n,3);
  int i=0;

  for (TTessel::Vertex_iterator v=tes->vertices_begin();v!=tes->vertices_end();v++){
    coords(i,0) = CGAL::to_double(v->point().x().exact());
    coords(i,1) = CGAL::to_double(v->point().y().exact());
    coords(i,2) = CGAL::to_double(v->degree());
    i++;
  }  

   Rcpp::List cnames = Rcpp::List::create(R_NilValue,
                                          Rcpp::CharacterVector::create("x0",
                                                                        "y1",
                                                                         "degree"));

  
  coords.attr("dimnames") = cnames;

  return coords;
}
*/

double min_edge_length(LineTes* tes){
  double mel=tes->edges_begin()->get_length();
  double r;
  for (LineTes::Edge_iterator e=tes->edges_begin();e!=tes->edges_end();e++){
    r=e->get_length();
    if (r<mel){
      mel=r;
    }
  } 
  return mel; 
}
bool is_an_X_vertex(LineTes::Vertex_handle v, LineTes *tes) {
  
  bool x=false;
  
  if (v->degree()==4 && is_point_inside_window(v->point(),tes)>0){
    x=true;
    LineTes::Halfedge_around_vertex_circulator e_x=v->incident_halfedges(); 
    for (int l=0;l<4;l++){
      x=x & e_x->source()->degree()!=1;
      e_x++;
    }
    
  }
    
  return x;
}

void  LineTes::remove_xvertices(double step){
  // replaces each X-vertex of the line tessellation by
  // two T-vertices by moving slightly the
  // endpoint of one of 4 adjacent segments.

  // step : the length of movement, if too large, "step"
  // is replaced by the minimal edge length

  double eps;
  std::vector<TTessel::Vertex_handle> XVertices; 
  TTessel::Vertex_handle v, vs ; 
  LineTes::Halfedge_around_vertex_circulator  curr; 
  LineTes::Halfedge_handle he, he_shadow;
  LineTes::Seg_handle s;
 
 // initializing a vector of X-vertices
   
   for (TTessel::Vertex_iterator v=vertices_begin();v!=vertices_end();v++){
        if (is_an_X_vertex(v,this)){
		 XVertices.push_back(v);
         }  
   } 
    
   
 while (XVertices.size()>0){

  v=XVertices[0];
  // a halfedge belonging to a segment to be disturbed
  curr=v->incident_halfedges(); 
  
  // a segment to be disturbed
  s=curr->segment();

  // if the direction of "curr" differs from the direction of s
  // we consider another halfedge adjacent to v, belonging to s and sharing its direction
  
  if (!curr->get_dir()){
        curr=curr->get_next_hf()->twin();
  }
  
  // we move (clockwise) to the next halfedge adjacent to v:
  // this halfedge contains the new endpoint of the perturbed segment

  curr++;

  // perturbation length;

  eps = 0.9*min_edge_length(this);
  if (step < eps){
    eps = step;
  }
  

   // moving an endpoint of s
   
  Segment new_s(s->pointSource(),
                curr->target()->point()+
                (eps/(curr->get_length()))*Vector(curr->target()->point(),curr->source()->point()));
  
  LineTes::Seg_handle new_s_tes=insert_segment(new_s);
  
  // removing the edges of "s" from "tes"
  
  he= s->halfedges_start();
  while (he->source()!=v){ 
         he_shadow=he;  
         vs = he->target();
         // if vs is an X-vertex it will be removed from "XVertices" list
         if (is_an_X_vertex(vs,this)){
	   for (int j=0;j<XVertices.size();j++){ // à simplifier
                  if (XVertices[j]==vs){
                      XVertices.erase(XVertices.begin()+j);	  
                      }
                 }
	   }
         he=he->get_next_hf();
         suppress_edge(he_shadow);
  }

  
  // any X-vertex preceeding "v" on "s" remains X-vertex after perturbation :
  // they should be added to "XVertices" list
         
         int n=new_s_tes->number_of_edges();
 
         
         curr = new_s_tes->halfedges_start();
         vs=curr->target();
          while (vs!=new_s_tes->halfedges_end()->target()){
	    if (is_an_X_vertex(vs,this)){
	        XVertices.push_back(vs);
	    }
           curr = curr->get_next_hf();
           vs=curr->target();     
	  }
	 
   
   
 }
 }






// To be moved to rlite_module.cpp
/*
Rcpp::NumericMatrix clipLineByPolygonR(double a, double b, LineTes* t){
  NT aNT=NT(a);
  NT bNT=NT(b);
  Polygon W=t->get_window();
  std::vector<Point2> getXY=clipLineByPolygon(aNT,bNT,W);
  Rcpp::NumericMatrix intersectionXY(2,4);
  for (int i=0;i<2;i++){
    for (int j=0;j<4;j++){
      intersectionXY(i,j)=CGAL::to_double(getXY[i][j]);
    }
  } 
  return intersectionXY;
}
*/

// Modifications to be moved to rlite_module.cpp
/*
RCPP_MODULE(lite){
  using namespace Rcpp;
  class_<LineTes>("LineTes")
    .method("getVerticesCoords",&getVerticesCoords,
	    "ka : return the  coordinates of all vertices together with their degree")
    .method("remove_xvertices",&remove_xvertices,
	    "ka : transform each X-vertex into 2 T-vertices")
    ;
  function("clipLineByPolygon",&clipLineByPolygonR,"toto");
}
*/
