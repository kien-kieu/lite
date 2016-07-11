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

/** \brief Compute the length of the shortest edge */
double LineTes::min_edge_length(){
  double mel=edges_begin()->get_length();
  double r;
  for (LineTes::Edge_iterator e=edges_begin();e!=edges_end();e++){
    r=e->get_length();
    if (r<mel){
      mel=r;
    }
  } 
  return mel; 
}

/** \brief Test whether a vertex is of X-type
 * \param v : handle to the vertex to be tested.
 * \param tes : pointer to the line tessellation where the vertex
 * comes from.
 * \return true if the vertex is an X-vertex and if it does not lie
 * on the boundary of the domain.
 *
 * An X-vertex is a vertex of degree 4 with two pairs of aligned incident
 * edges. 
 */
bool is_an_X_vertex(LineTes::Vertex_handle v) {
  if (v->degree()!=4) return false;
  LineTes::Halfedge_around_vertex_circulator e = v->incident_halfedges(),
    e0, e1, e2, e3;
  e0 = e; e1 = ++e; e2 = ++e; e3 = ++e;
  if (e0->get_next_hf()!=e2->twin()) return false;
  if (e1->get_next_hf()!=e3->twin()) return false;
  return true;
}
/** \brief Test whether a vertex is an irreducible X-vertex
 *
 * An X-vertex is said to be irreducible if it is still an X-vertex after
 * a removal of all I-vertices. I.e. if at least one of its incident edge 
 * ends with a I-vertex.
 */
bool is_an_irreducible_X_vertex(LineTes::Vertex_handle v) {
  if (!is_an_X_vertex(v)) return false;
  LineTes::Halfedge_around_vertex_circulator e_x = v->incident_halfedges(),
    done = e_x;
  do {
    if (e_x->source()->degree()==1) return false;
    e_x++;
  } while (e_x!=done);
  return true;
}
/** \brief Test whether a vertex lies along the domain boundary
 *
 * Return
 * - 0 if the vertex lies in the interior of the tessellated domain.
 * - 1 if the vertex lies on the boundary of the domain.
 */
int LineTes::is_on_boundary(Vertex_handle v) {
  Halfedge_around_vertex_circulator hc =  v->incident_halfedges(), 
    done = hc;
  do {
    int res = is_on_boundary(hc);
    if (res>0) return 1;
    hc++;
  } while (hc!=done);
  return 0;
}
/* For use by function remove_xvertices. Determine whether a vertex
   should be removed. A vertex is to be removed if it is an
   irreducible X-vertex lying inside the domain. */
bool rxv_select(LineTes::Vertex_handle v, LineTes* tes) {
  return is_an_irreducible_X_vertex(v) && tes->is_on_boundary(v)==0;
}
/** \brief Remove all internal and irreducible X-vertices from a line 
 * tessellation
 * \param step : lenghth of vertex move.
 * 
 * See function is_an_irreducible_X_vertex() for a definition of an irreducible 
 * X-vertex. The removal of an X-vertex generates two T-vertices. The procedure
 * for a removal 
 * is as follows:
 * - One of the segment passing through the vertex is broken into two halves
 *   at the vertex location.
 * - The end of a segment half lying at the vertex location is moved along the
 *   other incident segment a little bit.
 * The procedure is repeated until there is no internal irreducible X-vertex 
 * left in the 
 * tessellation. Note that this procedure can generate new I-vertices.
 * \sa remove_ivertices, remove_lvertices.
 */
void  LineTes::remove_xvertices(double step){
  double eps;
  std::vector<TTessel::Vertex_handle> XVertices; 
  TTessel::Vertex_handle v, vs ; 
  LineTes::Halfedge_around_vertex_circulator  curr; 
  LineTes::Halfedge_handle he, he_shadow;
  LineTes::Seg_handle s;
  
 // initializing a vector of X-vertices
   
  for (TTessel::Vertex_iterator v=vertices_begin();v!=vertices_end();v++){
    if (rxv_select(v,this)){
      XVertices.push_back(v);
    }  
  } 
  
  while (XVertices.size()>0){
    v=XVertices[0];
    // a halfedge belonging to a segment to be disturbed
    curr=v->incident_halfedges(); 
    // a segment to be disturbed
    s=curr->segment();
    // if the direction of "curr" differs from the direction of s we
    // consider another halfedge adjacent to v, belonging to s and
    // sharing its direction
    if (!curr->get_dir()){
      curr=curr->get_next_hf()->twin();
    }
    // we move (clockwise) to the next halfedge adjacent to v:
    // this halfedge contains the new endpoint of the perturbed segment
    curr++;
    // perturbation length;
    eps = 0.9*min_edge_length();
    if (step < eps){
      eps = step;
    }
    // moving an endpoint of s
    Vector move(curr->target()->point(),curr->source()->point());
    move = eps*move/curr->get_length();
    Segment new_s(s->pointSource(),curr->target()->point() + move);
    
    LineTes::Seg_handle new_s_tes = insert_segment(new_s);
    
    // removing the edges of "s" from "tes"
    
    he= s->halfedges_start();
    while (he->source()!=v){ 
      he_shadow=he;  
      vs = he->target();
      // if vs is an X-vertex it will be removed from "XVertices" list
      if (rxv_select(vs,this)){
	for (int j=0;j<XVertices.size();j++){ // à simplifier
	  if (XVertices[j]==vs){
	    XVertices.erase(XVertices.begin()+j);	  
	  }
	}
      }
      he=he->get_next_hf();
      suppress_edge(he_shadow);
    }
    
    /* vertices along the displaced half segment are tested for
       updating the list of vertices to be removed */

    int n=new_s_tes->number_of_edges();
    
    curr = new_s_tes->halfedges_start();
    vs=curr->target();
    while (vs!=new_s_tes->halfedges_end()->target()){
      if (rxv_select(vs,this)){
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
