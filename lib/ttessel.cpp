/* Line Tessellation (LiTe) library
   |||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
   |||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

/******************************************************************************/
/*             LINE BASED TESSELLATION CLASS  - METHODS                      */
/******************************************************************************/

#include "ttessel.h"

/******************************************************************************/
/*                              GLOBAL VARIABLES                              */
/******************************************************************************/

LineTes::Halfedge_handle NULL_HALFEDGE_HANDLE;
LineTes::Seg_handle NULL_SEG_HANDLE;
LineTes::Face_handle NULL_FACE_HANDLE;
CGAL::Random *rnd = 0;
//rnd = new CGAL::Random(seed);


/******************************************************************************/
/*                   METHODS FOR THE CLASS LineTes                            */
/******************************************************************************/

/** \brief Empty LineTes constructor
 *
 * Equivalent to the empty CGAL::Arrangement_2 constructor
 */
LineTes::LineTes() : Arrangement() {}
/** \brief Define the rectangular domain delimiting the tessellation
 *
 * The LineTes class is designed for representing tessellation of
 * bounded rectangular domains. The insert_window method should be
 * used as a first step where one defines the domain to be
 * tessellated.
 */
void LineTes::insert_window(Rectangle r) {
  Polygon p;
  for (int i=0;i<4;i++) {
    p.push_back(r[i]);
  }
  insert_window(p);
}
/** \brief Define the polygonal domain delimiting the tessellation
 *
 * Do not use this method at the moment as non-rectangular domains
 * are not fully supported by the LineTes class.
 */
void LineTes::insert_window(Polygon& p) {
  Halfedge_handle e;
  Vertex_handle v,v0;

  // Reset private member window
  window.clear();
  /* Start by removing all unnecessary vertices */
  Polygon::Vertex_circulator pv2 = p.vertices_circulator();
  pv2--;
  Polygon::Vertex_circulator pv1 = pv2++, pv0 = pv1++, pv_end = pv1;
  pv2++;
  do {
    if(!CGAL::collinear(*pv0,*pv1,*pv2)) {
      window.push_back(*pv1);
    }
    pv0++; pv1++; pv2++;
  } while (pv1!=pv_end);

  /* Test whether the polygon is convex */
  // if (!CGAL::is_convex_2(window.vertices_begin(),window.vertices_end())) {
  //   std::cerr << "Window must be a convex polygon" << std::endl;
  //   std::cerr << "Proposed window:" << std::endl;
  //   std::cerr << window << std::endl;
  //   exit(EXIT_FAILURE);
  // }

  /* Now insert the window edges in the arrangement */

  v = insert_in_face_interior(window[0],unbounded_face());
  v0 = v;

  for (int i=0;i<window.size()-1;i++) {
    if (window[i]<=window[i+1]) {
      e = insert_from_left_vertex(Curve(window[i],window[i+1]),v);
    } else {
      e = insert_from_right_vertex(Curve(window[i+1],window[i]),v);
    }
    v = e->target();
    e->set_length();
    e->twin()->set_length(e->get_length());
    set_junction(e,NULL_HALFEDGE_HANDLE);
    set_junction(NULL_HALFEDGE_HANDLE,e);
    e->set_dir(true);
    e->twin()->set_dir(false);
    Seg_handle s = new Seg;
    s->set_halfedge_handle(e);
    e->set_segment(s); e->twin()->set_segment(s);
    all_segments.add(s);
  } 
  // Last insertion is particular because vertices pre-exist at both ends of the added curve
  e = insert_at_vertices(Curve(window[window.size()-1],window[0]),v,v0);
  v = e->target();
  e->set_length();
  e->twin()->set_length(e->get_length());
  set_junction(e,NULL_HALFEDGE_HANDLE);
  set_junction(NULL_HALFEDGE_HANDLE,e);
  e->set_dir(true);
  e->twin()->set_dir(false);
  Seg_handle s = new Seg;
  s->set_halfedge_handle(e);
  e->set_segment(s); e->twin()->set_segment(s);
  all_segments.add(s);
  if (!is_valid()) {
    std::cerr << "arrangement is not valid: " << std::endl;
  }
}
/** \brief Split a tessellation cell from an edge to another one
 * 
 * \param e1 : handle of the halfedge where the splitting segment starts.
 * \param p1 : point on the edge where the splitting segment starts.
 * \param e2 : handle of the halfedge where the splitting segment ends.
 * \param p2 : point on the edge where the splitting segment ends.
 * \return  handle of one of the halfedges along the splitting segment.
 * \pre e1 should not bound the external face. If so, return 
 * NULL_HALFEDGE_HANDLE.
 */
LineTes::Halfedge_handle LineTes::split_from_edge(Halfedge_handle e1, Point2 p1,
						  Halfedge_handle e2, Point2 p2)
{
  Point2                     seg_deb,seg_fin; 
  Halfedge_handle            hf_it_deb,new_e;
  CGAL::Object               inter;

  if (is_on_boundary(e1)==2) 
    return NULL_HALFEDGE_HANDLE;

  Halfedge_handle e1_new(split_edge(e1,Curve(e1->source()->point(),p1),
				      Curve(p1,e1->target()->point()))); 
  Halfedge_handle e2_new(split_edge(e2,Curve(e2->source()->point(),p2),
				    Curve(p2,e2->target()->point())));
  new_e = add_edge(Curve(p1,p2),e1_new->target(),e2_new->target());

  return new_e;
}
/** \brief Split a tessellation cell from a vertex to an edge
 * \param v : handle of the vertex where the splitting segment starts.
 * \param e : handle of the edge where the splitting segment ends.
 * \param p: point where the splitting segment ends.
 * \return  handle of one of the halfedges along the splitting segment.
 */
LineTes::Halfedge_handle LineTes::split_from_vertex(Vertex_handle v,
						    Halfedge_handle e,
						    Point2 p) {
  bool merge_required = false;
  Curve c;
  if (v->degree()==1) {
    Point2 pprev = v->incident_halfedges()->source()->point();
    Point2 pmid = v->point();
    if (CGAL::collinear(pprev,pmid,p) && 
	CGAL::angle(pprev,pmid,p)==CGAL::OBTUSE) {
      merge_required = true;
      c = Curve(pprev,p);
    }
  }
  Halfedge_handle e2 = split_edge(e,Curve(e->source()->point(),p),
		  Curve(p,e->target()->point()));
  Halfedge_handle new_e = add_edge(Curve(v->point(),p),v,e2->target());
  if (merge_required) {
    new_e = merge_edge(v->incident_halfedges(),new_e,c);
  }
  return new_e;
}
/** \brief Suppress an edge from the tessellation
 * \param e : handle of the halfedge to be suppressed. 
 * \return a handle to the new face created by the edge removal.
 * \pre e should not lie on the domain boundary and it should be a the
 * end of an existing segment. If the condition is not fulfilled,
 * the tessellation is not modified and NULL_HALFEDGE_HANDLE is
 * returned.
 *
 * The vertices of the suppressed edge are removed if they are of
 * degree 2 and if the connect two consecutive edges on a segment.
 */
LineTes::Face_handle LineTes::suppress_edge(Halfedge_handle e) {
  if (is_on_boundary(e)>0 || (e->get_prev_hf()!=NULL_HALFEDGE_HANDLE && 
			      e->get_next_hf()!=NULL_HALFEDGE_HANDLE))
    return NULL_FACE_HANDLE;

  Face_handle f;
  Vertex_handle v1 = e->source() ,v2 = e->target();
  
  // Remove edge
  f = remove_edge(e);

  if (f==NULL_FACE_HANDLE)
    return NULL_FACE_HANDLE;

  if (v1->is_isolated())
    remove_isolated_vertex(v1);
  if (v1->degree()==2 && 
      v1->incident_halfedges()->get_next_hf()!=NULL_HALFEDGE_HANDLE) {
    merge_edge(v1->incident_halfedges(),
	       v1->incident_halfedges()->get_next_hf(),
	       Curve(v1->incident_halfedges()->source()->point(),
		     v1->incident_halfedges()->get_next_hf()->target()->point()));
  }

  if (v2->is_isolated())
    remove_isolated_vertex(v2);
  if (v2->degree()==2 && 
      v2->incident_halfedges()->get_next_hf()!=NULL_HALFEDGE_HANDLE) {
    merge_edge(v2->incident_halfedges(),
	       v2->incident_halfedges()->get_next_hf(),
	       Curve(v2->incident_halfedges()->source()->point(),
		     v2->incident_halfedges()->get_next_hf()->target()->point()));
  }

  return f;
}
/** \brief Clear a tessellation
 *
 * All segments are removed including the domain sides.
 */
void LineTes::clear() {
  Arrangement::clear();
  all_segments.clear();
}
/** \brief Check whether the LineTes object is valid
 * \param verbose : if true, details are sent to std::clog. Default to false.
 * return true if valid, false otherwise
 *
 * The following checks are performed:
 * - Validity as an arrangement.
 * - Validity of all halfedges (method check_all_halfedges).
 * - Validity of all segments (method check_all_segments).
 */
bool LineTes::is_valid(bool verbose) {
  if (!Arrangement::is_valid()) {
    if (verbose) {
      std::clog << "LineTes::is_valid: not a valid arrangement" << std::endl;
    }
    return false;
  }
  if (!check_all_halfedges(verbose)) {
    return false;
  }
  if (!check_all_segments(verbose)) {
    return false;
  }
  return true;
}
/** \brief Test whether all segments of a tessellation are valid
 * \return true if all segments are valid, false otherwise.
 * 
 * All segments are tested using the check_segment method.
 */
bool LineTes::check_all_segments(bool verbose) {
  Size n = 0;
  for (Seg_list_iterator si=segments_begin();si!=segments_end();si++) {
    if (!check_segment(*si,verbose)) {
      n++;
    }
  }
  if (n>0 && verbose) {
    std::clog << "LineTes::check_all_segments: found " << n;
    std::clog << " invalid segments" << std::endl;
  }
  return (n==0);
}
/** \brief Check all halfedges of a line tessellation
 *
 * \param verbose : if true, when a invalid halfedge is found, print 
 *                  some data on the standard output stream for logging.
 *
 * All halfedges of the line tessellation are inspected. Check is based
 * on the LineTes::check_halfedge method.
*/
bool LineTes::check_all_halfedges(bool verbose) {
  bool res = true;
  for(Halfedge_iterator e=halfedges_begin();e!=halfedges_end();e++) {
    if (!check_halfedge(e,verbose)) {
      if (verbose) {
	std::clog << "halfedge " << e->source()->point() << " ";
	std::clog << e->target()->point() << " is invalid" << std::endl;
      }
      res = false;
    }
  }
  return res;
}
/** \brief Test whether a segment of a tessellation is valid
 * \param s : a handle of the segment to be tested.
 * \param verbose : if true and the segment is non valid, some
 *                  information is sent to the standard output
 *                  for logging.
 * \return true if the segment is valid, false otherwise.
 *
 * Just check that the halfedge identifying the segment is indeed
 * a tessellation halfedge.
 */
bool LineTes::check_segment(Seg_handle s, bool verbose) {
  Halfedge_handle e = s->get_halfedge_handle();
  if (!halfedge_exists(e)) {
    std::clog << "Halfedge does not exist!" << std::endl;
    return false;
  }
  return true;
}
/** \brief Test whether a halfedge is valid
 * \param e : a handle of the halfedge to be checked.
 * \param verbose : if true and halfedge is found to be non valid,
 *                  some information is sent to the standard output
 *                  for logging.
 * \return true if the halfedge is valid, false otherwise.
 *
 * Several tests are performed: existence of the halfedge, validity of its
 * neighbours along its segment, consistency of the segment and direction 
 * attributes.
 */
bool LineTes::check_halfedge(Halfedge_handle e, bool verbose) {
  if (!halfedge_exists(e)) {
    if (verbose) {
      std::clog << "LineTes::check_halfedge: halfedge not found in halfedge list";
      std::clog << std::endl;
    }
    return false;
  }
  if (!check_halfedge_neighbours(e, verbose)) {
    return false;
  }
  if (!check_halfedge_segment(e, verbose)) {
    return false;
  }
  return true;
}
/** \brief Test whether a halfedge handle is indeed a halfedge handle of a 
 *  tessellation
 */
bool LineTes::halfedge_exists(Halfedge_handle e) {
  Halfedge_iterator ou_e = halfedges_begin();
  while (ou_e!=e && ou_e!=halfedges_end())
    ou_e++;
  if (ou_e==halfedges_end())
    return false;
  else
    return true;
}
/** \brief Check the neighbours of a halfedge along its segment
 * \param e : a handle of the halfedge to be checked.
 * \param verbose : if true, details are sent to std::clog. Default to false.
 * \return true if neighbours are valid, false otherwise.
 */
bool LineTes::check_halfedge_neighbours(Halfedge_handle e, bool verbose) {
  bool res = check_halfedge_prev_hf(e, verbose);
  res = res && check_halfedge_next_hf(e, verbose);
  return res;
}
/** \brief Check the previous neighbour of a halfedge along its segment
 * \param e : a handle of the halfedge to be checked.
 * \param verbose : if true, details are sent to std::clog. Default to false.
 * \return true if the previous neighbour is valid, false otherwise.
 */
bool LineTes::check_halfedge_prev_hf(Halfedge_handle e, bool verbose) {
  if (e->get_prev_hf()==NULL_HALFEDGE_HANDLE) {
    Halfedge_around_vertex_circulator ie0 = e->source()->incident_halfedges();
    Halfedge_around_vertex_circulator ie = ie0;
    do {
      if(ie!=e->twin()) {
	// check ie and e are not aligned
	if (are_aligned(ie->source()->point(),e->source()->point(),
			e->target()->point())) {
	  if (verbose) {
	    std::clog << "LineTes::check_halfedge_prev_hf method: prev_hf attribute";
	    std::clog << " is null, but found a previous neighbour along the ";
	    std::clog << "segment" << std::endl;
	  }
	  return false;
	}
      }
      ie++;
    } while(ie!=ie0);
  } else {
    // check connection
    if (e->get_prev_hf()->target()!=e->source()) {
      if (verbose) {
	std::clog << "LineTes::check_halfedge_prev_hf method: previous halfedge";
	std::clog << " is not connected to current halfedge" << std::endl;
      }
      return false;
    }
    // check e->get_prev_hf and e are aligned
    if(!are_aligned(e->get_prev_hf()->source()->point(),e->source()->point(),
		    e->target()->point())) {
      if (verbose) {
	std::clog << "LineTes::check_halfedge_prev_hf method: previous halfedge ";
	std::clog << "is not aligned with the current halfedged" << std::endl;
      }
      return false;
    }
  }
  return true;
}
/** \brief Check the next neighbour of a halfedge along its segment
 * \param e : a handle of the halfedge to be checked.
 * \param verbose : if true, details are sent to std::clog. Default to false.
 * \return true if the next neighbour is valid, false otherwise.
 */
bool LineTes::check_halfedge_next_hf(Halfedge_handle e, bool verbose) {
  return check_halfedge_prev_hf(e->twin());
}
/** \brief Check the segment attribute of a halfedge
 * \param e : handle of the halfedge to be checked.
 * \param verbose : if true, details are sent to std::clog. Default to false.
 * \return true if the segment  is valid, false otherwise.
 *
 * Features checked: alignement of the halfedge and its segment,
 * occurence of the segment in the tessellation segment list.
 */
bool LineTes::check_halfedge_segment(LineTes::Halfedge_handle e,
				     bool verbose) {
  Seg_handle s = e->get_segment();
  if (s==NULL_SEG_HANDLE) {
    if (verbose) {
      std::clog << "LineTes::check_halfdege_segment: halfedge has a null ";
      std::clog << "segment" << std::endl;
    }
    return false;
  }
  Point2 s0 = s->pointSource();
  Point2 s1 = s->pointTarget();
  // check e is collinear with the edge representing its segment
  if (!CGAL::collinear(s0,s1,e->source()->point()) ||
      !CGAL::collinear(s0,s1,e->source()->point())) {
    if (verbose) {
      std::clog << "LineTes::check_halfdege_segment: halfedge not aligned ";
      std::clog << "with its segment" << std::endl;
    }
    return false;
  }
  // check segment is in the segment list
  if (std::find(segments_begin(),segments_end(),s)==segments_end()) {
    if (verbose) {
      std::clog << "LineTes::check_halfdege_segment: halfedge segment ";
      std::clog << "not registered" << std::endl;
    }
    return false;
  }
  return true;
}
/** \brief Check the direction attribute of a halfedge
 * \param e : handle of the halfedge to be checked.
 * \param verbose : if true, details are sent to std::clog. Default to false.
 * \return true if the halfedge and its segment directions are consistent, 
 * false otherwise.
 */
bool LineTes::check_halfedge_dir(LineTes::Halfedge_handle e, bool verbose) {
  Vector ve(e->source()->point(),e->target()->point());
  Vector vs(e->segment()->pointSource(),e->segment()->pointTarget());
  if (CGAL::angle(ve,vs)!=CGAL::OBTUSE) {
    if (verbose) {
      std::clog << "LineTes::check_halfedge_dir: halfedge and segment ";
      std::clog << "directions not consistent" << std::endl;
    }
    return false;
  }
  return true;
}
/** \brief Return the number of sides of the tessellated domain
 */
LineTes::Size LineTes::number_of_window_edges() {
  return window.size();
}
/** \brief Return the perimeter of the tessellated domain
 */
double LineTes::get_window_perimeter() {
  double per=0.0;
  for(Polygon::Edge_const_iterator e=window.edges_begin();e!=window.edges_end();e++) {
    Vector v = (*e)[1]-(*e)[0];
    per +=  sqrt(CGAL::to_double(v.x()*v.x()+v.y()*v.y()));
  }
  return per;
}
/** \brief Print on standard output the segment points
 * \param endsOnly : see the documentation of the print method of the Seg class.
 */
void LineTes::print_all_segments(bool endsOnly) {
  for (Seg_list_iterator si=segments_begin();si!=segments_end();si++) {
    (*si)->print(endsOnly);
  }
}
/** \brief Test whether a halfedge lies along the domain boundary
 *
 * Test if an edge lies on the window boundary. Return
 * - 0 if the edge is internal.
 * - 1 if the halfedge bounds an internal (bounded) face.
 * - 2 if the halfedge bounds the external (unbounded) face.
 */
int LineTes::is_on_boundary(Halfedge_handle e){
  Face_handle f_e,f_e_twin;

  f_e = e->face();
  f_e_twin = e->twin()->face();

  if (!f_e->is_unbounded() && !f_e_twin->is_unbounded()) 
    return 0;
  else { // One of the faces bounded by the edge is unbounded
    if (!f_e->is_unbounded()) 
      return 1;
    else 
      return 2;
  }
}


LineTes::Halfedge_handle LineTes::split_edge(Halfedge_handle e,
					     X_monotone_curve_2 c1,
					     X_monotone_curve_2 c2) {
  /* Invoke the method split_edge of the class Planar_map and update
     the edge junctions, the dir and length fields of both new
     edges and the segment list. */

  Halfedge_handle e1, e2, e_next_hf;

  // Get next neighbour of e before split
  e_next_hf = e->get_next_hf();
  // Get segment supporting e before split
  Seg_handle s = e->segment();
  if (s->get_halfedge_handle()==e) {
    s->set_halfedge_handle(NULL_HALFEDGE_HANDLE);
  }
  // Get e direction before split
  bool e_dir = e->get_dir();

  // Do the split
  Halfedge_handle k2 = e;
  e1 = Arrangement::split_edge(e,c1,c2);
  // We believe that e, e1 and e2 have the same "directions".
  e2 = e1->next();
  
  // Set junctions
  set_junction(e1,e2);
  set_junction(e2,e_next_hf);

  // Set directions
  e1->set_dir(e_dir);
  e1->twin()->set_dir(!e1->get_dir());
  e2->set_dir(e1->get_dir());
  e2->twin()->set_dir(!e2->get_dir());

  // Set length
  e1->set_length();
  e2->set_length();
  e1->twin()->set_length(e1->get_length());
  e2->twin()->set_length(e2->get_length());

  // Set segment
  e1->set_segment(s); e1->twin()->set_segment(s);
  e2->set_segment(s); e2->twin()->set_segment(s);

  /* Update segment. */
  s->set_halfedge_handle(e1);

  return e1;
}


LineTes::Halfedge_handle LineTes::merge_edge(LineTes::Halfedge_handle e1,
					     LineTes::Halfedge_handle e2,
					     X_monotone_curve_2 cv) {
  /* Invoke the method merge_edge of the class Planar_map and update
     the edge junctions, the dir and length fields of the merged edge
     and the segment list. */

  Halfedge_handle hf,hf_avant,hf_apres;
  double new_length;
  
  // Get halfedges before and after e1,e2 and the sum of their length
  hf_avant = e1->get_prev_hf();
  hf_apres = e2->get_next_hf();
  new_length = e1->get_length()+e2->get_length();
  // Get segment supporting e1 and e2
  Seg_handle s = e1->segment();

  // Merge the two halfedges
  hf = Arrangement::merge_edge(e1,e2,cv); 
  
  // Set junctions
  set_junction(hf_avant,hf);
  set_junction(hf,hf_apres);
  
  // Set length
  hf->set_length(new_length);
  hf->twin()->set_length(hf->get_length()); 
 
  // Set directions
  hf->set_dir(e1->get_dir());
  hf->twin()->set_dir(!hf->get_dir());

  // Set segment
  hf->set_segment(s);
  hf->twin()->set_segment(s);

  // Update segment. 
  s->set_halfedge_handle(hf);
  
  return hf;
}

LineTes::Halfedge_handle LineTes::add_edge(X_monotone_curve_2 c,
					   Vertex_handle v1, 
					   Vertex_handle v2) {
  /* Invoke the method insert_at_vertices of the class Planar_map and
     update the edge junctions and set the fields dir and length of
     the new edge. */

  // Insert new edge
  Halfedge_handle e = insert_at_vertices(c,v1,v2);

  // Set length
  e->set_length();
  e->twin()->set_length(e->get_length());
  
  // Set junctions
  set_junction(e,compute_next_hf(e));
  set_junction(compute_prev_hf(e),e);

  // Set directions
  if (e->get_prev_hf()!=NULL_HALFEDGE_HANDLE || e->get_next_hf()!=NULL_HALFEDGE_HANDLE)
    e->set_dir();
  else
    e->set_dir(true);
  e->twin()->set_dir(!e->get_dir());

  /* Update segment list. Note that either an new segment is added or
     the halfedge e is added at one end of a segment. */

  Seg_handle s;
  if (e->get_next_hf()==NULL_HALFEDGE_HANDLE && e->get_prev_hf()==NULL_HALFEDGE_HANDLE) {
    s = new Seg;
    s->set_halfedge_handle(e);
    all_segments.add(s);
  } else if (e->get_next_hf()!=NULL_HALFEDGE_HANDLE)
    s = e->get_next_hf()->segment();
  else
    s = e->get_prev_hf()->segment();
  e->set_segment(s); e->twin()->set_segment(s);

  return e;
}


LineTes::Face_handle LineTes::remove_edge(Halfedge_handle e){
  /* Remove an edge from the planar make (invocation of the parent
     class method and update edge junctions). Return the
     handle to the face which has been merged.

     Condition : e should not lie on the window boundary. In such a
     case, the planar map is not modified and NULL_FACE_HANDLE is
     returned. 

     Note: vertices at the ends of the removed edge are not
     removed. If needed, removal of spurious vertices has to be
     performed after having invoked this method.*/

  Face_handle f;

  if (is_on_boundary(e)!=0) 
    return NULL_FACE_HANDLE;

  // Update junctions
  set_junction(e->get_prev_hf(),NULL_HALFEDGE_HANDLE);
  set_junction(NULL_HALFEDGE_HANDLE,e->get_next_hf());

  // Get segment before removal
  Seg_handle s = e->segment();
  
  // Update segment or segment list
  if (e->get_next_hf()==NULL_HALFEDGE_HANDLE && e->get_prev_hf()==NULL_HALFEDGE_HANDLE)
    all_segments.suppress(s);
  else if (e->get_prev_hf()!=NULL_HALFEDGE_HANDLE)
    s->set_halfedge_handle(e->get_prev_hf());
  else
    s->set_halfedge_handle(e->get_next_hf());

  // Remove the edge
  f = Arrangement::remove_edge(e,false,false);

  return f;
  }


/******************************************************************************/
/*                METHODS FOR THE CLASS LineTes::Seg                          */
/******************************************************************************/

/** \brief Empty constructor*/
LineTes::Seg::Seg() {}

/** \brief Return a Segment object identified by the provided halfedge handle*/
LineTes::Seg::Seg(Halfedge_handle hh): e(hh) {}

/** \brief Return the halfedge handle starting the segment
 *
 * The segment orientation is given by its identifying halfedge handle.
 */
LineTes::Halfedge_handle LineTes::Seg::halfedges_start() {
  Halfedge_handle e_curr = e;
  while(e_curr->get_prev_hf()!=NULL_HALFEDGE_HANDLE)
    e_curr = e_curr->get_prev_hf();
  return e_curr;
}

/** \brief Return the halfedge handle ending the segment*/
LineTes::Halfedge_handle LineTes::Seg::halfedges_end() {
  Halfedge_handle e_curr = e;
  while(e_curr->get_next_hf()!=NULL_HALFEDGE_HANDLE)
    e_curr = e_curr->get_next_hf();
  return e_curr;
}

/** \brief Return the number of edges along the segment*/
LineTes::Size LineTes::Seg::number_of_edges() {
  Halfedge_handle e = halfedges_start();
  Size i=1;
  while (e->get_next_hf()!=NULL_HALFEDGE_HANDLE) {
    e = e->get_next_hf();
    i++;
  }
  return i;
}

/** \brief Test on the number of edges of the segment*/
bool LineTes::Seg::number_of_edges_is_greater_than(unsigned int n) {
  Halfedge_handle e = halfedges_start();
  unsigned int i=1;
  while (e->get_next_hf()!=NULL_HALFEDGE_HANDLE) {
    i++;
    e = e->get_next_hf();
    if (i>n)
      return true;
  }
  return false;
}

/** \brief Set the halfedge handle identifying a segment
 *
 * If the direction attribute of the halfedged handle is false, the
 * opposite halfedge is used instead.
*/
void LineTes::Seg::set_halfedge_handle(Halfedge_handle hh) {
  if (hh==NULL_HALFEDGE_HANDLE) {
    e = hh;
    return;
  }
  if (hh->get_dir()) {
    e = hh;
  } else {
    e = hh->twin();
  }
}

/** \brief Return the point starting the segment*/
Point2 LineTes::Seg::pointSource() {
  return halfedges_start()->source()->point();
}

/** \brief Return the point ending the segment*/
Point2 LineTes::Seg::pointTarget() {
  return halfedges_end()->target()->point();
}

/** \brief Return the list of vertices (as points) lying along the segment*/
Points LineTes::Seg::list_of_points() {
  std::vector<Point2> lp;
  Halfedge_handle e = halfedges_start();
  lp.push_back(e->source()->point());
  while (e!=NULL_HALFEDGE_HANDLE) {
    lp.push_back(e->target()->point());
    e = e->get_next_hf();
  }
  return lp;
}

/** \brief Print on standard output the segment points
 *
 * \param endsOnly : if true, only endpoints are printed. Otherwise,
 * all vertices along the segments are printed.
 */
void LineTes::Seg::print(bool endsOnly) {
  Halfedge_handle e = halfedges_start();
  std::cout << "(" << CGAL::to_double (e->source()->point()[0]) << ", " 
	    << CGAL::to_double( e->source()->point()[1]) << ") - ";
  while (e->get_next_hf()!= NULL_HALFEDGE_HANDLE){
    if (!endsOnly) {
      std::cout << "(" << CGAL::to_double(e->target()->point()[0]) << ", " 
		<< CGAL::to_double( e->target()->point()[1]) << ") - ";
    }
    e=e->get_next_hf();
  }
  std::cout << "(" << CGAL::to_double (e->target()->point()[0]) << ", " 
	    << CGAL::to_double( e->target()->point()[1]) << ") ";
  std::cout<<std::endl;
}

/** \brief Insert a new line segment in a line tessellation
 *
 * \param iseg : line segment to be inserted in the tessellation.
 * \return handle to the newly created Seg object.
 *
 * Eventually the segment is clipped by the tessellation domain before
 * insertion.
 */
LineTes::Seg_handle LineTes::insert_segment(Segment iseg) {
  typedef std::list<Halfedge_handle> Edges;
  // Clip the segment to the window
  /* The window is a polygon. There is no way with CGAL to compute the
     intersection of a segment with a polygon even if the polygon is
     convex. The current approach is to convert the window into a
     Rectangle object. Therefore our implementation works only when
     the window is indeed an iso-rectangle. */
  Polygon w = get_window();
  /* CGAL::Bbox_2 bbox = w.bbox();
  Rectangle rw(bbox.xmin(),bbox.ymin(),bbox.xmax(),bbox.ymax()); // window as an iso-rectangle
  Segment cseg; // clipped segment
  CGAL::Object inter = CGAL::intersection(rw,iseg);
  if (!CGAL::assign(cseg,inter)) {
    std::cerr << "Insertion of segment, error when clipping segment to window";
    std::cerr << "Window " << rw << std::endl;
    std::cerr << "Segment " << iseg << std::endl;
    std::cerr << std::endl;
    } */
  Segment cseg = clip_segment_by_convex_polygon(iseg,w);
  CGAL::insert(*this,cseg); // insert segment into the arrangement
  // identify new edges generated by the segment insertion
  Edges new_edges;
  for (Edge_iterator e=this->edges_begin();e!=this->edges_end();e++) {
    if (e->get_length()==0.0)
      new_edges.push_back(e);
  }
  // set length and neighbours of the new edges
  for (Edges::iterator ei=new_edges.begin();ei!=new_edges.end();ei++) {
    Halfedge_handle e = *ei;
    e->set_length();
    e->twin()->set_length(e->get_length());
    set_junction(e,compute_next_hf(e));
    set_junction(compute_prev_hf(e),e);
  }
  /* There are also split edges with non-updated length, see below. */
  // set segment and dir of the new edges
  Seg_handle s = new Seg;
  while (new_edges.size()>0) {
    Halfedge_handle e = *(new_edges.begin());
    if ((e->get_next_hf()!=NULL_HALFEDGE_HANDLE && 
	 (e->get_next_hf()->get_segment()!=NULL_SEG_HANDLE)) ||
	(e->get_prev_hf()!=NULL_HALFEDGE_HANDLE &&
	 e->get_prev_hf()->get_segment()!=NULL_SEG_HANDLE)) {
      // halfedge comes from a split of a former edge
      e->set_dir();
      e->twin()->set_dir(!(e->get_dir()));
      // halfedge has at least a neighbour along its segment
      if (e->get_prev_hf()!=NULL_HALFEDGE_HANDLE) {
	e->set_segment(e->get_prev_hf()->get_segment());
      } else {
	e->set_segment(e->get_next_hf()->get_segment());
      }
      // Update neighbour lengths
      if (e->get_prev_hf()!=NULL_HALFEDGE_HANDLE) {
	e->get_prev_hf()->set_length();
	e->get_prev_hf()->twin()->set_length(e->get_prev_hf()->get_length());
      }
      if (e->get_next_hf()!=NULL_HALFEDGE_HANDLE) {
	e->get_next_hf()->set_length();
	e->get_next_hf()->twin()->set_length(e->get_next_hf()->get_length());
      }
      e->twin()->set_segment(e->get_segment());
      new_edges.remove(e);
    } else {
      // halfedge e lies on the inserted segment
      e->set_dir(true); e->twin()->set_dir(false);
      s->set_halfedge_handle(e);
      e->set_segment(s); e->twin()->set_segment(s);
      all_segments.add(s);
      new_edges.remove(e);
      Halfedge_handle es = e->get_next_hf();
      while (es!=NULL_HALFEDGE_HANDLE) {
	es->set_segment(s); es->twin()->set_segment(s);
	es->set_dir(true); es->twin()->set_dir(false);
	new_edges.remove(es);
	es = es->get_next_hf();
      }
      es = e->get_prev_hf();
      while (es!=NULL_HALFEDGE_HANDLE) {
	es->set_segment(s); es->twin()->set_segment(s);
	es->set_dir(true); es->twin()->set_dir(false);
	new_edges.remove(es);
	es = es->get_prev_hf();
      }
    }
  }
  return s;
}

/** \brief Remove all I-vertices from a line tessellation
 * \param imax : maximal number of I-vertices to be removed. If zero, 
 * (default value) all I-vertices are removed.
 * \param verbose : If true, information about I-vertex processing is 
 * sent to the standard output. Default to false.
 *
 * An I-vertex is a tessellation vertex of degree 1. Removal is
 * performed either by shortening or lengthening of the incident
 * edge. Shortening is preferred if the removed length is less than
 * the added length.
 *
 * The sequence of removals is ordered by length variation. At each
 * iteration, the I-vertex to be removed is the one inducing the
 * smallest length variation.
 */
void LineTes::remove_ivertices(Size imax, bool verbose) {
  if (imax==0)
    imax = std::numeric_limits<Size>::max();
  double total_removed_length = 0.0, total_added_length = 0.0;
  double total_internal_length = 0.0;
  if (verbose) { // compute internal length
    for (Seg_list_iterator s=segments_begin();s!=segments_end();s++) {
      if (!is_on_boundary(*s)) { // internal segment
	NT l2= CGAL::squared_distance((*s)->pointSource(),
				      (*s)->pointTarget());
	total_internal_length += sqrt(CGAL::to_double(l2));
      }
    }
  }
  Size counter = 1;
  while (true) { // loop on vertice removal
    // Find I-vertices and length changes induced by their removal
    std::vector<IVertex> ivertices;
    for (Vertex_iterator v=vertices_begin();v!=vertices_end();v++) {
      if (v->degree()==1) {
	IVertex ivertex;
	ivertex.v = v;
	// find incident edge
	Halfedge_handle e = v->incident_halfedges();
	ivertex.e = e;
	double lm = e->get_length(); // removed length
	// Computation of the added length
	Rayon r(v->point(),Vector(e->source()->point(),e->target()->point()));
	// find the halfedge hit by ray r and the intersection
	/* Not that easy. The target halfedge may be on the outer ccb
	   or on a "hole" boundary. Approach below: visit all
	   halfedges (on ccb or along holes), compute intersection, if
	   any keep the closest one. */
	Face_handle f = e->face();
	std::vector<Halfedge_handle> all_edges;
	Ccb_halfedge_circulator hc=f->outer_ccb();
	Halfedge_handle h0 = hc;
	do {
	  if (hc->source()!=v && hc->target()!=v)
	    all_edges.push_back(hc);
	  hc++;
	} while (hc!=h0);
	for (Hole_iterator ho=f->holes_begin();ho!=f->holes_end();ho++) {
	  hc = *ho;
	  h0 = hc;
	  do {
	  if (hc->source()!=v && hc->target()!=v)
	    all_edges.push_back(hc);
	  hc++;
	  } while (hc!=h0);
	}
	double lp = std::numeric_limits<double>::max();
	Point2 inter_location;
	for (std::vector<Halfedge_handle>::iterator e=all_edges.begin();e!=all_edges.end();e++) {
	  CGAL::Object inter = CGAL::intersection(r,
	 					  Segment((*e)->source()->point(),
	 						  (*e)->target()->point())); 
	  if (CGAL::assign(inter_location,inter)) {
	    double len = sqrt(CGAL::to_double(CGAL::squared_distance(v->point(),
								inter_location)));
	    if (len<lp) {
	      ivertex.p = inter_location;
	      ivertex.esplit = *e;
	      lp = len;
	    }
	  }
	}
	if (lm<=lp) {
	  ivertex.changed_length = lm;
	  ivertex.shorten = true;
	} else {
	  ivertex.changed_length = lp;
	  ivertex.shorten = false;
	}
	ivertices.push_back(ivertex);
      } // end test degree is 1
    } // end loop on all vertices
    if (ivertices.size()==0) {
      if (verbose) 
	std::cout << "no I-vertex left" << std::endl;
      counter--;
      break;
    }
    // Remove the I-vertex with smallest length change
    std::vector<IVertex>::iterator ivertex = min_element(ivertices.begin(),
							 ivertices.end(),
							compare_ivertices);
    if (verbose) {
      std::cout << "iteration " << counter << std::endl;
      std::cout << "processing I-vertex " << ivertex->v->point() << std::endl;
      std::cout << "length variation " << ivertex->changed_length << std::endl;
      if (ivertex->shorten) {
	std::cout << "shortening " << std::endl;
      } else {
	std::cout << "lengthening " << std::endl;
      }
    }
    if (ivertex->shorten) { // shortening is preferred
      suppress_edge(ivertex->e);
      total_removed_length += ivertex->changed_length;
    } else {               // lengthening is preferred
      split_from_vertex(ivertex->v,ivertex->esplit,ivertex->p);
      total_added_length += ivertex->changed_length;
    } // end alternative shortening/lengthening
    if (counter>imax) {
      if (verbose) {
	std::cout << "maximal number of iterations reached" << std::endl;
	std::cout << "still " << ivertices.size() << " I-vertices left";
	std::cout << std::endl;
      }
      break;
    } else {
      counter++;
    }
  } // end loop on I-vertex removal
  if (verbose) {
    std::cout << std::endl;
    std::cout << counter << " I-vertices removed" << std::endl;
    std::cout << "initial total (internal) length " << total_internal_length;
    std::cout << std::endl;
    std::cout << "removed length " << total_removed_length << std::endl;
    std::cout << "added length " << total_added_length << std::endl;
  }
} // end method

/** \brief Remove all L-vertices from a line tessellation
 * \param imax : maximal number of I-vertices to be removed. If zero, 
 * (default value) all I-vertices are removed.
 * \param verbose : if true, information about I-vertex processing is 
 * sent to the standard output. Default to false.
 *
 * An L-vertex is a tessellation vertex of degree 2 with both incident
 * edges which are not aligned. Removal is performed by
 * lengthening one of the incident edges. The edge to be lengthened is the one that induces the shortest length variation.
 *
 * The sequence of removals is ordered by length variation. At each
 * iteration, the L-vertex to be removed is the one inducing the
 * smallest length variation.
 */
void LineTes::remove_lvertices(Size imax, bool verbose) {
  if (imax==0)
    imax = std::numeric_limits<Size>::max();
  double total_added_length = 0.0;
  double total_internal_length = 0.0;
  if (verbose) { // compute internal length
    for (Seg_list_iterator s=segments_begin();s!=segments_end();s++) {
      if (!is_on_boundary(*s)) { // internal segment
	NT l2= CGAL::squared_distance((*s)->pointSource(),
				      (*s)->pointTarget());
	total_internal_length += sqrt(CGAL::to_double(l2));
      }
    }
  }
  Size counter = 1;
  while (true) { // loop on vertice removal
    // Find L-vertices and length changes induced by their removal
    // We use the IVertex structure for storing data about L-vertices
    // Field shorten is not used.
    std::vector<IVertex> lvertices;
    for (Vertex_iterator v=vertices_begin();v!=vertices_end();v++) {
      if (v->degree()==2 && is_point_inside_window(v->point(),this)>0) {
	IVertex lvertex1, lvertex2;
	Halfedge_around_vertex_circulator e1 = v->incident_halfedges(), e2 = e1;
	e2++;
	lvertex1.v = e1->target();
	lvertex1.changed_length = precompute_lengthening(e1,&(lvertex1.esplit),
						       &lvertex1.p);
	lvertex2.v = e2->target();
	lvertex2.changed_length = precompute_lengthening(e2,&(lvertex2.esplit),
						       &lvertex2.p);
	if (compare_ivertices(lvertex1,lvertex2)) {
	  lvertices.push_back(lvertex1);
	} else {
	  lvertices.push_back(lvertex2);
	}
      } // end test degree is 2
    } // end loop on all vertices
    if (lvertices.size()==0) {
      if (verbose) 
	std::cout << "no L-vertex left" << std::endl;
      counter--;
      break;
    }
    // Remove the L-vertex with smallest length change
    std::vector<IVertex>::iterator lvertex = min_element(lvertices.begin(),
							 lvertices.end(),
							 compare_ivertices);
    if (verbose) {
      std::cout << "iteration " << counter << std::endl;
      std::cout << "processing L-vertex " << lvertex->v->point() << std::endl;
      std::cout << "length variation " << lvertex->changed_length << std::endl;
    }
    split_from_vertex(lvertex->v,lvertex->esplit,lvertex->p);
    total_added_length += lvertex->changed_length;
    if (counter>imax) {
      if (verbose) {
	std::cout << "maximal number of iterations reached" << std::endl;
	std::cout << "still " << lvertices.size() << " L-vertices left";
	std::cout << std::endl;
      }
      break;
    } else {
      counter++;
    }
  } // end loop on L-vertex removal
  if (verbose) {
    std::cout << std::endl;
    std::cout << counter << " L-vertices removed" << std::endl;
    std::cout << "initial total (internal) length " << total_internal_length;
    std::cout << std::endl;
    std::cout << "added length " << total_added_length << std::endl;
  }
} // end method

/** \brief Test whether a line tessellation is a T-tessellation
 * \param verbose : if true, details are sent to std::clog. Default to false.
 * \param out : output stream to be written if verbose is true.
 * \return true if is a T-tessellation, false otherwise.
 * 
 * Check whether all internal vertices are T-vertices and whether all
 * vertices on the domain boundary are either of degree 2 or
 * T-vertices.
 */
bool LineTes::is_a_T_tessellation(bool verbose, std::ostream& out) {
  for (Vertex_iterator v=vertices_begin();v!=vertices_end();v++) {
    if (is_point_inside_window(v->point(),this)) {
      if (!is_a_T_vertex(v)) {
	if (verbose) {
	  out << "Internal vertex " << v->point() << " is not T";
	  out << std::endl;
	}
	return false;
      }
    } else { // vertex on the window boundary
      if (v->degree()!=2 && !is_a_T_vertex(v)) {
	if (verbose) {
	  out << "Boundary vertex " << v->point();
	  out << " is not of degree 2";
	  out << " and is not a T-vertex" << std::endl;
	}
	return false;
      }
    }
  } // end loop on vertices
  return true;
}

/** \brief Write a LineTes object to an output stream
 * \param os : output stream.
 *
 * Output is done based on a format without loss of numerical
 * precision. Therefore it is possible to recover a LineTes object
 * from the output. The output format is as follows:
 * - A coordinate is output as two integers (numerator and denominator).
 * - A point is output as its two Cartesian coordinates.
 * - A line segment is output as its two end points.
 * - First, vertices of the tessellated domain are sent to the output 
 *   stream (single line).
 * - Subsequently, internal segments are sent to the output stream
 *   (one per line).
 */	
void LineTes::write(std::ostream& os) {
  Polygon w = get_window();
  for (Size i=0;i<w.size();i++) {
    os << w[i].x().exact() << " " << w[i].y().exact();
    if (i!=(w.size()-1))
      os << " ";
  }
  os << std::endl;
  for (Seg_list_iterator si=segments_begin();si!=segments_end();si++) {
    Seg_handle s = *si;
    if (!(is_on_boundary(s))) {
      os << s->pointSource().x().exact() << " ";
      os << s->pointSource().y().exact() << " ";
      os << s->pointTarget().x().exact() << " ";
      os << s->pointTarget().y().exact() << std::endl;
    }
  }
}
/** \brief Read a LineTes object from an input stream
 * \param is : input stream.
 *
 * Data provided by the input stream should be formatted as described in
 * LineTes::write documentation.
 */
void LineTes::read(std::istream& is) {
  std::string line;
  // Read the coordinates of the domain corners
  if(!std::getline(is,line)) {
    std::cerr << "Could not read the domain coordinates" << std::endl;
    std::cerr << "Reading aborted" << std::endl;
    return;
  }
  std::istringstream iss(line);
  std::vector<NT> wcoords;
  while(true) {
    CGAL::Gmpq buf;
    if (iss>>buf) {
      NT coord(buf);
      wcoords.push_back(coord);
    } else {
      break;
    }
  }
  // Now read segments
  std::vector<std::vector<CGAL::Gmpq> > scoords;
  while (std::getline(is,line)) {
    std::istringstream iss_bis(line);
    CGAL::Gmpq x0, y0, x1, y1;
    std::vector<CGAL::Gmpq> buf(4,0);
    if (iss_bis >> buf[0] >> buf[1] >> buf[2] >> buf[3]) {
      scoords.push_back(buf);
    } else {
      std::cerr << "Unexpected input when reading a line tessellation";
      std::cerr << std::endl;
      std::cerr << "Input line" << std::endl;
      std::cerr << line << std::endl;
      std::cerr << "Reading aborted" << std::endl;
      return;
    }
  }
  // Fill the LineTes object
  clear();
  Polygon w;
  if (wcoords.size()%2 != 0) {
    std::cerr << "Uneven number of coordinates on the first line" << std::endl;
    return;
  }
  for (Size i=0;i!=wcoords.size();i+=2) {
    w.push_back(Point2(wcoords[i],wcoords[i+1]));
  }
  insert_window(w);
  for (Size i=0;i!=scoords.size();i++) {
    NT x0(scoords[i][0]);
    NT y0(scoords[i][1]);
    NT x1(scoords[i][2]);
    NT y1(scoords[i][3]);
    Segment s(Point2(x0,y0),Point2(x1,y1));
    insert_segment(s);
  }
}
/******************************************************************************/
/*                METHODS FOR THE BASE CLASS LineTes_halfedge                 */
/******************************************************************************/

/** \brief Basic constructor
 *
 * Invoke the parent constructor, set length to zero and next and previous
 * halfedges along the segment as null halfedges.
 */
LineTes_halfedge::LineTes_halfedge() : CGAL::Arr_halfedge_base<Curve>(),
				       length(0),
				       next_hf(NULL_HALFEDGE_HANDLE),
				       prev_hf(NULL_HALFEDGE_HANDLE),
                                       segment_ptr(NULL_SEG_HANDLE) {}
 
/** \brief Update the length attribute of a halfedge
 */
void LineTes_halfedge::set_length() {
  CGAL::Vector_2<Kernel> v(curve());
  length = sqrt(CGAL::to_double(v.x()*v.x()+v.y()*v.y()));
}
/** \brief Set the next halfedge along the segment
 */  
void LineTes_halfedge::set_next_hf(LineTes::Halfedge_handle e) {
  next_hf = e;
}
/** \brief Set the previous halfedge along the segment
 */
void LineTes_halfedge::set_prev_hf(LineTes::Halfedge_handle e) {
  prev_hf = e;
}
/** \brief Update the direction attribute of a halfedge
 * \pre The halfedge should have either a next or a previous halfedge on 
 * the segment.
 *
 * The direction is inherited from the neighbour halfedges on the segment. 
 * Priority given to the next halfedge. 
 */
void LineTes_halfedge::set_dir() {
  if (next_hf!= NULL_HALFEDGE_HANDLE) 
    dir = next_hf->get_dir();
  else if (prev_hf!= NULL_HALFEDGE_HANDLE)
    dir = prev_hf->get_dir();
  else
    std::clog << "LineTes_halfedge::set_dir : no neighbour found" 
	      << std::endl;
}
/** \brief Set the direction attribute of a halfedge
 *
 * The convention is to set the direction to true if the halfedge and 
 * the halfedge identifying the segment have the same direction.
 */
void LineTes_halfedge::set_dir(bool d) {
  dir = d;
}




/******************************************************************************/
/*                     METHODS FOR THE CLASS TTessel                          */
/******************************************************************************/

/** \brief Empty constructor
 *
 * Invoke the empty constructor of the parent class LineTes.
 */
TTessel::TTessel() : LineTes() {}

/** \brief Generate a TTessel object from a LineTes Object
 * \param lt : the LineTes object to copy
 *
 * The LineTes object should be a valid T-tessellation.
 */
TTessel::TTessel(LineTes& lt) {
  TTessel();
  if (lt.is_a_T_tessellation()) {
    Polygon p = lt.get_window();
    insert_window(p);
    for (Seg_list_iterator s=lt.segments_begin();s!=lt.segments_end();s++)
      insert_segment(Curve((*s)->pointSource(),(*s)->pointTarget()));
    int_length = 0.0;
    for (Seg_list_iterator s=segments_begin();s!=segments_end();s++) {
      double seglen = 
	sqrt(CGAL::to_double(CGAL::squared_distance((*s)->pointSource(),
						    (*s)->pointTarget())));
      if (is_on_boundary(*s)) { 
	int_length += seglen;
	continue;
      } else
	int_length += 2*seglen;
	
      if ((*s)->number_of_edges()==1) {
	non_blocking_segments.add(*s);
      } else {
	blocking_segments.add(*s);
      }
    }
  } else {
    std::cout << "warning: input is not a T-tesselation, returning a void TTessel object" << std::endl;
  }
}

/** \brief Define the rectangular domain to be tessellated.
 */
void TTessel::insert_window(Rectangle r) {
  LineTes::insert_window(r);
  int_length = 0;
  for (int i=0;i<4;i++) {
    Vector v(r[i],r[i+1]);
    int_length += sqrt(CGAL::to_double(v.x()*v.x()+v.y()*v.y()));
  }
}

/** \brief Define the polygonal domain to be tessellated
 *
 * Do not use this method at the moment as non-rectangular domains
 * are not fully supported by the LineTes class.
 */
void TTessel::insert_window(Polygon& p) {
  LineTes::insert_window(p);
  int_length = 0;
  for (Polygon::Edge_const_iterator e=p.edges_begin();e!=p.edges_end();e++) 
    int_length += sqrt(CGAL::to_double(e->squared_length()));
}

/** \brief Split a T-tessellation
 *
 * A cell of the current tessellation is split by a new
 * edge.
 */
TTessel::Halfedge_handle TTessel::update(Split spl) {
  // Do the split.
  Halfedge_handle e = split_from_edge(spl.get_e1(),spl.get_p1(),spl.get_e2(),
				      spl.get_p2());
  // Update info
  int_length += 2*e->get_length();
  non_blocking_segments.add(e->segment());
  Seg_handle s = e->next()->segment();
  if (!is_on_boundary(s) && !s->number_of_edges_is_greater_than(2)) {
    non_blocking_segments.suppress(s);
    blocking_segments.add(s);
  }
  s = e->twin()->next()->segment();
  if (!is_on_boundary(s) && !s->number_of_edges_is_greater_than(2)) {
    non_blocking_segments.suppress(s);
    blocking_segments.add(s);
  }
  return e;
}

/** \brief Merge a T-tessellation
 *
 * Two cells separated by a non-blocking segment are merged
 * by the deletion of that segment.
 */
TTessel::Face_handle TTessel::update(Merge me) {

  Face_handle f;
  Seg_handle s,s1,s2 ;
  Halfedge_handle e;

  e = me.get_e();
  s = e->segment();

  // Update blocking and non-blocking segment lists before merge

  non_blocking_segments.suppress(s);

  s1 = e->next()->segment();
  if (!is_on_boundary(s1) && !s1->number_of_edges_is_greater_than(2)) { 
    non_blocking_segments.add(s1);  
    blocking_segments.suppress(s1);
  }

  s2 = e->twin()->next()->segment();
  if (!is_on_boundary(s2) && !s2->number_of_edges_is_greater_than(2)) {
    non_blocking_segments.add(s2);  
    blocking_segments.suppress(s2);
  }
  // Get length of suppressed edge before merge
  double suppressed_length = e->get_length();

  // Do the merge
  f = suppress_edge(e);
  
  // Update internal length
  int_length -= 2*suppressed_length;

  return f;
}

/** \brief Flip a T-tessellation
 *
 * A segment is shortened and the segment blocked by the former is lengthened (???).
 */
TTessel::Halfedge_handle TTessel::update(Flip flip) {

  Halfedge_handle e1 = flip.get_e1(); // edge to be suppressed
  Seg_handle s1 = e1->segment(); // segment to be shortened
  Seg_handle s1_i = e1->next()->segment(); /* segment which
						       blocks the
						       segment to be
						       shortened */
  // Find the edge which must be prolonged
  Halfedge_around_vertex_circulator e2 = e1->source()->incident_halfedges();
  while (e2->segment()==s1)
    e2++;
  Seg_handle s2 = e2->segment(); // segment to be prolonged

  // Prolong the segment and get the new edge
  Halfedge_handle new_e = split_from_vertex(e1->source(),flip.get_e2(),
					    flip.get_p2());
  int_length += 2*new_e->get_length();
  int_length -= 2*e1->get_length();
  // Get the segment incident to the prolonged segment
  Seg_handle s2_i = new_e->next()->segment();

  // Shorten the segment
  suppress_edge(e1);

  // Update the lists of blocking and non-blocking segments

  /* Prolonged segment s2. If s2 has only two edges, s2 was
     non-blocking before the flip and it is blocking now. */
  if (!s2->number_of_edges_is_greater_than(2)) {
    non_blocking_segments.suppress(s2);
    blocking_segments.add(s2);
  }
  /* Shortened segment s1. If s1 has only one edge, s1 is non-blocking
     now and it was blocking before the flip. */
  if (!s1->number_of_edges_is_greater_than(1)) {
    blocking_segments.suppress(s1);
    non_blocking_segments.add(s1);
  }
  /* If the segment incident to s1 before the flip is incident to s2
     after the flip, its status does not change. */
  if (s2_i!=s1_i) {
    /* Segment incident to s2. If s2_i has only two edges, s2_i was
       non-blocking before the flip and it is blocking now. */
    if (!is_on_boundary(s2_i) && !s2_i->number_of_edges_is_greater_than(2)) {
      non_blocking_segments.suppress(s2_i);
      blocking_segments.add(s2_i);
    }
    /* Segment incident to s1. If s1_i has only one edge, s1_i is non-blocking
       now and it was blocking before the flip. */
    if (!is_on_boundary(s1_i) && !s1_i->number_of_edges_is_greater_than(1)) {
      blocking_segments.suppress(s1_i);
      non_blocking_segments.add(s1_i);
    }
  }

  return new_e;
}
/** \brief Empty the T-tessellation
 */
void TTessel::clear() {
  LineTes::clear();
  blocking_segments.clear();
  non_blocking_segments.clear();
}
/** \brief Check validity of a T-tessellation
 * param verbose : if true, details are sent to std::clog. Default to false.
 * \return true if valid, false otherwise.
 *
 * The following checks are performed:
 * - Validity as a LineTes object.
 * - All internal vertices are T-vertices.
 * - All non-blocking segments consist of one edge only.
 * - All blocking segments consist of two edges at least.
 */
bool TTessel::is_valid(bool verbose) {
  if (!LineTes::is_valid(verbose))
    return false;
  if (!is_a_T_tessellation(verbose))
    return false;
  // check non-blocking segments
  for (Seg_sublist_iterator s=non_blocking_segments_begin();
       s!=non_blocking_segments_end();s++) {
    if ((*s)->number_of_edges()!=1) {
      if (verbose) {
	std::clog << "TTessel::is_valid: segment " << (*s)->pointSource();
	std::clog << " " << (*s)->pointTarget() << " considered as ";
	std::clog << "non-blocking has " << (*s)->number_of_edges();
	std::clog << " edge(s)" << std::endl;
      }
      return false;
    }
  }
  // check blocking segments
  for (Seg_sublist_iterator s=blocking_segments_begin();
       s!=blocking_segments_end();s++) {
    if ((*s)->number_of_edges()<2) {
      if (verbose) {
	std::clog << "TTessel::is_valid: segment " << (*s)->pointSource();
	std::clog << " " << (*s)->pointTarget() << " considered as ";
	std::clog << "blocking has " << (*s)->number_of_edges();
	std::clog << " edge(s)" << std::endl;
      }
      return false;
    }
  }
  return true;
}
/** \brief Propose a split applicable to the T-tessellation*/
TTessel::Split TTessel::propose_split() {
  Halfedge_handle e = alt_length_weighted_random_halfedge();
  return Split(e,rnd->get_double(0,1),acos(1-2*rnd->get_double(0,1)));
}
/** \brief Propose a merge applicable to the T-tessellation*/
TTessel::Merge TTessel::propose_merge() {
  Size n = (Size) rnd->get_int(0,non_blocking_segments.size());
  Seg_handle s = non_blocking_segments[n];
  return Merge(s->get_halfedge_handle());
}
/** \brief Propose a flip applicable to the T-tessellation*/
TTessel::Flip TTessel::propose_flip() {
  Size n = (Size) rnd->get_int(0,blocking_segments.size());
  Seg_handle s = blocking_segments[n];
  Halfedge_handle e;
  int side = rnd->get_int(0,2); // 0 or 1
  if (side==0) {
    e = s->halfedges_start()->twin();
  } else {
    e = s->halfedges_end();
  }
  return Flip(e);
}

/** \brief Return a random sample of splits
 * \param n : sample size
 * \return a list of Split objects
 */
TTessel::Split_list TTessel::split_sample(Size n) {
  // Split_list sample(n);
  // for(Size i =0; i < n; i++) {
  //   sample[i] = propose_split();
  // }
  // return sample;
  std::vector<Halfedge_handle> e = alt_length_weighted_random_halfedge(n);
  Split_list splits(n);
  for (unsigned int i=0;i!=e.size();i++) {
    splits[i] = Split(e[i],rnd->get_double(0,1),acos(1-2*rnd->get_double(0,1)));
  }
  return splits;
}

/** \brief Return a Poisson sample of splits
 *
 * The Poissonnian process of splits has an intensity measure
 * proportional to the uniform measure on the space of splits.
 * \param intensity  : intensity of the Poissonnian process as a double.
 * \return a list of Split objects
 */
TTessel::Split_list TTessel::poisson_splits(double intensity) {
  Size split_number;
  double poisson_par = intensity*get_total_internal_length()/CGAL_PI;
  // Knuth algorithm for Poisson generation
  double L = exp(-poisson_par);
  int k = 0;
  double p = 1.0;
  do {
    k++;
    p *= rnd->get_double(0,1);
  } while (p>L);
  split_number = k-1;
  Split_list res = split_sample(split_number);
  return res;
}
/** \brief Return all possible merges
 * \return a list of Merge objects
 */
TTessel::Merge_list TTessel::all_merges() {
  Merge_list merges(non_blocking_segments.size());
  for(Size i = 0; i < merges.size(); i++) {
    Seg_handle s = non_blocking_segments[i];
    merges[i] = Merge(s->get_halfedge_handle());
  }
  return merges;
}

/** \brief Return all possible flips
 * \return a list of Flip objects
 */
TTessel::Flip_list TTessel::all_flips() {
  Size n = blocking_segments.size();
  Flip_list flips(2*n);
  for(Size i = 0; i < n; i++) {
    Seg_handle s = blocking_segments[i];
    Halfedge_handle e;
    e = s->halfedges_start()->twin();
    flips[i] = Flip(e);
    e = s->halfedges_end();
    flips[i+n] = Flip(e);
  }
  return flips;
}
/** \brief Return the number of internal blocking segments.
 */
TTessel::Size TTessel::number_of_blocking_segments() {
  return blocking_segments.size();
}
/** \brief Return the number of internal non-blocking segments.
 */
TTessel::Size TTessel::number_of_non_blocking_segments() {
  return non_blocking_segments.size();
}
/** \brief Return an iterator pointing to the first internal 
 * blocking segment.
 */
TTessel::Seg_sublist_iterator TTessel::blocking_segments_begin() {
  return blocking_segments.begin();
}
/** \brief Return an iterator pointing to the past-the-end internal
 * blocking segment.
 */
TTessel::Seg_sublist_iterator TTessel::blocking_segments_end() {
  return blocking_segments.end();
}
/** \brief Return an iterator pointing to the first internal 
 * non-blocking segment.
 */
TTessel::Seg_sublist_iterator TTessel::non_blocking_segments_begin() {
  return non_blocking_segments.begin();
}
/** \brief Return an iterator pointing to the past-the-end internal
 * non-blocking segment.
 */
TTessel::Seg_sublist_iterator TTessel::non_blocking_segments_end() {
  return non_blocking_segments.end();
}

/** \brief Return all faces of the tessellation as polygons
 */
Polygons TTessel::all_faces() {
  Polygons res;
  for (Face_iterator f = faces_begin(); f!= faces_end(); f++) {
    if (!f->is_unbounded()) {
      Polygon p = face2poly(f);
      res.push_back(p);
    }
  }
  return res;
}
// std::vector<PolygonCoordinates> TTessel::all_faces() {
//   std::vector<PolygonCoordinates> res;
//   for (Face_iterator f = faces_begin(); f!= faces_end(); f++) {
//     std::vector<double> x;
//     std::vector<double> y;
//     PolygonCoordinates c;
//     Polygon p = face2poly(f);
//     for(unsigned int j=0; j!=p.size(); j++) {
//       x.push_back(CGAL::to_double(p[j][0]));
//       y.push_back(CGAL::to_double(p[j][1]));
//     }
//     c["x"] = x;
//     c["y"] = y;
//     res.push_back(c);
//   }
//   return res;
// }

/** \brief Print on standard output all blocking segments
 *
 * Segments are printed using the LineTes::Seg::print method.
 * Coordinates of all vertices along each segment are printed. */
void TTessel::print_blocking_segments() {
  for (Seg_sublist_iterator si=blocking_segments_begin();si!=blocking_segments_end();si++) {
    (*si)->print();
  }
}
/** \brief Print on standard output all non-blocking segments
 *
 * Segments are printed using the LineTes::Seg::print method.
 * Coordinates of all vertices along each segment are printed. */
void TTessel::print_non_blocking_segments() {
  for (Seg_sublist_iterator si=non_blocking_segments_begin();si!=non_blocking_segments_end();si++) {
    (*si)->print();
  }
}


TTessel::Halfedge_handle  TTessel::length_weighted_random_halfedge() {
  /* Draw a halfedge bounding an internal face with probability
     proportional to its length. Rejection algorithm. */

  Halfedge_handle e;
  bool found = false;
  double til = total_internal_length();
  while (!found) {
    e = random_halfedge();
    double rd = rnd->get_double(0,til);
    if (rd<e->get_length())
      found = true;
  }
  return e;
}

TTessel::Halfedge_handle TTessel::alt_length_weighted_random_halfedge() {
  /* Draw a halfedge bounding an internal face with probability
     proportional to its length. Point selection algorithm which seems
     faster than the rejection algorithm. */
  double til = total_internal_length();
  double cumlen = 0.0;
  double select = rnd->get_double(0,til);
  Halfedge_handle res = halfedges_end();
  for (Halfedge_iterator e=halfedges_begin();e!=halfedges_end();e++) {
    if (e->face()->is_unbounded())
      continue;
    cumlen += e->get_length();
    if (cumlen>select) {
      res = e;
      break;
    }
  }
  return res;
}

std::vector<TTessel::Halfedge_handle> TTessel::alt_length_weighted_random_halfedge(unsigned int n) {
  double til = total_internal_length();
  double cumlen = 0.0;
  std::vector<double> select(n);
  for (unsigned int i=0;i!=select.size();i++) {
    select[i]= rnd->get_double(0,til);
  }
  std::sort(select.begin(),select.end(),std::greater<double>());
  std::vector<Halfedge_handle> res;
  for (Halfedge_iterator e=halfedges_begin();e!=halfedges_end();e++) {
    if (e->face()->is_unbounded())
      continue;
    cumlen += e->get_length();
    while (cumlen>select.back()) {
      res.push_back(e);
      select.pop_back();
      if (select.size()==0) 
	break;
    }
    if (select.size()==0) 
      break;
  }
  return res;
}

TTessel::Halfedge_handle TTessel::random_halfedge() {
  /* Draw at random a halfedge bounding an internal face. */
  
  bool found = false;
  Size n = number_of_halfedges();
  Halfedge_iterator e;
  while (!found) {
    e = halfedges_begin();
    Size ri = (Size) rnd->get_int(0,n);
    for (Size i=0;i!=ri;i++)
      e++;
    if (!e->face()->is_unbounded())
      found = true;
  }
  return e; 
}
/** \brief Print the tessellation in a format supported by califlopp 
 * \param out : ouput stream
 */
void TTessel::printRCALI(std::ostream& out){
  /* Print to  output the tessellation in a format supported
   by califlopp */
  out <<  number_of_faces() - number_of_unbounded_faces() << std::endl;
  int faces_counter=0;
  for (Face_iterator f = faces_begin();f!=faces_end();f++) {
    if (! f->is_unbounded())
      {
	faces_counter +=1;
	Ccb_halfedge_circulator e=f->outer_ccb();
	Halfedge_handle e_first=e;
	int vertex_counter=0;
	
	do {
              if (e->get_next_hf()==NULL_HALFEDGE_HANDLE || (e->get_next_hf()!=e->next())) 
		vertex_counter+=1 ;
	      e++;
	}
	while (e!=e_first);
	out <<  faces_counter << " " << faces_counter << " " << vertex_counter << std::endl;
	
	
	e=f->outer_ccb();
	e_first=e;
	out << e->target()->point()[0];
	e++;
	while (e!=e_first){
	  if (e->get_next_hf()==NULL_HALFEDGE_HANDLE || (e->get_next_hf()!=e->next())) 
	    out << " "  <<  e->target()->point()[0] ;
	  e++;
	}
        out << std::endl;
	
	
        e=f->outer_ccb();
        e_first=e;
	
        out << e->target()->point()[1];
        e++;
        while (e!=e_first){
	  if (e->get_next_hf()==NULL_HALFEDGE_HANDLE || (e->get_next_hf()!=e->next())) 
	    out << " "  <<  e->target()->point()[1] ;
	  e++;
	}
       out << std::endl;
       
      } 
  }
}

/******************************************************************************/
/*                     METHODS FOR THE CLASS TTessel::Modification            */
/******************************************************************************/
/** \brief Return the list of tessellation elements that are modified*/
ModList TTessel::Modification::modified_elements(){
  ModList m;
return m;}



/******************************************************************************/
/*                     METHODS FOR THE CLASS TTessel::Split                   */
/******************************************************************************/

/** \brief Default constructor that creates an empty Split object
 */
TTessel::Split::Split() {}
/** \brief Constructor that creates an effective Split object
 * \param e : halfedge handle bounding the cell to be split where the 
 *            splitting segment starts
 * \param x : relative position along the halfedge defining where the
 *            splitting segments starts. When s is closed to 0, the splitting 
 *            segment starts near the source of the halfedge. When s is closed 
 *            to 1, the splitting segment is near the target of the halfedge
 * \param angle : angle between the halfedge and the splitting segment. Must
 *                be between 0 and pi
 * \return A Split object
 */
TTessel::Split::Split(Halfedge_handle e, double x, double angle) : 
  Modification(), e1(e) {

  // Computation of the splitting line

  Vector e1_vector(e1->source()->point(),e1->target()->point());
  double e1_angle = atan2(CGAL::to_double(e1_vector[1]),
			  CGAL::to_double(e1_vector[0]));
  double l_angle = e1_angle+angle;
  double a = sin(l_angle);
  double b = -cos(l_angle);
  double u,v,w;
  u = a*CGAL::to_double(e1->source()->point()[0])+
    b*CGAL::to_double(e1->source()->point()[1]);
  v = a*CGAL::to_double(e1->target()->point()[0])+
    b*CGAL::to_double(e1->target()->point()[1]);
  /*if (u>v) { // exchange u and v
    w = u;
    u = v;
    v = w;
    }*/
  w = u+x*(v-u);
  NT a2(a), b2(b), c2(-w); 
  Line l(a2,b2,c2);
  // Computation of p1
  CGAL::Object inter;
  inter = CGAL::intersection(Segment(e1->source()->point(),
				     e1->target()->point()),l);
  if (!CGAL::assign(p1,inter)) {
    std::cerr << "Split constructor: l does not hit e1" << std::endl;
    std::cerr << "e1: " << e1->source()->point();
    std::cerr << " " << e1->target()->point() << std::endl;
    std::cerr << "l: " << l << std::endl;
  }

  // Computation of e2 and p2

  TTessel::Halfedge_handle ee = e1->next();
  Point2 pp;
  NT minDist = -1;
  Line l_dir;
  if (l.has_on_negative_side(e1->target()->point())) {
    l_dir = l;
  } else {
    l_dir = l.opposite();
  }
  Rayon r(p1,l_dir);

  while (ee!=e1) {
    inter = CGAL::intersection(Segment(ee->source()->point(),
				       ee->target()->point()),r);
    if (CGAL::assign(pp,inter)) {
      NT distp1p2 = CGAL::squared_distance(p1,pp);
      if (minDist>distp1p2 || minDist==-1) {
      	minDist = distp1p2;
      	e2 = ee;
      	p2 = pp;
      }
    }
    ee = ee->next(); 
  }
  if (e2==e1)
    std::cerr << "Split constructor : intersection of l with e2 not found"
	      << std::endl;
}
/** \brief Compute elements of a T-tessellation that are modified by a split
 */ 
ModList TTessel::Split::modified_elements() {

  ModList modifs;

  // No suppressed vertices

  // Added vertices

  modifs.add_vertices.push_back(get_p1());
  modifs.add_vertices.push_back(get_p2());

  // Suppressed edges 

  Halfedge_handle e = get_e1();
  modifs.del_edges.push_back(Segment(e->source()->point(),
				     e->target()->point()));
  e = get_e2();
  modifs.del_edges.push_back(Segment(e->source()->point(),
				     e->target()->point()));
  // Added edges

  modifs.add_edges.push_back(Segment(get_e1()->source()->point(),
				     get_p1()));
  modifs.add_edges.push_back(Segment(get_p1(),
				     get_e1()->target()->point()));
  modifs.add_edges.push_back(Segment(get_e2()->source()->point(),
				     get_p2()));
  modifs.add_edges.push_back(Segment(get_p2(),
				     get_e2()->target()->point()));
  modifs.add_edges.push_back(Segment(get_p1(),get_p2()));

  // Suppressed face

  Polygon poly;
  TTessel::Ccb_halfedge_circulator e_circ=get_e1()->ccb();
  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=get_e1());
  modifs.del_faces.push_back(poly);

  // Add faces

  poly.erase(poly.vertices_begin(),poly.vertices_end());
  poly.push_back(get_p1());
  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=get_e2());
  poly.push_back(get_p2());
  modifs.add_faces.push_back(poly);
  // Another one
  poly.erase(poly.vertices_begin(),poly.vertices_end());
  poly.push_back(get_p2());
  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=get_e1());
  poly.push_back(get_p1());
  modifs.add_faces.push_back(poly);

  // Segment s1 is suppressed
  std::vector<Point2> seg = get_e1()->segment()->list_of_points();
  modifs.del_segs.push_back(seg);
  // Segment s1 with one extra point p1 is added ???
  if (get_e1()->get_dir())
    seg.insert(std::find(seg.begin(),seg.end(),
			 get_e1()->target()->point()),get_p1());
  else
    seg.insert(std::find(seg.begin(),seg.end(),
			 get_e1()->source()->point()),get_p1());
  modifs.add_segs.push_back(seg);

  // Idem for s2
  seg.clear();
  seg = get_e2()->segment()->list_of_points();
  modifs.del_segs.push_back(seg);
  // Segment s2 with one extra point p2 is added
  if (get_e2()->get_dir())
    seg.insert(std::find(seg.begin(),seg.end(),
			 get_e2()->target()->point()),get_p2());
  else
    seg.insert(std::find(seg.begin(),seg.end(),
			 get_e2()->source()->point()),get_p2());
  modifs.add_segs.push_back(seg);

  // New segment
  seg.clear();
  seg.push_back(get_p1()); seg.push_back(get_p2());
  modifs.add_segs.push_back(seg);

  return modifs;
}
/** \brief Check a split is valid
 *
 * The following tests are performed
 * - Both split edges bound the same face.
 * - The split edges should not be aligned.
 * - The new vertices should not be at the end of the split edges.
 *
 * If the split is not valid, messages are sent to the standard
 * output for logging.
*/
bool TTessel::Split::is_valid() {
  // Check that e1 and e2 bound the same face
  Halfedge_handle e = e1->next();
  while (e!=e1)
    if (e==e2)
      break;
    else
      e = e->next();
  if (e==e1) {
    std::clog << "Split::is_valid : e1 and e2 not on the same face"
	      << std::endl;
    return false;
  }
  // Check that e1 and e2 are not aligned
  Line e1_l(e1->source()->point(),e1->target()->point());
  Line e2_l(e2->source()->point(),e2->target()->point());
  if (e1_l==e2_l || e1_l==e2_l.opposite()) {
    std::clog << "Split::is_valid: e1 and e2 are aligned" << std::endl;
    std::clog << "e1: " << e1->source()->point() << " ";
    std::clog << e1->target()->point() << std::endl;
    std::clog << "e2: " << e2->source()->point() << " ";
    std::clog << e2->target()->point() << std::endl;
    return false;
  }
  // Check that p1 is on e1 but not at its ends
  Segment s(e1->source()->point(),e1->target()->point());
  if (!s.has_on(p1)) {
    std::clog << "Split::is_valid : p1 is not on e1" << std::endl;
    return false;
  } else if (p1==e1->source()->point()||p1==e1->source()->point()) {
    std::clog << "Split::is_valid : p1 is at one end of e1" << std::endl;
    return false;
  }
  // Check that p2 is on e2 but not at its ends
  s = Segment(e2->source()->point(),e2->target()->point());
  if (!s.has_on(p2)) {
    std::clog << "Split::is_valid : p2 is not on e2" << std::endl;
    return false;
  } else if (p2==e2->source()->point()||p2==e2->source()->point()) {
    std::clog << "Split::is_valid : p2 is at one end of e2" << std::endl;
    return false;
  }
  return true;
}


/******************************************************************************/
/*                     METHODS FOR THE CLASS TTessel::Merge                   */
/******************************************************************************/

/** \brief Default constructor that generates an empty Merge object.
 */
TTessel::Merge::Merge() {}
/** \brief Constructor that generates an effective Merge object
 * \param he : halfedge handle defining the edge to be removed
 */
TTessel::Merge::Merge(Halfedge_handle he) : e(he) {}
/** \brief Compute the elements of a T-tessellation that are modified by a merge
 */
ModList TTessel::Merge::modified_elements() {

  ModList modifs;

  // Suppressed vertices
  
  modifs.del_vertices.push_back(p1());
  modifs.del_vertices.push_back(p2());
  

  // No added vertices
  
  // Suppressed edges

 modifs.del_edges.push_back(Segment(get_e()->source()->point(),get_e()->target()->point()));
 modifs.del_edges.push_back(Segment(p2(),e2()->target()->point()));
 modifs.del_edges.push_back(Segment(e2()->get_prev_hf()->source()->point(),p2()));
 modifs.del_edges.push_back(Segment(p1(),e1()->target()->point()));
 modifs.del_edges.push_back(Segment(e1()->get_prev_hf()->source()->point(),p1()));

 // Added edges 

 modifs.add_edges.push_back(Segment(e2()->get_prev_hf()->source()->point(),e2()->target()->point()));
 modifs.add_edges.push_back(Segment(e1()->get_prev_hf()->source()->point(),e1()->target()->point()));


  // Suppressed faces

  Polygon poly;
  TTessel::Ccb_halfedge_circulator e_circ=e1()->ccb();

  // Face 1  
  poly.push_back(p1());
  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=get_e()->twin());
  modifs.del_faces.push_back(poly);
  
  // Face 2
  
  poly.erase(poly.vertices_begin(),poly.vertices_end()); 
  e_circ=e2()->ccb();
  poly.push_back(p2());

  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=get_e());

  modifs.del_faces.push_back(poly);
 
   // Added face
   
  poly.erase(poly.vertices_begin(),poly.vertices_end()); 
  e_circ=e1()->ccb();
  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=e2()->get_prev_hf());
  e_circ=e2()->ccb();
  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=e1()->get_prev_hf());

  modifs.add_faces.push_back(poly);
  

  // Suppressed and added segments

  // Seg 1 
   
  std::vector<Point2> seg = e1()->segment()->list_of_points();
  Halfedge_handle e;

  modifs.del_segs.push_back(seg);

  // Segment without p1 is added

  seg.erase(std::find(seg.begin(),seg.end(),p1()));
  modifs.add_segs.push_back(seg);
 
  // Seg 2
  
  seg.clear();
  seg = e2()->segment()->list_of_points();
  modifs.del_segs.push_back(seg);

   // Segment without p2 is added
  
  seg.erase(std::find(seg.begin(),seg.end(),p2()));
  modifs.add_segs.push_back(seg);
  
  // Seg 
  
  seg.clear();
  seg.push_back(p1());
  seg.push_back(p2());

  modifs.del_segs.push_back(seg);

  // Added segments

  return modifs;
}
/** \brief Check a merge is valid
 *
 * The merge is considered as valid if the edge to be removed 
 * is the only one on its segment.*/
bool TTessel::Merge::is_valid() {
  if (e->get_prev_hf()!=NULL_HALFEDGE_HANDLE || e->get_next_hf()!=NULL_HALFEDGE_HANDLE) {
    std::clog << "Merge::is_valid: the segment supporting e has more than";
    std::clog << " edge" << std::endl;
    return false;
  }
  return true;
}

 
/******************************************************************************/
/*                     METHODS FOR THE CLASS TTessel::Flip                    */
/******************************************************************************/

/** \brief Default constructor that returns an empty Flip object
 */
TTessel::Flip::Flip(){}
/** \brief Constructor that generates an effective Flip object
 * \param he : halfedge handle defining the edge (lying at the end of a blocking
 *             segment) that is removed by the flip.
 */
TTessel::Flip::Flip(Halfedge_handle he) : e1(he) {

  // Find the edge which must be extended
  Seg_handle s1 = e1->segment(); // segment to be shortened
  Halfedge_around_vertex_circulator e3 = e1->source()->incident_halfedges();
  while (e3->segment()==s1)
    e3++;

  // Find a halfedge along the face to be split
  Point2 pt1(e3->source()->point()); 
  Point2 pt2(e3->target()->point());
  Point2 pt3(e3->next()->target()->point());

  /* hf_it_deb : halfedge starting at the end of e1 and bounding the face
     to be split. */
  Halfedge_handle hf_it_deb;
  if (left_turn(pt1,pt2,pt3))  
    hf_it_deb = e3->next()->twin()->next();
  if (right_turn(pt1,pt2,pt3)) 
    hf_it_deb = e3->next();

  /*if (is_on_boundary(hf_it_deb)!=0) {
    std::cerr << "Flip constructor: e1 is incident to the boundary"
    << std::endl;
    return;
  }*/

  // Point at the end of the segment to be extended
  Point2 seg_deb(e1->source()->point());
  // Halfline
  Rayon r(seg_deb,Vector(e3->source()->point(),e3->target()->point())); 

  // Search of the halfedge e2 hit by the ray r along the face to be split
  NT minDist = -1;
  Point2 pp;
  Halfedge_handle ee = hf_it_deb->next();
  while (ee != hf_it_deb->prev()) {
    CGAL::Object inter;
    inter = CGAL::intersection(r,Segment(ee->source()->point(),
					 ee->target()->point()));
    if (CGAL::assign(pp,inter)) {
      NT distsegdebpp = CGAL::squared_distance(seg_deb,pp);
      if (minDist>distsegdebpp || minDist==-1) {
      	minDist = distsegdebpp;
      	e2 = ee;
	p2 = pp;
      }
    }
    ee = ee->next(); 
  }    
  if (minDist<0)
    std::cerr << "Flip constructor: e2 and p2 not found" << std::endl;
}
/** \brief Return the list of tessellation elements modified by a flip
 *
 * Tessellation elements are vertices, edges, faces and segments. Possible
 * modifications are removal and creation, see ModList.*/
ModList TTessel::Flip::modified_elements() {

  ModList modifs;

  // Suppressed vertices
  modifs.del_vertices.push_back(get_e1()->target()->point());

  // Added vertices
  modifs.add_vertices.push_back(get_p2());

 
  // Suppressed edges :
  
  Segment s_split,s1_split,s2_split,s_merge,s1_merge,s2_merge;
  Point2 p_merge;
  
  s_merge =  Segment(get_e1()->source()->point(),get_e1()->target()->point());
  s1_merge = Segment(get_e1()->next()->source()->point(),get_e1()->next()->target()->point());
  s2_merge = Segment(get_e1()->next()->get_prev_hf()->target()->point(),get_e1()->next()->get_prev_hf()->source()->point());
  p_merge = s1_merge[0];
  
  s_split = Segment(get_e1()->source()->point(),get_p2());
  s1_split = Segment(get_p2(),get_e2()->target()->point());
  s2_split = Segment(get_p2(),get_e2()->source()->point());
  
  modifs.del_edges.push_back(s_merge);
  modifs.del_edges.push_back(s1_merge); 
  modifs.del_edges.push_back(s2_merge);
  
  if (s1_split[1]!=p_merge && s2_split[1]!=p_merge) {
  	modifs.del_edges.push_back(Segment(s1_split[1],s2_split[1]));
  }
  
  
  // Added edges:
  
  modifs.add_edges.push_back(s_split);
  
  
  if (s1_split[1]!=p_merge && s2_split[1]!=p_merge)  {
  	modifs.add_edges.push_back(s1_split);
  	modifs.add_edges.push_back(s2_split);
  	modifs.add_edges.push_back(Segment(s1_merge[1],s2_merge[1]));
  }
  
  if (s1_split[1]==p_merge){
    modifs.add_edges.push_back(s2_split);
    modifs.add_edges.push_back(Segment(s1_split[0],s1_merge[1]));
  }
  
  if (s2_split[1]==p_merge){
    modifs.add_edges.push_back(s1_split);
    modifs.add_edges.push_back(Segment(s2_split[0],s2_merge[1]));
  }
  
  
  
    
      // 2 suppressed faces

  Polygon poly;
  TTessel::Ccb_halfedge_circulator e_circ = get_e1()->ccb();
  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=get_e1());
  modifs.del_faces.push_back(poly);
  poly.erase(poly.vertices_begin(),poly.vertices_end());
  e_circ=get_e1()->twin()->ccb();
  do {
    poly.push_back(e_circ->target()->point());
  } while (++e_circ!=get_e1()->twin());
  modifs.del_faces.push_back(poly);

  // 2 added faces

  // Extended face
   
  poly.erase(poly.vertices_begin(),poly.vertices_end());
  e_circ = get_e1()->ccb();
  
  do { 
  	poly.push_back(e_circ->target()->point());
    e_circ++;
  } while (e_circ!=get_e2() && e_circ!=get_e1());
  poly.erase(poly.vertices_begin());
  
  bool flag = e_circ==get_e2();
  poly.push_back(get_p2());
  if (flag)
    e_circ = get_e1()->twin()->ccb();
  else
    e_circ = get_e2()->ccb();
    
  while (e_circ->target()->point()!=get_e1()->target()->point()) {
    poly.push_back(e_circ->target()->point());
    e_circ++;
  }
  modifs.add_faces.push_back(poly);
  
  // Shortened face
  
  
  poly.erase(poly.vertices_begin(),poly.vertices_end());
  poly.push_back(get_p2());
  if (flag)
    e_circ = get_e2()->ccb();
  else
    e_circ = get_e1()->twin()->ccb();
  do {
    poly.push_back(e_circ->target()->point());
    e_circ++;
  } while (e_circ!=get_e1() && e_circ!=get_e2()); 
  
  
  modifs.add_faces.push_back(poly);
  
  // Segments

  std::vector<Point2> seg;

  // Segment to be shortened
  seg = get_e1()->segment()->list_of_points();
  modifs.del_segs.push_back(seg);
  if (get_e1()->get_dir())
    seg.pop_back();
  else
    seg.erase(seg.begin());
  modifs.add_segs.push_back(seg);

  // Segment to be extended
  seg.clear();
  Seg_handle s1 = get_e1()->segment(); // segment to be shortened
  Halfedge_around_vertex_circulator e3 = 
    get_e1()->source()->incident_halfedges();
  while (e3->segment()==s1)
    e3++;
  seg = e3->segment()->list_of_points();
  modifs.del_segs.push_back(seg);
  if (e3->get_dir())
    seg.push_back(get_p2());
  else
    seg.insert(seg.begin(),get_p2());
  modifs.add_segs.push_back(seg);

  // Segment to be split
  seg.clear();
  seg = get_e2()->segment()->list_of_points();
  modifs.del_segs.push_back(seg); // Segment before split
  // Compute the segment after split
  if (get_e2()->get_dir())
      seg.insert(std::find(seg.begin(),seg.end(),get_e2()->target()->point()),
		 get_p2());
    else
      seg.insert(std::find(seg.begin(),seg.end(),get_e2()->source()->point()),
		 get_p2());
  if (get_e1()->next()->segment()==get_e2()->segment()) {
    // Segments to be split and merged are the same
    // Remove the junction at the merged edges and exit the function
    seg.erase(std::find(seg.begin(),seg.end(),get_e1()->target()->point()));
    modifs.add_segs.push_back(seg);
    return modifs;
  }
  modifs.add_segs.push_back(seg);
  
  // Segment to be merged
  seg.clear();
  seg = get_e1()->next()->segment()->list_of_points();
  modifs.del_segs.push_back(seg); // Segment before merging
  seg.erase(std::find(seg.begin(),seg.end(),get_e1()->target()->point()));
  modifs.add_segs.push_back(seg); // Merged segment

  return modifs;
}

/******************************************************************************/
/*                   METHODS FOR THE CLASS CatItems                           */
/******************************************************************************/

/** \brief Test whether a CatItems object is empty
 */
template <typename T>
bool CatItems<T>::IsEmpty() {
  return (vertices.size()+edges.size()+faces.size()+segs.size())==0;
}
  
/******************************************************************************/
/*                   METHODS FOR THE CLASS CatVector                         */
/******************************************************************************/

/** \brief Fill a CatVector object with the elements of a FVector object
 *
 * \param v : the CatVector to be filled
 * \param x : the input object
 *
 * The lengths of arguments must fit.
 */
void fill(CatVector& v,FVector& x) {
  int j = 0;
  for (unsigned int i=0;i<v.vertices.size();i++,j++)
    v.vertices[i] = x[j];
  for (unsigned int i=0;i<v.edges.size();i++,j++)
    v.edges[i] = x[j];
  for (unsigned int i=0;i<v.faces.size();i++,j++)
    v.faces[i] = x[j];
  for (unsigned int i=0;i<v.segs.size();i++)
    v.segs[i] = x[j];
}
/** \brief Generate a vector of doubles from a CatVector object
 */ 
std::vector<double> asVectorOfDoubles(const CatVector& v) {
  std::vector<double> res;
  for (unsigned int i=0;i<v.vertices.size();i++)
    res.push_back(v.vertices[i]);
  for (unsigned int i=0;i<v.edges.size();i++)
    res.push_back(v.edges[i]);
  for (unsigned int i=0;i<v.faces.size();i++)
    res.push_back(v.faces[i]);
  for (unsigned int i=0;i<v.segs.size();i++)
    res.push_back(v.segs[i]);
  return res;
}
/** \brief Generate a FVector object from CatVector object
 */
FVector asFVector(const CatVector& v) {
  std::vector<double> buf = asVectorOfDoubles(v);
  FVector res(buf.begin(),buf.end());
  return res;
}
/** \brief Print a CatVector to an output stream
 */
std::ostream& operator<<(std::ostream &os, const CatVector &p) {
  if (p.vertices.size()>0) {
    for(unsigned int i=0;i!=p.vertices.size();i++) {
      os << p.vertices[i] << " ";
    }
    os << " (vertices) ";
  }
  if (p.edges.size()>0) {
    for(unsigned int i=0;i!=p.edges.size();i++) {
      os << p.edges[i] << " ";
    }
    os << " (edges) ";
  }
  if (p.faces.size()>0) {
    for(unsigned int i=0;i!=p.faces.size();i++) {
      os << p.faces[i] << " ";
    }
    os << " (faces) ";
  }
  if (p.segs.size()>0) {
    for(unsigned int i=0;i!=p.segs.size();i++) {
      os << p.segs[i] << " ";
    }
    os << " (segs) ";
  }
  return os;
}

/** \brief Compute the sum of components of a CatVector object
 */
double sum(const CatVector &p) {
  typedef std::vector<double> V;
  double res = 0.0;
  for(V::const_iterator v=p.vertices.begin(); v!=p.vertices.end();v++) {
    res += *v;
  }
  for(V::const_iterator v=p.edges.begin(); v!=p.edges.end();v++) {
    res += *v;
  }
  for(V::const_iterator v=p.faces.begin(); v!=p.faces.end();v++) {
    res += *v;
  }
  for(V::const_iterator v=p.segs.begin(); v!=p.segs.end();v++) {
    res += *v;
  }
  return res;
}
/** \brief Generate a FMatrix object from a CatMatrix object
 */
FMatrix asFMatrix(const CatMatrix& m) {
  std::vector<FVector> buf;
  for (unsigned int i=0;i<m.vertices.size();i++) 
    buf.push_back(asFVector(m.vertices[i]));
  for (unsigned int i=0;i<m.edges.size();i++) 
    buf.push_back(asFVector(m.edges[i]));
  for (unsigned int i=0;i<m.faces.size();i++) 
    buf.push_back(asFVector(m.faces[i]));
  for (unsigned int i=0;i<m.segs.size();i++) 
    buf.push_back(asFVector(m.segs[i]));
  FMatrix res(buf);
  return res;
}
/** \brief Outer product of two CatVector objects
 */
CatMatrix outer(const CatVector &v1, const CatVector &v2) {
  CatMatrix m;
  for(unsigned int i=0;i!=v1.vertices.size();i++)
    m.vertices.push_back(v1.vertices[i]*v2);
  for(unsigned int i=0;i!=v1.edges.size();i++)
    m.edges.push_back(v1.edges[i]*v2);
  for(unsigned int i=0;i!=v1.faces.size();i++)
    m.faces.push_back(v1.faces[i]*v2);
  for(unsigned int i=0;i!=v1.segs.size();i++)
    m.segs.push_back(v1.segs[i]*v2);
  return m;
}


/******************************************************************************/
/*                     METHODS FOR THE CLASS Energy                           */
/******************************************************************************/

Energy::Energy(): value(0) {}

/** \brief Define the T-tessellation whose energy is to be monitored
 *
 * Current energy is computed from scratch.
 */
void Energy::set_ttessel(TTessel *t) {
  ttes = t;
  value=0;

  for (Size i=0;i!=theta.vertices.size();i++) {
    for (TTessel::Vertex_iterator v=ttes->vertices_begin();v!=ttes->vertices_end();
	 v++)
      value += theta.vertices[i]*features.vertices[i](v->point(),ttes);
  }
    
  for (Size i=0;i!=theta.edges.size();i++) {
    for (TTessel::Edge_iterator e=ttes->edges_begin();e!=ttes->edges_end();
	 e++)
      value += theta.edges[i]*features.edges[i](Segment(e->source()->point(),
							e->target()->point()),ttes);
  }
      
  for (Size i=0;i!=theta.faces.size();i++) {
    for (TTessel::Face_iterator f=ttes->faces_begin();f!=ttes->faces_end();
	 f++){
      if (!f->is_unbounded()) // only internal faces
	    value += theta.faces[i]*features.faces[i](face2poly(f),ttes);
	 
    }
  }

  for (Size i=0;i!=theta.segs.size();i++) {
    for(TTessel::Seg_list_iterator s=ttes->segments_begin();s!=ttes->segments_end();s++) {
      TTessel::Halfedge_handle e = (*s)->halfedges_start();
      if (ttes->is_on_boundary(e)!=0)
	continue;
      value += theta.segs[i]*features.segs[i]((*s)->list_of_points(),ttes);
    }
  }
}

/** \brief Compute the statistic variation 
 *
 * Compute the statistic variation for a given modification
 * (either a split or merge or flip). The statistic is the vector of the
 * the functions \f$\phi_i\f$'s. For a given modification \f$u\f$, variations
 * \f[
 * \phi_i(uT)-\phi_i(T)
 * \f]
 * are computed. Since functions \f$\phi_i\f$'s are sums over basic 
 * tessellation components (vertices, edges etc...), variations can be computed
 * from the lists of added and removed components under the modification 
 * \f$u\f$.
 *
 * \param modif : either a Split or Merge or Flip object
 * \return component-by-component statistic variation as a Parameter object
 */
CatVector Energy::statistic_variation(TTessel::Modification& modif){ 
  double loc_var=0;
  ModList ml = modif.modified_elements();
  CatVector stat_var;

  stat_var.vertices = std::vector<double>(theta.vertices.size(),0.0);
  stat_var.edges = std::vector<double>(theta.edges.size(),0.0);
  stat_var.faces = std::vector<double>(theta.faces.size(),0.0);
  stat_var.segs = std::vector<double>(theta.segs.size(),0.0);

  // Vertices
  for (unsigned int i=0;i!=theta.vertices.size();i++){
    loc_var=0;
    for (std::vector<Point2>::iterator i_pt=ml.del_vertices.begin(); 
	 i_pt!=ml.del_vertices.end();i_pt++){
      loc_var -=(features.vertices[i])(*i_pt,ttes);
    } 
    for (std::vector<Point2>::iterator i_pt=ml.add_vertices.begin(); 
	 i_pt!=ml.add_vertices.end();i_pt++){
      loc_var +=(features.vertices[i])(*i_pt,ttes);
    }
    stat_var.vertices[i] = loc_var;
  }
  // Edges
  for (unsigned int i=0;i!=theta.edges.size();i++){
    loc_var=0;
    for (std::vector<Segment>::iterator i_ed=ml.del_edges.begin(); 
	 i_ed!=ml.del_edges.end();i_ed++){
      loc_var -=(features.edges[i])(*i_ed,ttes);
    } 
    for (std::vector<Segment>::iterator i_ed=ml.add_edges.begin(); 
	 i_ed!=ml.add_edges.end();i_ed++){
      loc_var +=(features.edges[i])(*i_ed,ttes);
    }
    stat_var.edges[i] = loc_var;
  }

  // Faces
  for (unsigned int i=0;i!=theta.faces.size();i++){
    loc_var=0;
    for (std::vector<Polygon>::iterator i_fc=ml.del_faces.begin(); 
	 i_fc!=ml.del_faces.end();i_fc++){
      loc_var -=(features.faces[i])(*i_fc,ttes);
    } 
    for (std::vector<Polygon>::iterator i_fc=ml.add_faces.begin(); 
	 i_fc!=ml.add_faces.end();i_fc++){
      loc_var += (features.faces[i])(*i_fc,ttes);
    }
    stat_var.faces[i] = loc_var;
  }

  // Segments
  for (unsigned int i=0;i!=theta.segs.size();i++){
    loc_var=0;
    for (std::vector<std::vector<Point2> >::iterator i_sg=ml.del_segs.begin(); 
	 i_sg!=ml.del_segs.end();i_sg++){
      loc_var -=(features.segs[i])(*i_sg,ttes);
    } 
    for (std::vector<std::vector<Point2> >::iterator i_sg=ml.add_segs.begin(); 
	 i_sg!=ml.add_segs.end();i_sg++){
      loc_var +=(features.segs[i])(*i_sg,ttes);
    }
    stat_var.segs[i] = loc_var;
  }

  return stat_var;
}

/** \brief Compute the energy variation 
 *
 * If \f$u\f$ is the considered 
 * modification (either a split or merge or flip),
 * \f[
 * \sum_i \theta_i (\phi_i(uT)-\phi_i(T))
 * \f]
 * is returned.
 *
 * \param modif : either a Split or Merge or Flip object
 * \return the energy variation
 */
double  Energy::variation(TTessel::Modification& modif){ 
  double var=0;
  CatVector stat_var;

  stat_var = statistic_variation(modif);
  for (unsigned int i=0;i!=theta.vertices.size();i++){
    var += stat_var.vertices[i]*theta.vertices[i];
  }
  for (unsigned int i=0;i!=theta.edges.size();i++){
    var += stat_var.edges[i]*theta.edges[i];
  }
  for (unsigned int i=0;i!=theta.faces.size();i++){
    var += stat_var.faces[i]*theta.faces[i];
  }
  for (unsigned int i=0;i!=theta.segs.size();i++){
    var += stat_var.segs[i]*theta.segs[i];
  }
  return var;
}
/** \brief Add a theta component for a vertex feature in the energy formula 
 *
 * /param tv : the numerical value of the theta component to be added */
void  Energy::add_theta_vertices(double tv){
  theta.vertices.push_back(tv);
}
/** \brief Add a theta component for a edge feature in the energy formula 
 *
 * /param te : the numerical value of the theta component to be added */
void  Energy::add_theta_edges(double te){
  theta.edges.push_back(te);
}
/** \brief Add a theta component for a face feature in the energy formula 
 *
 * /param tv : the numerical value of the theta component to be added */
void  Energy::add_theta_faces(double tf){
  theta.faces.push_back(tf);
}
/** \brief Add a theta component for a segment feature in the energy formula 
 *
 * /param ts : the numerical value of the theta component to be added */
void  Energy::add_theta_segs(double ts){
  theta.segs.push_back(ts);
}
/** \brief Remove all theta components associated with vertex features*/
void  Energy::del_theta_vertices(){
  theta.vertices.clear();
}
/** \brief Remove all theta components associated with edge features*/
void  Energy::del_theta_edges(){
  theta.edges.clear();
}
/** \brief Remove all theta components associated with face features*/
void  Energy::del_theta_faces(){
  theta.faces.clear();
}
/** \brief Remove all theta components associated with segment features*/
void  Energy::del_theta_segs(){
  theta.segs.clear();
}
/** \brief Add a vertex feature in the energy formula
 *
 * \param fv : a pointer to a function that takes as input arguments
 *             the vertex location and the T-tessellation. The function
 *             should return the contribution to the energy of the given
 *             vertex.*/
void  Energy::add_features_vertices(double (*fv)(Point2,TTessel*) ){
  features.vertices.push_back(fv);
}
/** \brief Add an edge feature in the energy formula
 *
 * \param fe : a pointer to a function that takes as input arguments
 *             an edge (as a line segment) and the T-tessellation. The function
 *             should return the contribution to the energy of the given
 *             edge.*/
void  Energy::add_features_edges(double (*fe)(Segment,TTessel*) ){
  features.edges.push_back(fe);
}
/** \brief Add a face feature in the energy formula
 *
 * \param ff : a pointer to a function that takes as input arguments
 *             a face (as a polygon) and the T-tessellation. The function
 *             should return the contribution to the energy of the given
 *             face.*/
void  Energy::add_features_faces(double (*ff)(Polygon,TTessel*) ){
  features.faces.push_back(ff);
}
/** \brief Add a segment feature in the energy formula
 *
 * \param fs : a pointer to a function that takes as input arguments
 *             a segment and the T-tessellation. The function
 *             should return the contribution to the energy of the given
 *             segment.*/
void  Energy::add_features_segs(double (*fs)(std::vector<Point2>,TTessel*) ){
  features.segs.push_back(fs);
}
/** \brief Remove all vertex features from the energy formula*/
void  Energy::del_features_vertices(){
  features.vertices.clear();
}
/** \brief Remove all edge features from the energy formula*/
void  Energy::del_features_edges(){
  features.edges.clear();
}
/** \brief Remove all face features from the energy formula*/
void  Energy::del_features_faces(){
  features.faces.clear();
}
/** \brief Remove all segment features from the energy formula*/
void  Energy::del_features_segs(){
  features.segs.clear();
}

/******************************************************************************/
/*                       METHODS FOR THE CLASS SMFChain                       */
/******************************************************************************/

/** \brief Constructor for class SMFChain
 *
 * \param e : pointer to an Energy object defining the distribution to be
 *            simulated.
 * \param prob_split : probability of proposing a split
 * \param prob_merge : probability of proposing a merge
 *
 * The probability of proposing a flip is complementary to the sum of
 * probabilities of proposing a split or a merge. Therefore, the sum of
 * the latter probabilities should be less than 1.
 */
SMFChain::SMFChain(Energy *e, double prob_split, double prob_merge) : 
  engy(e), p_split(prob_split), p_merge(prob_merge) {
  if (p_split+p_merge>1)
  std::cerr << "SMFChain constructor: sum of probabilities > 1" << std::endl;
  p_split_merge = p_split+p_merge;
  p_flip = 1.0-p_split_merge;
}
/** \brief Set the probabilities of proposing a split and a merge
 * \param s : probability to propose a split
 * \param m : probability to propose a merge
 *
 * The sum s+m should be smaller than 1. The probability to propose
 * a flip is 1-s-m.*/  
void SMFChain::set_smf_prob(double s,double m) {
  p_split = s;
  p_merge = m;
  if (p_split+p_merge>1)
    std::cerr << "SMFChain constructor: sum of probabilities > 1" << std::endl;
  p_split_merge = p_split+p_merge;
  p_flip = 1.0-p_split_merge;
}
/** \brief Draw a type of modification (split or merge or flip)*/
SMFChain::ModType SMFChain::propose_modif_type() {
  double r = rnd->get_double(0,1);
  if (r<=p_split)
    return SPLIT;
  if (r<=p_split_merge)
    return MERGE;
  else
    return FLIP;
}
/** \brief Compute the Hastings ratio for a given split
 *
 * \param s : the split for which the Hastings ratio is computed
 * \param e_var : the address of a double where to store the energy
 *                variation
 * \return the value of the Hastings ratio*/
double SMFChain::Hasting_ratio(TTessel::Split& s,double *e_var) {

  double  r;
  int delta_nb=1;
  TTessel *tesl = engy->get_ttessel();
  r = p_merge/p_split;
  // A split may change a non-blocking segment into a blocking one
  TTessel::Seg_handle s_split = s.get_e1()->segment();
  if (s_split->number_of_edges()==1 && !tesl->is_on_boundary(s_split))
    delta_nb--;
  s_split = s.get_e2()->segment();
  if (s_split->number_of_edges()==1 && !tesl->is_on_boundary(s_split))
    delta_nb--;
  r *= (1/CGAL_PI)*tesl->get_total_internal_length()/
    (tesl->number_of_non_blocking_segments()+delta_nb);
  *e_var=engy->variation(s);
  r *= exp(-*e_var);
  return r;
}
/** \brief Compute the Hastings ratio for a given merge
 *
 * \param m : the merge for which the Hastings ratio is computed
 * \param e_var : the address of a double where to store the energy
 *                variation
 * \return the value of the Hastings ratio*/
double SMFChain::Hasting_ratio(TTessel::Merge& m,double *e_var) {
  double r;
  TTessel *tesl = engy->get_ttessel();
  r = p_split/p_merge;
  r *=(CGAL_PI/1)* tesl->number_of_non_blocking_segments()/
    (tesl->get_total_internal_length()-2*m.get_e()->get_length());
  *e_var =engy->variation(m);
  r *= exp(-*e_var);
  return r;
}
/** \brief Compute the Hastings ratio for a given flip
 *
 * \param f : the flip for which the Hastings ratio is computed
 * \param e_var : the address of a double where to store the energy
 *                variation
 * \return the value of the Hastings ratio*/
double SMFChain::Hasting_ratio(TTessel::Flip& f,double *e_var) {
  
  double r;
  TTessel *tesl = engy->get_ttessel();
  
  // Computation of the variation of the number of blocking segments
  // due to the flip
  int delta_blocking = 0;
 
  // Shortened segment
  TTessel::Seg_handle s_short = f.get_e1()->segment();
  if (s_short->number_of_edges()==2)
    delta_blocking--;

  // Extended segment
  /* Find the edge which must be prolonged. ENCORE ET ENCORE ! ECRIRE
     UNE METHODE. */
  TTessel::Halfedge_around_vertex_circulator e3;
  e3 = f.get_e1()->source()->incident_halfedges();
  while (e3->segment()==s_short)
    e3++;
  TTessel::Seg_handle s_exten = e3->segment();
  if (s_exten->number_of_edges()==1)
    delta_blocking++;

  // Incident segments
  TTessel::Seg_handle incident_s_short = f.get_e1()->next()->segment();
  TTessel::Seg_handle incident_s_exten = f.get_e2()->segment();
  if (incident_s_short!=incident_s_exten) {
    if (incident_s_short->number_of_edges()==2 && !tesl->is_on_boundary(incident_s_short))
      delta_blocking--;
    if (incident_s_exten->number_of_edges()==1 && !tesl->is_on_boundary(incident_s_exten))
      delta_blocking++;
  } // if only one incident segment, remains blocking

  Size sb =  tesl->number_of_blocking_segments(); 
 
  if (sb ==(Size) -delta_blocking) {
    r = sb;
    std::cerr << "Hasting_ratio(Flip): predicted number of blocking segments";
    std::cerr << " after flip = 0" << std::endl;
    std::cerr << "   Current number of segments: ";
    std::cerr << tesl->number_of_segments()-tesl->number_of_window_edges() << std::endl;
    std::cerr << "   Current number of non-blocking segments: ";
    std::cerr << tesl->number_of_non_blocking_segments() << std::endl;
    std::cerr << "   Current number of blocking segments: ";
    std::cerr << sb << std::endl;
    //tesl->display();
  } else
    r = 1.0*sb/(sb+delta_blocking);
  *e_var = engy->variation(f);
  r *= exp(-*e_var);
  return r;
}
/** \brief Proceed to update(s) of the SMF chain
 * \param nb_step : number of updates to be done
 * \return some statistics on the proposals and their
 * acceptation*/
ModifCounts SMFChain::step(unsigned long int nb_step) {
  TTessel::Split s;
  TTessel::Merge m;
  TTessel::Flip  f;
  double         delta_engy;
  TTessel        *tesl = engy->get_ttessel();
  ModifCounts    mc={0,0,0,0,0,0};
    
  for(unsigned long int i=0;i!=nb_step;i++) {
    ModType mt = propose_modif_type();
    double x = rnd->get_double(0,1);
    double acs_coef, r;
    switch(mt) {
    case SPLIT:
      mc.proposed_S +=1;
      s = tesl->propose_split();
      r = Hasting_ratio(s,&delta_engy);
      if (r>1 || x<r){
        mc.accepted_S +=1;
	tesl->update(s);
	engy->add_value(delta_engy);
	
      }
      break;
    case MERGE:
      mc.proposed_M +=1;
      if(tesl->number_of_non_blocking_segments()>0) {
	m = tesl->propose_merge();
        r = Hasting_ratio(m,&delta_engy);
	if (r>1 || x<r){
          mc.accepted_M +=1;
	  tesl->update(m);
	  engy->add_value(delta_engy);
	}
      }
      break;
    case FLIP:
      mc.proposed_F +=1;
      if (tesl->number_of_blocking_segments()>0) {
	f = tesl->propose_flip();
	r = Hasting_ratio(f,&delta_engy);
	if (r>1 || x<r) {
          mc.accepted_F +=1;
	  tesl->update(f);
	  engy->add_value(delta_engy);
	}
      }
      break;
    default:
      std::cerr << "INVALID MODIFICATION" <<std::endl;
      break;
    }
  }
  return mc;
}

/******************************************************************************/
/*                     METHODS FOR THE CLASS PseudoLikDiscrete                */
/******************************************************************************/
/** \brief Constructor initializing the model to be considered
 *
 * Some components of the discrete approximation of the pseudolikelihood
 * related to merges and flips are
 * precomputed.*/
PseudoLikDiscrete::PseudoLikDiscrete(Energy *e) : engy(e) {
  typedef TTessel::Merge_list Merges;
  typedef TTessel::Flip_list  Flips;
  TTessel* tes = engy->get_ttessel();
  // Computation of sum_stat_merges
  sum_stat_merges = 0.0*engy->get_theta();
  Merges merges = tes->all_merges();
  for (Merges::iterator mit=merges.begin(); mit!=merges.end();mit++) {
    TTessel::Merge & m = *mit;
    sum_stat_merges -= engy->statistic_variation(m);
  }
  // Computation of stat_flips and sum_stat_flips
  sum_stat_flips = 0.0*engy->get_theta();
  Flips flips = tes->all_flips();
  for(Flips::iterator fit=flips.begin(); fit!=flips.end(); fit++) {
    TTessel::Flip & f = *fit;
    CatVector p = engy->statistic_variation(f);
    stat_flips.push_back(-1.0*p);
    sum_stat_flips -= p;
  }
}
/** \brief Increment the sample of dummy splits
 *
 * \param n : the number of dummy splits to be added.
 *
 * Terms of the approximation of the pseudolikelihood related to splits
 * are updated.
*/
void PseudoLikDiscrete::AddSplits(Size n) {
  typedef TTessel::Split_list Splits;
  TTessel* tes = engy->get_ttessel();
  Splits splits = tes->split_sample(n);
  for(Splits::iterator sit=splits.begin(); sit!=splits.end(); sit++) {
    TTessel::Split & s = *sit;
    CatVector p = engy->statistic_variation(s);
    stat_splits.push_back(-1.0*p);
  }
}

/** \brief Add a given dummy split
 *
 * \note Is this really useful?
 */
void PseudoLikDiscrete::AddGivenSplit(TTessel::Split s) {
  CatVector p = engy->statistic_variation(s);
  stat_splits.push_back(-1.0*p);
}

/** \brief Remove all dummy splits
 */
void PseudoLikDiscrete::ClearSplits() {
  stat_splits.clear();
}

/** \brief Return log-pseudo-likelihood value
 * \param theta : parameter for which the likelihood
 * is computed.
 * \return the value of the discrete approximation of the log-pseudolikelihood.
 */
double PseudoLikDiscrete::GetValue(CatVector theta) {
  typedef std::vector<CatVector> VPar;
  double lpl, buf;
  lpl = -sum(theta*sum_stat_merges);
  buf = 0.0;
  for(VPar::iterator p=stat_splits.begin(); p!=stat_splits.end(); p++) {
    buf += exp(sum(theta*(*p)));
  }
  TTessel* tes = engy->get_ttessel();
  lpl -= tes->get_total_internal_length()/(CGAL_PI*stat_splits.size())*buf;
  lpl -= sum(theta*sum_stat_flips);
  buf = 0.0;
  for(VPar::iterator p=stat_flips.begin(); p!=stat_flips.end(); p++) {
    buf += exp(sum(theta*(*p)));
  }
  lpl -= buf;
  return lpl;
}

/** \brief Return pseudo-likelihood gradient
 * \param theta : parameter where the gradient of the log-pseudolikelihood
 * has to be computed.
 * \return the value of the gradient of the discrete approximation of the
 * log-pseudolikelihood
 */
CatVector PseudoLikDiscrete::GetGradient(CatVector theta) {
  typedef std::vector<CatVector> VPar;
  CatVector g, buf;
  g = -1.0*sum_stat_merges;
  buf = 0.0*g;
  for(VPar::iterator p=stat_splits.begin(); p!=stat_splits.end(); p++) {
    buf += exp(sum(theta*(*p)))*(*p);
  }
  TTessel* tes = engy->get_ttessel();
  g -= tes->get_total_internal_length()/(CGAL_PI*stat_splits.size())*buf;
  g -= sum_stat_flips;
  buf = 0.0*buf;
  for(VPar::iterator p=stat_flips.begin(); p!=stat_flips.end(); p++) {
    buf += exp(sum(theta*(*p)))*(*p);
  }
  g -= buf;
  return g;
} 

/** \brief Return the Hessian of the log-pseudolikelihood
 */
CatMatrix PseudoLikDiscrete::GetHessian(CatVector theta) {
  CatMatrix H = 0.0*outer(theta,theta); // null matrix
  for (unsigned int i=0;i<stat_splits.size();i++)
    H -= exp(sum(theta*stat_splits[i]))*outer(stat_splits[i],stat_splits[i]);
  // the sum is normalized in order to approximate the integral
  TTessel* tes = engy->get_ttessel();
  H =  tes->get_total_internal_length()/(CGAL_PI*stat_splits.size())*H;
  for (unsigned int i=0;i<stat_flips.size();i++)
    H -= exp(sum(theta*stat_flips[i]))*outer(stat_flips[i],stat_flips[i]);
  return H;
}
/** \brief Send T-tessellation pseudolikelihood data to an output stream
 *
 * \param os : output stream
 * \param pld : the pseudolikelihood of a T-tessellation
*/
std::ostream& operator<<(std::ostream &os, const PseudoLikDiscrete &pld) {
  os << "sum of merge statistics: " << pld.sum_stat_merges << std::endl;
  os << "split statistics" << std::endl;
  for (unsigned int i=0; i!=pld.stat_splits.size();i++) {
    os << pld.stat_splits[i] << std::endl;
  }
  os << "flip statistics" << std::endl;
  for (unsigned int i=0; i!=pld.stat_flips.size();i++) {
    os << pld.stat_flips[i] << std::endl;
  }
  os << "sum of flip statistics: " << pld.sum_stat_flips << std::endl;
  return os;
}
/******************************************************************************/
/*                    METHODS FOR PSEUDOLIKNOIS CLASS                         */
/******************************************************************************/

/** \brief Generate a PLInferenceNOIS object 
 * \param eng : a pointer to an Energy object defining the model to be fitted
 * \param store : whether the intermediate estimates and log-pseudolikelihood
 * approximations should be stored
 * \param lambda : stepsize to be used when updating the estimate
 */
PLInferenceNOIS::PLInferenceNOIS(Energy* eng, bool store, 
					   double lambda) : 
  PseudoLikDiscrete(eng), store_path(store), stepsize(lambda) {
  // add as many dummy splits as there are possible merges
  AddSplits(GetEnergy()->get_ttessel()->all_merges().size());
  if (store_path) {
    estimates.push_back(eng->get_theta());
    values.push_back(GetValue(eng->get_theta()));
  }
}

/** \brief Proceed to several iterations of the estimation algorithm
 *
 * \param n : number of iterations
 *
 */
void PLInferenceNOIS::Step(unsigned int n) {
  Energy* e = GetEnergy();
  for (unsigned int i=0;i!=n;i++) {
    CatVector currTheta = e->get_theta();
    CatVector currGrad = GetGradient(currTheta);
    /*if (previous_theta.IsEmpty()) { // first (Newton) move*/
      CatMatrix H = GetHessian(currTheta);
      FMatrix FH = asFMatrix(H);
      FVector FG = asFVector(currGrad);
      double one = 1.0;
      FMatrix FHinv = LA::inverse(FH,one);
      FVector Feps = -FHinv*FG;
      CatVector shift(currTheta);
      fill(shift,Feps);
      e->set_theta(currTheta+stepsize*shift);
      // Barzilai and Borwein method for computing step size
      // sometimes fails and yield oscillating trajectories. Instead,
      // we use the Newton method.
      /*} else { // standard Barzilai and Borwein move
      CatVector deltaTheta = currTheta-previous_theta;
      CatVector previous_grad = GetGradient(previous_theta);
      CatVector deltaGrad = currGrad-previous_grad;
      double epsilon = -sum(deltaTheta*deltaGrad)/sum(deltaGrad*deltaGrad);
      std::cout << "debug: epsilon " << epsilon << std::endl;
      e->set_theta(currTheta+epsilon*GetGradient(currTheta));
      }*/
    previous_theta = currTheta;
    if (store_path) {
      estimates.push_back(e->get_theta());
    }
    AddSplits(e->get_ttessel()->all_merges().size());
    if (store_path) {
      values.push_back(GetValue(e->get_theta()));
    }
  }
}
/** \brief Access to the current parameter estimate
 */
CatVector PLInferenceNOIS::GetEstimate() {
  Energy* e = GetEnergy();
  CatVector param = e->get_theta();
  return param;
}
/** \brief Iterate the estimation algorithm until convergence
 *
 * \param tol : target relative precision. Must be between 0 and 1.
 * \param nmax : maximal number of iterations
 *
 * \return the number of performed iterations
 *
 * Proceed to iterations of the inference algorithm (see documentation
 * of method Step) until the stopping criterion is met. One stops when
 * either nmax iterations have been performed or when the
 * log-pseudolikelihood approximation (lp) has been reduced by a
 * factor less than tol*(abs(lp)+tol).
 */
unsigned int PLInferenceNOIS::Run(double tol, unsigned int nmax) {
  if (nmax==0)
    return 0;
  unsigned int i = 0;
  double current_lpl;
  do {
    previous_lpl = GetValue(GetEstimate());
    Step(1);
    current_lpl = GetValue(GetEstimate());
    i++;
  } 
  while (fabs(current_lpl-previous_lpl)>=tol*(fabs(current_lpl)+tol) &&
	 i<nmax);
  return i;
}

/******************************************************************************/
/*                           FUNCTION BODIES                                  */
/******************************************************************************/

/** \brief Test whether three points are aligned
 * \param p : first point.
 * \param q : intermediate point.
 * \param r : last point.
 * \param verbose : if true, when points are found not aligned, a 
 *                  message explaining whny is printed to the standard 
 *                  output for logging. Default to false.
 * \return true if the points are aligned, false otherwise.
 * If q is not strictly between p and r, the points are considered as not aligned.
 */
bool are_aligned(Point2 p, Point2 q, Point2 r,bool verbose) {
  if (!CGAL::collinear(p,q,r)) {
    if (verbose)
      std::clog << "are_aligned function: points are not collinear" << std::endl;
    return false;
  }
  if (p==q || q==r) {
    if (verbose)  {
      std::clog << "are_aligned function: the intermediate point is equal ";
      std::clog << "to an end point" << std::endl;
    }
    return false;
  }
  if (CGAL::angle(p,q,r)==CGAL::ACUTE) {
    if (verbose) {
      std::clog << "are_aligned function: the intermediate point is not ";
      std::clog << "between the end points" << std::endl;
    }
    return false;
  }
  return true;
}
/** \brief Clip a segment by a convex polygon
 * \param S : segment to be clipped.
 * \param P : clipping polygon.
 * \return clipped segment, that is the intersection between the segment and
 * the polygon.
 * \pre The clipping polygon must be convex.
 */
Segment clip_segment_by_convex_polygon(Segment S, Polygon P) {
  typedef CGAL::Extended_cartesian<NT> EKernel;
  typedef CGAL::Nef_polyhedron_2<EKernel> Nef;
  std::vector<Nef::Point> Pvertices, Svertices;
  for (Size i=0;i!=P.size();++i) 
    Pvertices.push_back(Nef::Point(P[i].x(),P[i].y()));
  Nef nefP(Pvertices.begin(),Pvertices.end());
  Svertices.push_back(Nef::Point(S[0].x(),S[0].y()));
  Svertices.push_back(Nef::Point(S[1].x(),S[1].y()));
  Nef nefS(Svertices.begin(),Svertices.end());
  Nef nefInter = nefS.intersection(nefP);
  std::vector<Point2> Intervertices;
  Nef::Explorer e = nefInter.explorer();
  for (Nef::Explorer::Vertex_const_iterator v=e.vertices_begin();
       v!=e.vertices_end();v++) {
    if (e.is_standard(v)) {
      Point2 p(e.point(v).x(),e.point(v).y());
      Intervertices.push_back(p);
    }
  }
  return Segment(Intervertices[0],Intervertices[1]);
}
/** \brief Clip a segment by a polygon
 * \param S : segment to be clipped.
 * \param P : clipping polygon.
 * \return clipped segment, that is the longest intersection between the segment and
 * the polygon.
 * 
 * When the clipping polygon is not convex, its intersection with the segment
 * may consists of several segments. In such a case, the longest clipped segment is returned. If there are several clipped segments with maximal length, an arbitrarily chosen one is returned.
 */
Segment clip_segment_by_polygon(Segment S, Polygon P) {
  typedef CGAL::Extended_cartesian<NT> EKernel;
  typedef CGAL::Nef_polyhedron_2<EKernel> Nef;
  std::vector<Nef::Point> Pvertices, Svertices;
  for (Size i=0;i!=P.size();++i) 
    Pvertices.push_back(Nef::Point(P[i].x(),P[i].y()));
  Nef nefP(Pvertices.begin(),Pvertices.end());
  Svertices.push_back(Nef::Point(S[0].x(),S[0].y()));
  Svertices.push_back(Nef::Point(S[1].x(),S[1].y()));
  Nef nefS(Svertices.begin(),Svertices.end());
  Nef nefInter = nefS.intersection(nefP);
  Segment longestSegment;
  Nef::Explorer ex = nefInter.explorer();
  NT maxLength = 0;
  for (Nef::Explorer::Halfedge_const_iterator e=ex.halfedges_begin();
       e!=ex.halfedges_end();e++) {
    if (ex.is_standard(ex.source(e)) && ex.is_standard(ex.target(e))) {
      Point2 p1 = Point2(ex.point(ex.source(e)).x(),
			 ex.point(ex.source(e)).y());
      Point2 p2 = Point2(ex.point(ex.target(e)).x(),
			 ex.point(ex.target(e)).y());
      NT len = CGAL::squared_distance(p1,p2);
      if (len>maxLength) {
	maxLength = len;
	longestSegment = Segment(p1,p2);
      }
    }
  }
  return longestSegment;
}
/** \brief Predict added length when lengthening an edge in an arrangement
 * \param[in] e : halfedge that would be lengthened.
 * \param[out] ehit : halfedge where lengthening would end.
 * \param[out] p : location of the end of lengthening.
 * \return added length.
 *
 * Predict what would happen if an edge in an arrangement is extended
 * to the next existing edge.
 */
 double precompute_lengthening(Arrangement::Halfedge_handle e,
			       Arrangement::Halfedge_handle* ehit,
			       Point2* p) {
   Arrangement::Vertex_handle v = e->target();
   Rayon r(v->point(),Vector(e->source()->point(),e->target()->point()));
   // find the halfedge hit by ray r and the intersection
   /* Not that easy. The target halfedge may be on the outer ccb
      or on a "hole" boundary. Approach below: visit all
      halfedges (on ccb or along holes), compute intersection, if
      any keep the closest one. */
   // First collect all halfedges
   Arrangement::Face_handle f = e->face();
   std::vector<Arrangement::Halfedge_handle> all_edges;
   Arrangement::Ccb_halfedge_circulator hc=f->outer_ccb();
   Arrangement::Halfedge_handle h0 = hc;
   do {
     // ignore the halfedge to be extended
     if (hc->source()!=v && hc->target()!=v)
       all_edges.push_back(hc);
     hc++;
   } while (hc!=h0);
   for (Arrangement::Hole_iterator ho=f->holes_begin();
	ho!=f->holes_end();ho++) {
     hc = *ho;
     h0 = hc;
     do {
       if (hc->source()!=v && hc->target()!=v)
	 all_edges.push_back(hc);
       hc++;
     } while (hc!=h0);
   } // end collection of all halfedges
   double lp = std::numeric_limits<double>::max();
   Point2 inter_location;
   for (std::vector<Arrangement::Halfedge_handle>::iterator ei=all_edges.begin();
	ei!=all_edges.end();ei++) {
     CGAL::Object inter = CGAL::intersection(r,
					     Segment((*ei)->source()->point(),
						     (*ei)->target()->point())); 
     if (CGAL::assign(inter_location,inter)) {
       double len = sqrt(CGAL::to_double(CGAL::squared_distance(v->point(),
								inter_location)));
       if (len<lp) {
	 *p = inter_location;
	 *ehit = *ei;
	 lp = len;
       }
     }
   }
   return lp;
 }

/** \brief Compute the previous halfedge along a segment
 *
 * \param e: halfedge whose neighbour is to be computed
 * \return previous neighbour, NULL_HALFEDGE_HANDLE if e is the first halfedge on 
 * its segment.
 */
LineTes::Halfedge_handle compute_prev_hf(LineTes::Halfedge_handle e) {
  
  LineTes::Halfedge_around_vertex_circulator hf_circ = 
    e->source()->incident_halfedges();
  LineTes::Halfedge_handle e_prev=NULL_HALFEDGE_HANDLE;
  unsigned int i;
  
  for (i=0;i<e->source()->degree();i++) {
    if (hf_circ != e->twin()) {
      if (CGAL::collinear(e->target()->point(),e->source()->point(),
			  hf_circ->source()->point()) &&
	  CGAL::angle(e->target()->point(),e->source()->point(),
		      hf_circ->source()->point())==CGAL::OBTUSE) {
	e_prev = hf_circ;
	break;
      } 	    	
    }
    hf_circ++;
  } 
  
  return e_prev;
}

/** \brief Compute the next halfedge along a segment
 *
 * \param e: halfedge whose neighbour is to be computed
 * \return next neighbour, NULL_HALFEDGE_HANDLE if e is the last halfedge on 
 * its segment.
 */
LineTes::Halfedge_handle compute_next_hf(LineTes::Halfedge_handle e) {

  LineTes::Halfedge_around_vertex_circulator hf_circ = 
    e->target()->incident_halfedges();
  LineTes::Halfedge_handle e_next=NULL_HALFEDGE_HANDLE;
  unsigned int i;

  for (i=0;i<e->target()->degree();i++) {
    if (hf_circ != e) {
      if (CGAL::collinear(e->source()->point(),e->target()->point(),
			  hf_circ->source()->point()) &&
	  CGAL::angle(e->source()->point(),e->target()->point(),
		      hf_circ->source()->point())==CGAL::OBTUSE) {
	e_next = hf_circ->twin();
	break;  
      }
    }
    hf_circ++;
  } 
  return e_next;
}

/** \brief Set neighbour relationships between two halfedges
 * \param e1 : halfedge.
 * \param e2 : other halfedge supposed to follow e1 on the same segment.
 * \note{Either e1 or e2 can be null.}
 */
void set_junction(LineTes::Halfedge_handle e1,LineTes::Halfedge_handle e2) {
  /* Mise à jour de la relation de voisinage entre 2 demi-arêtes contigues 
     et leurs jumelles.
     e1 : première demi-arête.
     e2 : deuxième demi-arête.
     Pré-requis : e2 doit succéder à e1.
     Note : l'une des demi-arêtes peut être "nulle". */
  if (e1!=NULL_HALFEDGE_HANDLE) {
    e1->set_next_hf(e2);
    if (e2!=NULL_HALFEDGE_HANDLE)
      e1->twin()->set_prev_hf(e2->twin());
    else
      e1->twin()->set_prev_hf(NULL_HALFEDGE_HANDLE);
  }
  if (e2!=NULL_HALFEDGE_HANDLE) {
    e2->set_prev_hf(e1);
    if (e1!=NULL_HALFEDGE_HANDLE)
      e2->twin()->set_next_hf(e1->twin());
    else
      e2->twin()->set_next_hf(NULL_HALFEDGE_HANDLE);
  }
}

/** \brief Test whether a halfedge handle is valid
 *
 * Just check whether the halfedge handle can be found within
 * the halfedge iterators.
 * \param tesl : the line tessellation to be searched.
 * \param e : halfedge to be checked.
 */
bool exist_halfedge(LineTes& tesl,LineTes::Halfedge_handle e) {
  bool exists = false;
  for (LineTes::Halfedge_iterator et=tesl.halfedges_begin();
       et!=tesl.halfedges_end();et++) {
    if (et==e) {
      exists = true;
      break;
    }
  }
  return exists;
}
/** \brief Find a tessellation halfedge ending at two given points
 *
 * \param tesl : a tessellation.
 * \param p1 : a point.
 * \param p2 : another point.*/
LineTes::Halfedge_handle find_halfedge(LineTes& tesl,Point2 p1,Point2 p2) {
  LineTes::Halfedge_iterator et;
  bool found = false;
  for (et=tesl.halfedges_begin();
       et!=tesl.halfedges_end();et++) {
    if (et->source()->point()==p1 && et->target()->point()==p2) {
      found = true;
      break;
    }
  }
  if (found)
    return et;
  else
    return NULL_HALFEDGE_HANDLE;
}

/*
unsigned long int number_of_internal_vertices(TTessel& t)
{
  TTessel::Face_handle f=t.unbounded_face();
  TTessel::Ccb_halfedge_circulator e;
 

  TTessel::Halfedge_handle e_it;

  for (e_it= t.halfedges_begin();e_it!=t.halfedges_end();e_it++){
    TTessel::Halfedge_const_handle e_curr=e_it;
    if (f->is_halfedge_on_inner_ccb(e_curr)) break;
  }

  e = e_it->ccb();

 do{
   compt ++;
   e++;
 }
 while(e!=e_it);

 return t.number_of_vertices()-compt;
}
*/

/** \brief Test whether a vertex of a line tessellation is a T-vertex
 *  \param v : vertex to be tested.
 *  \param verbose : if true, information about I-vertex processing is 
 * sent to the standard output. Default to false.
*/
bool is_a_T_vertex(LineTes::Vertex_handle v, bool verbose) {
  if (v->degree()!=3) {
    if (verbose) {
      std::clog << "is_a_T_vertex function: vertex ";
      std::clog << v->point() << " not a T-vertex because ";
      std::clog << "its degree not equal to 3" << std::endl;
    }
    return false;
  }
  LineTes::Halfedge_around_vertex_circulator e = v->incident_halfedges();
  LineTes::Halfedge_handle e0 = e;
  int n = 0; // for storing the number of terminal incident halfedges
  do {
    if (e->get_next_hf()==NULL_HALFEDGE_HANDLE) {
      n++;
    }
    e++;
  } while(e!=e0);
  if (n!=1) {
    if (verbose) {
      std::clog << "is_a_T_vertex function: vertex ";
      std::clog << v->point() << " not a T-vertex because ";
      if (n==3) {
	std::clog << "it is a Y-vertex" << std::endl;
      } else {
	std::clog << "of an unexpected reason" << std::endl;
	std::clog << "number of terminal incident halfedges: " << n;
	std::clog << std::endl;
	std::clog << "incident halfedges:" << std::endl;
	do {
	  std::clog << "\t" << e->source()->point() << " ";
	  std::clog << e->target()->point() << " ";
	  if (e->get_next_hf()==NULL_HALFEDGE_HANDLE) {
	    std::clog << "terminal";
	  } else {
	    std::clog << "not terminal";
	  }
	  std::clog << std::endl;
	  e++;
	} while(e!=e0);
      }
    }
    return false;
  }
  return true;
}
/** \brief Return the total number of internal vertices
 * \param t : the T-tessellation.
 * \return the number of internal vertices.
 *
 * A vertex is internal if it does not lie on the boundary of
 * the tessellated domain.
 */
unsigned long int number_of_internal_vertices(TTessel& t){
  unsigned long int compt=0;
  	
  TTessel::Seg_list_iterator s=t.segments_begin();
  for (int i=0;i<t.number_of_window_edges();i++){
    compt += (*s)->number_of_edges();
    s++;
  }
  return t.number_of_vertices()-compt;

}

/******************************************************************************/
/*                         FONCTIONS CALCULATING FEATURES                     */
/******************************************************************************/

/** \brief Test whether a point is in the interior of the tessellated domain
 * \return a double resulting from conversion from a Boolean.
 */
double is_point_inside_window(Point2 pt,LineTes* t){
  return (double)t->get_window().has_on_bounded_side(pt);
}
/** \brief Test whether a point is in the interior of the tessellated domain
 * \return a double resulting from conversion from a Boolean.
 */
double is_point_inside_window(Point2 pt,TTessel* t){
  // repeated code, any way to avoid that without slowing computation too much?
  return (double)t->get_window().has_on_bounded_side(pt);
}


/** \brief Return the length of an internal T-tessellation segment
 *
 * If the segment is not internal, zero is returned.
 */ 
double edge_length(Segment seg, TTessel* t){
  Point2 p1=seg.vertex(0), p2=seg.vertex(1);
  std::vector<Point2> segVertices(2);
  segVertices[0] = p1; segVertices[1] = p2;
  if (is_segment_internal(segVertices,t)>0)
    return 0;
  else 
    return sqrt(CGAL::to_double(seg.squared_length()));
}
/** \brief Return 1
 *
 * \param s : a tessellation segment as a vector of its 2 ends.
 * \param t : the tessellation to be considered.
 *
 * Silly function that can be used by an Energy object for specifying
 * the number of tessellation segments as a feature.*/
double seg_number(std::vector<Point2> s, TTessel* t){
  return 1;
}

/** \brief Test whether a segment of a T-tessellation is internal
 *
 * \param s : a vector containing the two ends of the segment to be
 *            tested. The ends must be Point2 objects.
 * \param t : a pointer to the T-tessellation.
 * \return 0 if the segment is internal, 1 otherwise.
 * \pre the segment must be included in the domain of the T-tessellation.
 *
 * A segment is considered as internal if it is not fully contained 
 * in the domain boundary.
 */
double is_segment_internal(std::vector<Point2> s, TTessel* t){
  Point2 p1 = s.front(), p2 = s.back();
  if (is_point_inside_window(p1,t)>.5 || is_point_inside_window(p2,t)>.5)
    return 1.0;
  // Both segment ends are on the boundary
  Polygon dom = t->get_window();
  // Test whether both lie on the same domain side
  for (Polygon::Edge_const_iterator e=dom.edges_begin();
       e!=dom.edges_end();e++) {
    if (e->has_on(p1) && e->has_on(p2)) {
      return 0.0;
    }
  }
  return 1.0;
}

/** \brief Return 1
 *
 * \param f: a polygon representing a tessellation face. 
 * \param t: the tessellation to be considered.
 * \ingroup features
 *
 * Silly function that can be used by an Energy object for 
 * specifying the number of tessellation face as a feature.*/
double face_number(Polygon f, TTessel* t){
  return 1;
}

/** \brief Return the squared face area of a tessellation face
 *
 * \param f : tessellation face as a polygon.
 * \param t : the tessellation to be considered.
 * \ingroup features*/
double face_area_2(Polygon f, TTessel* t){
    return to_double(f.area()*f.area());
}
/** \brief Return the perimeter of a tessellation face
 *
 * \param f : tessellation face as a polygon.
 * \param t : the tessellation to be considered.
 * \ingroup features*/ 
double face_perimeter(Polygon f, TTessel* t){
  
  double l=0;

  for (int i=0;i!=f.size();i++){
    l += sqrt(CGAL::to_double(squared_distance(f[i],f[(i+1)%f.size()])));
  }
  return l;
}
/** \brief Return the shape parameter of a tessellation face.
 *
 * \param f : tessellation face as a polygon
 * \param t : the tessellation to be considered
 *
 * The shape parameter is the ratio between the perimeter
 * and 4 times the squared root of the area. */
double face_shape(Polygon f, TTessel* t){
  return face_perimeter(f,t)/(4*sqrt(CGAL::to_double(f.area()))); 
}

/** \brief Return the angle between two planar vectors
 *
 * \param v1 : first vector.
 * \param v2 : second vector.*/
double angle_between_vectors(Vector v1,Vector v2){
  NT n1_2 = v1*v1;
  NT n2_2 = v2*v2;
  if (v1*v1==0 || v2*v2==0) {
    std::cerr << "angle_between_vectors: one of the vectors is null" << std::endl;
  }
  NT p =  v1*v2;
  if (p*p==n1_2*n2_2)
    if (p>0)
      return 0.;
    else
      return CGAL_PI;
  return  acos(CGAL::to_double(p)/sqrt(CGAL::to_double(n1_2*n2_2)));
}
/** \brief Measure the deviation of a T-tessellation face from a rectangle
 * \param f : T-tessellation face as a polygon.
 * \param t : the T-tessellation to be considered.
 * \ingroup features
 * 
 * Only vertices where incident edges form an acute angle contribute 
 * to the measure. The additive contribution of such a vertex \f$v\f$ 
 * is
 * \f[
 * \frac{\pi}{2}-\phi(v),
 * \f]
 * where \f$\phi(v)\f$ is the (acute) angle between the edges incident
 * to \f$v\f$. The above contribution is close to zero if the incident
 * edges are almost perpendicular. It is close to \f$\pi/2\f$ if the 
 * edges are almost aligned. 
*/
double face_sum_of_angles(Polygon f, TTessel* t){
  Vector v1,v2;
  Point2 p;
  double tot_angle=0;

  for(int i=0;i!=f.size();i++) {
    int j=(i-1)%((int)f.size());
    if (j<0) j=j+f.size();
    p=f[i];
    //if (is_point_inside_window(p,t)>.5)
      {
	v1 = Vector(f[i],f[j]);
        v2 = Vector(f[i],f[(i+1)%((int)f.size())]);
	double a = CGAL_PI/2-angle_between_vectors(v1,v2);
	if (a>=0)
	  tot_angle += a;
      }
}
  return tot_angle;
}

/** \brief Returns a Polygon object representing a face of the
 *     tessellation. 
 *
 * Only vertices joining non-aligned edges are specified during the
 * Polygon construction.
 **/
Polygon face2poly(TTessel::Face_handle f)
{
  Polygon poly;
  TTessel::Ccb_halfedge_circulator e;
  
  e =f->outer_ccb();
  TTessel::Halfedge_handle e_first=e;

  do {
    if (e->get_next_hf()==NULL_HALFEDGE_HANDLE || (e->get_next_hf()!=e->next())) 
    poly.push_back(e->target()->point());
    e++;
  }
  while (e!=e_first);

    return(poly);
}
/** \brief Return the smallest angle in a tessellation face
 *
 * \param f : tessellation face as a polygon.
 * \param t : the tessellation to be considered.
 * \ingroup features
 *
 * Angles between consecutive edges are considered. This function 
 * can be used by an Energy object for specifying the sum of face
 * smallest angle as a feature. That feature is of interest for 
 * penalizing faces with small angles.
*/
double min_angle(Polygon f, TTessel* t){
  Vector v1,v2;
  Point2 p;
  double min_angle;
 
  min_angle = angle_between_vectors(Vector(f[f.size()-1],f[0]),
				    Vector(f[0],f[1]));
  for(int i=0;i!=f.size();i++) {
    int j=(i-1)%((int)f.size());
    if (j<0) j=j+f.size();
    p=f[i];
    v1 = Vector(f[i],f[j]);
    v2 = Vector(f[i],f[(i+1)%((int)f.size())]);
    double a = angle_between_vectors(v1,v2);
    if (a<min_angle)
      min_angle = a;
  }
  return min_angle;
}
/** \brief Return the sum of squared areas of all faces of a T-tessellation
 *
 * \param t : the tessellation to be considered.
 */
double sum_of_faces_squared_areas(TTessel* t){
double sa = 0.;
Polygon poly;
      for (TTessel::Face_iterator f = t->faces_begin();f!=t->faces_end();f++) {
      	if (!f->is_unbounded()) { 
        	poly = face2poly(f);
	    	sa += face_area_2(face2poly(f),t);
      	}
       }
return(sa);
}
/** \brief Return the sum of smallest angles on all faces of a T-tessellation
 *
 * \param t : the tessellation to be considered.
 */
double sum_of_min_angles(TTessel* t){
  double sma = 0.;
  Polygon poly;
  for (TTessel::Face_iterator f = t->faces_begin();f!=t->faces_end();f++) {
    if (!f->is_unbounded()) { 
      poly = face2poly(f);
      sma += min_angle(face2poly(f),t);
    }
  }
  return(sma);
}
/** \brief Return the sum of obtuse angles on all faces of a T-tessellation
 *
 * \param t : the tessellation to be considered.
 */ 
double sum_of_angles_obt(TTessel* t){
  double sao = 0.;
  Polygon poly;
  for (TTessel::Face_iterator f = t->faces_begin();f!=t->faces_end();f++) {
    if (!f->is_unbounded()) { 
      poly = face2poly(f);
      sao +=face_sum_of_angles(poly,t);
    }
  }
  return(sao);
}
/** \brief Return the squared number of edges on a segment of a 
 * T-tessellation
 *
 * \param s :  a T-tessellation segment as a vector of its vertices.
 * \param t : the T-tessellation to be considered.
 * \ingroup features
 *
 * Function that can be used by an Energy object for 
 * specifying the sum of squared number of edges on all segments as
 * a feature. This feature can be used for controlling the variability
 * of segment sizes (here number of edges per segment).*/
double segment_size_2(std::vector<Point2> s, TTessel* t){
  /* Return the squared number of edges on the segment if internal,
     otherwise 0. */
  Point2 p1 = s.front(), p2 = s.back();
  std::vector<Point2> segVertices(2);
  segVertices[0] = p1; segVertices[1] = p2;
  if (is_segment_internal(segVertices,t)>0)
    return 0;
  else 
    return (s.size()-1)*(s.size()-1);
}
/** \brief Return the sum of squared numbers of edges on all segments of a
 * T-tessellation
 *
 * \param t : the tessellation to be considered.
 */
double sum_of_segment_squared_sizes(TTessel* t) {
  double s4 = 0;
  for (TTessel::Seg_list_iterator si=t->segments_begin();si!=t->segments_end();si++)
    s4 += segment_size_2((*si)->list_of_points(),t);
  return s4;
}

// FIN KA : 26-06-06