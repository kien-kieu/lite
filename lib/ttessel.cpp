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
  */
LineTes::LineTes() : Arrangement() {
  unbounded_face()->set_data(false);
}
/** \brief Define the the tessellation domain as a rectangle
 *
 * The insert_window method should be used as a first step when
 * handling line tessellations.
 */
void LineTes::insert_window(Rectangle r) {
  Polygon p;
  for (int i=0;i<4;i++) {
    p.push_back(r[i]);
  }
  insert_window(p);
}
/** \bried Define the tessellation domain as a polygon without hole
 */
void LineTes::insert_window(Polygon& p) {
  HPolygon hp(p);
  HPolygons hps;
  hps.push_back(hp);
  insert_window(hps);
}
/** \brief Define the tessellation domain
 *
 * The tessellation domain can be a union of disjoint polygons with
 * holes.
 */
void LineTes::insert_window(HPolygons& p) {
  Halfedge_handle e;
  Vertex_handle v,v0;
  Face_handle f;

  // Reset private member window
  window.clear();
  for (HPolygons::const_iterator pi=p.begin();pi!=p.end();pi++) {
    HPolygon sp = simplify(*pi);
    window.push_back(sp);
  }
  
  /* Now insert the window edges in the arrangement */
  for (HPolygons::iterator pi=window.begin();pi!=window.end();pi++) {
    Polygons b = boundaries(*pi);
    for (Polygons::iterator bi=b.begin();bi!=b.end();bi++) {
      v = insert_point(*this,(*bi)[0]);
      v0 = v;
      for (int i=0;i<bi->size()-1;i++) {
	if ((*bi)[i]<=(*bi)[i+1]) {
	  e = insert_from_left_vertex(Curve((*bi)[i],(*bi)[i+1]),v);
	} else {
	  e = insert_from_right_vertex(Curve((*bi)[i+1],(*bi)[i]),v);
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
      e = insert_at_vertices(Curve((*bi)[bi->size()-1],(*bi)[0]),v,v0);
      v = e->target();
      e->set_length();
      e->twin()->set_length(e->get_length());
      set_junction(e,NULL_HALFEDGE_HANDLE);
      set_junction(NULL_HALFEDGE_HANDLE,e);
      e->set_dir(true);
      e->twin()->set_dir(false);
      if (bi->area()<0) { 
	// the inserted polygon is a hole
	e->twin()->face()->set_data(false);
      } else { // the inserted polygon is an outer boundary
	e->face()->set_data(true);
      }
      Seg_handle s = new Seg;
      s->set_halfedge_handle(e);
      e->set_segment(s); e->twin()->set_segment(s);
      all_segments.add(s);
    }
  }
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

  // Both faces apart new edge are inside the domain
  new_e->face()->set_data(true);
  new_e->twin()->face()->set_data(true);

  return new_e;
}
/** \brief Split a tessellation cell from a vertex to an edge
 * \param v : handle of the vertex where the splitting segment starts.
 * \param e : handle of the edge where the splitting segment ends.
 * \param p: point where the splitting segment ends.
 * \return  handle of one of the halfedges along the splitting segment.
 * \pre e should not bound the external face. If so, return 
 * NULL_HALFEDGE_HANDLE.
 */
LineTes::Halfedge_handle LineTes::split_from_vertex(Vertex_handle v,
						    Halfedge_handle e,
						    Point2 p) {
  if (is_on_boundary(e)==2) 
    return NULL_HALFEDGE_HANDLE;

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
  // Both faces apart the new edge are inside the domain
  new_e->face()->set_data(true);
  new_e->twin()->face()->set_data(true);

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
  f->set_data(true); // new face inside the domain
  return f;
}
/** \brief Clear a tessellation
 * \param remove_window : if true (default), the domain is removed. Otherwise, 
 * only internal segments are removed. 
 */
void LineTes::clear(bool remove_window) {
  Arrangement::clear();
  all_segments.clear();
  unbounded_face()->set_data(false);
  if (remove_window) {
    window.clear();
  } else {
    HPolygons dom = get_window();
    insert_window(dom);
  }
}
/** \brief Return the number of internal segments in a line tessellation
 *
 * A segment is internal if it does no lie along the domain boundary.
 */
Size LineTes::number_of_internal_segments() {
  Size count = 0;
  for (Seg_list_iterator si=segments_begin();si!=segments_end();si++) {
    if (!is_on_boundary(*si)) count++;
  }
  return count;
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
  Size nb = 0;
  for (HPolygons::const_iterator hpi=window.begin();hpi!=window.end();hpi++) {
    Polygons b = boundaries(*hpi);
    for (Polygons::const_iterator pi=b.begin();pi!=b.end();pi++) {
      nb += pi->size();
    }
  }
  return nb;
}
/** \brief Return the perimeter of the tessellated domain
 */
double LineTes::get_window_perimeter() { // ???
  double per=0.0;
  for (HPolygons::const_iterator hpi=window.begin();hpi!=window.end();hpi++) {
    Polygons b = boundaries(*hpi);
    for (Polygons::const_iterator pi=b.begin();pi!=b.end();pi++) {
      for(Polygon::Edge_const_iterator e=pi->edges_begin();e!=pi->edges_end();
	  e++) {
	per += sqrt(CGAL::to_double(CGAL::squared_distance((*e)[0],(*e)[1])));
      }
    }
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
 * - 1 if the halfedge bounds an internal face.
 * - 2 if the halfedge bounds the external face.
 */
int LineTes::is_on_boundary(Halfedge_handle e){
  Face_handle f_e,f_e_twin;

  f_e = e->face();
  f_e_twin = e->twin()->face();

  if (f_e->data() && f_e_twin->data()) 
    return 0;
  else { // One of the faces bounded by the edge is inside the domain
    if (f_e->data()) // The face bounded by the halfedge is inside
      return 1;
    else  // The face bounded by the halfedge is outside
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
  HPolygons w = get_window();
  Segment cseg = clip_segment_by_polygon(iseg,w);
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
	Face_handle f = e->face();
	Halfedge_handle esplit;
	ivertex.p = ray_exit_face(r,f,esplit);
	ivertex.esplit = esplit;
	double lp = sqrt(CGAL::to_double(CGAL::squared_distance(v->point(),
								ivertex.p)));

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
 * - A coordinate is output as two integers (numerator and denominator
 *   separated by a slash).
 * - A point is output as its two Cartesian coordinates.
 * - A line segment is output as its two end points.
 * - First, vertices of the tessellated domain are sent to the output 
 *   stream. For each holed polygon of the domain, vertices of the outer 
 *   boundary are output anti-clockwise (single line). Then vertices 
 *   along holes are output clockwise (one line per hole). 
 * - Subsequently, internal segments are sent to the output stream
 *   (one per line).
 */	
void LineTes::write(std::ostream& os) {
  HPolygons w = get_window();
  for (Size i=0;i<w.size();i++) {
    Polygons b = boundaries(w[i]);
    for (Size j=0;j!=b.size();j++) {
      for (Size k=0;k!=b[j].size();k++) {
      os << b[j][k].x().exact() << " " << b[j][k].y().exact();
      if (k!=(b[j].size()-1))
	os << " ";
      }
      os << std::endl;
    }
  }
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
  HPolygons win;
  HPolygon wel; // window element
  Polygon weloub; // window element outer boundary
  Polygons hb; // holes
  bool got_first_wel = false;
  bool done_with_domain = false;
  std::vector<Segment> segs;
  std::string line;
  // Clear the LineTes object
  clear();
  // Read the coordinates of the domain corners
  while (std::getline(is,line)) {
    std::istringstream iss(line);
    std::vector<NT> wcoords;
    CGAL::Gmpq buf;
    while (iss>>buf) {
	NT coord(buf);
	wcoords.push_back(coord);
    }      
    if (wcoords.size()%2 != 0) {
      std::cerr << "Uneven number of coordinates on the current line";
      std::cerr << std::endl;
      return;
    }
    if (wcoords.size()>4) { // coordinates of domain vertices
      Polygon ccb; // a boundary (internal or external)
      for (Size i=0;i!=wcoords.size();i+=2) {
	ccb.push_back(Point2(wcoords[i],wcoords[i+1]));
      }
      if (ccb.area()>0) { // outer ccb
	if (got_first_wel) { 
	  /* An outer boundary and holes have been read already. They
	     are stored in weloub and hb. Build a holed polygon and push
	     it in window. */
	  HPolygon wel(weloub,hb.begin(),hb.end()); // window element
	  win.push_back(wel);
	} else {
	  got_first_wel = true;
	}
	hb.clear();
	weloub = ccb;
      } else { 
	if (ccb.area()<0) { // hole
	  hb.push_back(ccb);
	} else { // undefined orientation
	  std::cerr << "Invalid orientation of a domain ccb" << std::endl;
	}
      }
    } else { // segment
      if (!done_with_domain) {
	HPolygon wel(weloub,hb.begin(),hb.end()); // window element
	win.push_back(wel);
	insert_window(win);
	done_with_domain = true;
      }
      NT x0(wcoords[0]);
      NT y0(wcoords[1]);
      NT x1(wcoords[2]);
      NT y1(wcoords[3]);
      Segment s(Point2(x0,y0),Point2(x1,y1));
      insert_segment(s);
    }
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
    HPolygons p = lt.get_window();
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

/** \brief Define domain to be tessellated as a rectangle.
 */
void TTessel::insert_window(Rectangle r) {
  Polygon p;
  for (int i=0;i<4;i++) {
    p.push_back(r[i]);
  }
  insert_window(p);
}
/** \brief Define the domain to be tessellated as a polygon without hole.
 */
void TTessel::insert_window(Polygon& p) {
  HPolygon hp(p);
  HPolygons uhp;
  uhp.push_back(hp);
  insert_window(uhp);
}
/** \brief Define the domain to be tessellated.
 *
 * The domain is a union of holed polygons.
 */
void TTessel::insert_window(HPolygons& uhp) {
  LineTes::insert_window(uhp);
  int_length = 0;
  for (HPolygons::const_iterator hp=uhp.begin();hp!=uhp.end();hp++) {
    Polygons b = boundaries(*hp);
    for (Polygons::const_iterator p=b.begin();p!=b.end();p++) {
      for (Polygon::Edge_const_iterator e=p->edges_begin();
	   e!=p->edges_end();e++) 
	int_length += sqrt(CGAL::to_double(e->squared_length()));
    }
  }
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
  Seg_handle s_at_ends[2] = {e->next()->segment(),
			     e->twin()->next()->segment()};
  for (unsigned i=0;i<2;i++) {
    Seg_handle s = s_at_ends[i];
    if (is_on_boundary(s)==0 && !s->number_of_edges_is_greater_than(2)) {
      non_blocking_segments.suppress(s);
      blocking_segments.add(s);
    }
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
  if (is_on_boundary(s1)==0 && !s1->number_of_edges_is_greater_than(2)) { 
    non_blocking_segments.add(s1);  
    blocking_segments.suppress(s1);
  }

  s2 = e->twin()->next()->segment();
  if (is_on_boundary(s2)==0 && !s2->number_of_edges_is_greater_than(2)) {
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
 * \param remove_window : if true (default), the domain is removed too.
 * otherwise only internal segments are removed.
 */
void TTessel::clear(bool remove_window) {
  LineTes::clear(remove_window);
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
HPolygons TTessel::all_faces() {
  HPolygons res;
  for (Face_iterator f = faces_begin(); f!= faces_end(); f++) {
    if (f->data()) {
      HPolygon p = face2poly(f);
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
    if (!e->face()->data())
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
    if (!e->face()->data())
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
    if (e->face()->data())
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
    if (f->data())
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

  // TTessel::Halfedge_handle ee = e1->next();
  // Point2 pp;
  // NT minDist = -1;
  Line l_dir;
  if (l.has_on_negative_side(e1->target()->point())) {
    l_dir = l;
  } else {
    l_dir = l.opposite();
  }
  Rayon r(p1,l_dir);
  Halfedge_handle exit_edge;
  Face_handle f = e1->face();
  p2 = ray_exit_face(r,f,exit_edge);
  e2 = exit_edge;

  // bool p2_found(false);
  // while (ee!=e1) {
  //   inter = CGAL::intersection(Segment(ee->source()->point(),
  // 				       ee->target()->point()),r);
  //   std::cout << "Split constructor: computing intersection of ray with edge" << std::endl;//debug 
  //   if (CGAL::assign(pp,inter)) {
  //     std::cout << "Split constructor: ray hits the current edge" << std::endl;//debug
  //     NT distp1p2 = CGAL::squared_distance(p1,pp);
  //     if (minDist>distp1p2 || minDist==-1) {
  //     	minDist = distp1p2;
  //     	e2 = ee;
  //     	p2 = pp;
  // 	if (!p2_found) p2_found = true;
  //     }
  //   }
  //   ee = ee->next(); 
  //   std::cout << "Split constructor: advancing along e1 ccb" << std::endl;//debug
  // }
  // std::cout << "Split constructor: exit from loop along e1 ccb" << std::endl;//debug 
  // if (!p2_found)
  //   std::cerr << "Split constructor : intersection of l with e2 not found"
  // 	      << std::endl;
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
  TTessel::Face_handle del_face = get_e1()->face();
  HPolygon del_poly = face2poly(del_face);
  del_poly = simplify(del_poly);
  modifs.del_faces.push_back(del_poly);

  // Added faces
  /* 4 possible cases
     case 1: the splitting line segment connects two points along 
             the outer boundary.
     case 2: the splitting line segment connects two points along
             the same inner boundary (bounding a hole).
     case 3: the splitting line segment connects a point on the 
             outer boundary to a point on an inner boundary.
     case 4: the splitting line segment connects two points on two
             different inner boundaries.
  */
  Halfedge_handle he1(get_e1()), he2(get_e2());
  Point2 pt1(get_p1()), pt2(get_p2());
  TTessel::Face_handle f = he1->face();
  HPolygon fp = face2poly(f,false);

  Polygons fp_borders = boundaries(fp);
  PECirc pe1, pe2;
  Size i1 = fp_borders.size(), i2 = i1, j;
  i1 = find_edge_in_polygons(fp_borders, 
			     Segment(he1->source()->point(), 
				     he1->target()->point()), pe1);
  i2 = find_edge_in_polygons(fp_borders, 
			     Segment(he2->source()->point(), 
				     he2->target()->point()), pe2);
  HPolygons add_faces = hpolygon_insert_edge(fp_borders, i1, i2, pe1, pe2, 
					  pt1, pt2);
  for (Size i=0;i<add_faces.size();i++)
    modifs.add_faces.push_back(simplify(add_faces[i]));

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

  modifs.del_edges.push_back(Segment(get_e()->source()->point(),
				     get_e()->target()->point()));
  modifs.del_edges.push_back(Segment(p2(),e2()->target()->point()));
  modifs.del_edges.push_back(Segment(e2()->get_prev_hf()->source()->point(),
				     p2()));
  modifs.del_edges.push_back(Segment(p1(),e1()->target()->point()));
  modifs.del_edges.push_back(Segment(e1()->get_prev_hf()->source()->point(),
				     p1()));

  // Added edges 

  modifs.add_edges.push_back(Segment(e2()->get_prev_hf()->source()->point(),
				     e2()->target()->point()));
  modifs.add_edges.push_back(Segment(e1()->get_prev_hf()->source()->point(),
				     e1()->target()->point()));

  // Face modifications
  HPolygon f = face2poly(e->face(),false), f_twin, hpoly_buf;
  bool single_face = e->face()==e->twin()->face();
  Polygons f_borders = boundaries(f), f_twin_borders;
  PECirc pe, pe_twin;
  Size i = find_edge_in_polygons(f_borders, 
				 Segment (e->source()->point(),
					  e->target()->point()),pe), i_twin;
  hpoly_buf = simplify(f);
  modifs.del_faces.push_back(hpoly_buf);
  if (single_face) {
    i_twin = find_edge_in_polygons(f_borders, 
				      Segment (e->target()->point(),
					       e->source()->point()),pe_twin);
    hpoly_buf = hpolygon_remove_edge(f_borders,f_borders,
				     i,i_twin,
				     pe,pe_twin,
				     single_face);
    hpoly_buf = simplify(hpoly_buf);
    modifs.add_faces.push_back(hpoly_buf);    
  } else {
    f_twin = face2poly(e->twin()->face(),false);
    f_twin_borders = boundaries(f_twin);
    i_twin = find_edge_in_polygons(f_twin_borders, 
				      Segment (e->target()->point(),
					       e->source()->point()),pe_twin);
    hpoly_buf = simplify(f_twin);
    modifs.del_faces.push_back(hpoly_buf);
    hpoly_buf = hpolygon_remove_edge(f_borders,f_twin_borders,
				     i,i_twin,
				     pe,pe_twin,
				     single_face);
    hpoly_buf = simplify(hpoly_buf);
    modifs.add_faces.push_back(hpoly_buf);
  }
   
  // segments
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
  Rayon r(e3->target()->point(),
	  Vector(e3->source()->point(),e3->target()->point()));
  // Find the face where the new edge is to be inserted
  Face_handle f;
  if (e1->face()==e1->get_prev_hf()->face()) {
    f = e1->face();
  } else {
    f = e1->twin()->face();
  }
  // Set up private members e2 (split edge) and p2 (new vertex on e2)
  p2 = ray_exit_face(r,f,e2);
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

  // Faces
  HPolygon f2 = face2poly(get_e2()->face(),false);
  Polygons f2_borders = boundaries(f2);
  PECirc pe_p1, pe2; // pe_p1 edge of f2 containaing p1
  bool f1_equals_f1_minus = get_e1()->face()==get_e1()->get_prev_hf()->face();
  Segment seg_buf(get_e2()->source()->point(),get_e2()->target()->point());
  Size i2 = find_edge_in_polygons(f2_borders,seg_buf,pe2);
  if (f1_equals_f1_minus) {
    seg_buf = Segment(get_e1()->source()->point(),
		      get_e1()->target()->point());
  } else {
    seg_buf = Segment(get_e1()->twin()->source()->point(),
		      get_e1()->twin()->target()->point());
  }
  Size i1 = find_edge_in_polygons(f2_borders,seg_buf,pe_p1);
  HPolygons del_faces;
  del_faces.push_back(f2);
  HPolygons add_faces_insert = hpolygon_insert_edge(f2_borders,i1,i2,
						    pe_p1,pe2,
						    get_e1()->source()->point(),
						    get_p2());
  std::vector<bool> add_faces_insert_match(add_faces_insert.size(),false);
  bool add_face_f1 = false, add_face_f1_twin = false;
  Polygons poly; /* poly face bounded by e1 in the (virtual)
		    tessellation where the line segment p1p2 has been
		    inserted */
  poly = boundaries(add_faces_insert[0]);
  PECirc pe1;
  seg_buf = Segment(get_e1()->source()->point(),get_e1()->target()->point());
  try {
    i1 = find_edge_in_polygons(poly,seg_buf,pe1);
    add_faces_insert_match[0] = true;
  } catch (std::domain_error const& ex0) {
    if (add_faces_insert.size()>1) {
      poly = boundaries(add_faces_insert[1]);
      try {
	i1 = find_edge_in_polygons(poly,seg_buf,pe1);
	add_faces_insert_match[1] = true;
      } catch (std::domain_error const &ex1) {
	// e1 not in add_faces_insert[0] neither in add_faces_insert[1]
	poly = boundaries(face2poly(get_e1()->face(),false));
	add_face_f1 = true;
	try {
	  i1 = find_edge_in_polygons(poly,seg_buf,pe1);
	} catch (std::domain_error const &ex2) {
	  throw std::domain_error("Flip::modified_elements, edge e1 was not "
				  "found anywhere");
	}
      }
    } else { // only one polygon in add_faces_insert
      poly = boundaries(face2poly(get_e1()->face(),false));
      add_face_f1 = true;
      try {
	i1 = find_edge_in_polygons(poly,seg_buf,pe1);
      } catch (std::domain_error const &ex2) {
	throw std::domain_error("Flip::modified_elements, edge e1 was not "
				"found anywhere");
      }
    }
  }
  Polygons poly_twin; /* face bounded by the twin of e1 in the virtual
			 tessellation where the line segment p1p2 has
			 been inserted. */
  seg_buf = seg_buf.opposite();
  Size i1_twin;
  PECirc pe1_twin;
  HPolygon add_face_remove;
  try {
    i1_twin = find_edge_in_polygons(poly,seg_buf,pe1_twin);
    add_face_remove = hpolygon_remove_edge(poly,poly,i1,i1,pe1,pe1_twin,true);
  } catch (std::domain_error const& ex0) {
    /* twin of e1 not in poly */
    if (!add_face_f1) {
      poly_twin = boundaries(face2poly(get_e1()->twin()->face()));
      i1_twin = find_edge_in_polygons(poly_twin,seg_buf,pe1_twin);
      add_face_f1_twin = true;
    } else {
      /* poly_twin in add_faces_insert */
      try {
	poly_twin = boundaries(add_faces_insert[0]);
	i1_twin = find_edge_in_polygons(poly_twin,seg_buf,pe1_twin);
	add_faces_insert_match[0] = true;
      } catch (std::domain_error const& ex1) {
	if (add_faces_insert.size()==1)
	  throw std::domain_error("Flip::modified_elements, twin of e1 not in "
				  "poly neither in add_faces_insert which is "
				  "of length 1");
	poly_twin = boundaries(add_faces_insert[1]);
	add_faces_insert_match[1] = true;
	try {
	  i1_twin = find_edge_in_polygons(poly_twin,seg_buf,pe1_twin);
	} catch (std::domain_error const& ex2) {
	  throw std::domain_error("Flip::modified_elements, poly is the "
				  "face bounded by e1 but twin of e1 not "
				  "found in add_faces_insert");
	}
      }
    }
    add_face_remove = hpolygon_remove_edge(poly,poly_twin,i1,i1_twin,
					   pe1,pe1_twin,false);
  }
  modifs.del_faces.push_back(simplify(f2));
  if (add_face_f1)
    modifs.del_faces.push_back(simplify(HPolygon(poly[0],poly.begin()+1,
						 poly.end())));
  if (add_face_f1_twin) // poly_twin != poly and poly!=e1->face()
    modifs.del_faces.push_back(simplify(HPolygon(poly_twin[0],
						 poly_twin.begin()+1,
						 poly_twin.end())));
  for (Size k=0;k<add_faces_insert.size();k++) 
    if (!add_faces_insert_match[k])
      modifs.add_faces.push_back(simplify(add_faces_insert[k]));
  modifs.add_faces.push_back(simplify(add_face_remove));
    
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
      if (f->data()) // only internal faces
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
    for (HPolygons::iterator i_fc=ml.del_faces.begin(); 
	 i_fc!=ml.del_faces.end();i_fc++){
      loc_var -=(features.faces[i])(*i_fc,ttes);
    } 
    for (HPolygons::iterator i_fc=ml.add_faces.begin(); 
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
void  Energy::add_features_faces(double (*ff)(HPolygon,TTessel*) ){
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
/** \brief Test whether a point is in the interior of a holed polygon
 * \param pt : point to be tested.
 * \param poly : polygon with holes.
 * \note Return false if the point lies on the boundary of the polygon.
 */
bool is_inside(Point2 pt,HPolygon &poly) {
  if (CGAL::bounded_side_2(poly.outer_boundary().vertices_begin(),
			   poly.outer_boundary().vertices_end(),
			   pt)!=CGAL::ON_BOUNDED_SIDE) {
    return false;
  }
  for (HPolygon::Hole_const_iterator hi=poly.holes_begin();
       hi!=poly.holes_end();hi++) {
    if (CGAL::bounded_side_2(hi->vertices_begin(),hi->vertices_end(),pt)
	!=CGAL::ON_UNBOUNDED_SIDE) {
      return false;
    }
  }
  return true;
}
/** \brief Test whether a point is in the interior of a union
 * of holed polygons
 * \param pt : point to be tested.
 * \param polys : polygons with holes.
 * \note Return false if the point lies on the boundary of one of 
 * the polygons.
 */
bool is_inside(Point2 pt,HPolygons& polys) {
  for (Size i=0;i!=polys.size();i++) {
    if (is_inside(pt,polys[i])) {
      return true;
    }
  }
  return false;
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
 * \exception std::domain_error S does not hit P.
 * 
 * When the clipping polygon is not convex, its intersection with the segment
 * may consists of several segments. In such a case, the longest clipped segment is returned. If there are several clipped segments with maximal length, an arbitrarily chosen one is returned.
 */
Segment clip_segment_by_polygon(Segment S, Polygon P) 
  throw (std::domain_error const&) 
{
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
  if (nefInter.is_empty()) {
    throw std::domain_error("function clip_segment_by_polygon failed to return"
			    " a segment because the input segment does not hit"
			    " the polygon");
  } 
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
/** \brief Clip a segment by a polygon with holes
 * \param S : segment to be clipped.
 * \param P : clipping polygon.
 * \return clipped segment, that is the longest intersection between the segment and
 * the polygon.
 * \exception std::domain_error S does not hit P.
 * 
 * When the clipping polygon is not convex, its intersection with the segment
 * may consists of several aligned segments and isolated points. In such a case, the longest clipped segment is returned. If there are several clipped segments with maximal length, an arbitrarily chosen one is returned.
 */
Segment clip_segment_by_polygon(Segment S, HPolygon P)
  throw (std::domain_error const&) 
{
  typedef CGAL::Extended_cartesian<NT> EKernel;
  typedef CGAL::Nef_polyhedron_2<EKernel> Nef;
  std::vector<Nef::Point> Svertices, Pvertices;
  for (Size i=0;i<P.outer_boundary().size();++i) 
    Pvertices.push_back(Nef::Point(P.outer_boundary()[i].x(),
				   P.outer_boundary()[i].y()));
  Nef nefP(Pvertices.begin(),Pvertices.end());
  Nef::Explorer ex_nefP = nefP.explorer();
  for (HPolygon::Hole_const_iterator h=P.holes_begin();
       h!=P.holes_end();h++) {
    std::vector<Nef::Point> Hvertices;
    for (Size j=0;j<(h->size());++j) 
      Hvertices.push_back(Nef::Point((*h)[j].x(),(*h)[j].y()));
    Nef PHole(Hvertices.rbegin(),Hvertices.rend());
    nefP = nefP.difference(PHole);
  }
  Svertices.push_back(Nef::Point(S[0].x(),S[0].y()));
  Svertices.push_back(Nef::Point(S[1].x(),S[1].y()));
  Nef nefS(Svertices.begin(),Svertices.end());
  Nef nefInter = nefS.intersection(nefP);
  if (nefInter.is_empty()) {
    throw std::domain_error("function clip_segment_by_polygon failed to return"
			    " a segment because the input segment does not hit"
			    " the (holed) polygon");
  } 
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
/** \brief Clip a segment by a union of holed polygons
 * \param S : segment to be clipped.
 * \param P : clipping polygons.
 * \return clipped segment, that is the longest intersection between the 
 * segment and the polygon.
 * \exception std::domain_error S does not hit P.
 * The intersection between the segment and the polygons may consists
 * of several aligned segments and isolated points. In such a case,
 * the longest clipped segment is returned. If there are several
 * clipped segments with maximal length, an arbitrarily chosen one is
 * returned.
 */
Segment clip_segment_by_polygon(Segment S, HPolygons P) 
  throw (std::domain_error const&) 
{
  Segment res, buf;
  NT len_max = -1;
  for (HPolygons::const_iterator hpi=P.begin();hpi!=P.end();hpi++) {
    try {
      buf = clip_segment_by_polygon(S,*hpi);
    } catch(std::domain_error const& exception) {
      continue;
    }
    Point2 p1(buf.source());
    Point2 p2(buf.target());
    NT len = CGAL::squared_distance(p1,p2);
    if (len>len_max) {
      res = buf;
      len_max = len;
    }
  }
  if (len_max<0)
    throw std::domain_error("function clip_segment_by_polygon failed to return"
			    " a segment because the input segment does not hit"
			    " any of the (holed) polygons");
  return res;
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
  HPolygons domain = t->get_window();
  return (double)is_inside(pt,domain);
}
/** \brief Test whether a point is in the interior of the tessellated domain
 * \return a double resulting from conversion from a Boolean.
 */
double is_point_inside_window(Point2 pt,TTessel* t){
  // repeated code, any way to avoid that without slowing computation too much?
  HPolygons domain = t->get_window();
  return (double)is_inside(pt,domain);
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
  HPolygons dom = t->get_window();
  // Test whether both lie on the same domain side
  for (Size i=0;i!=dom.size();i++) {
    Polygons borders;
    borders.push_back(dom[i].outer_boundary());
    for (HPolygon::Hole_const_iterator hi=dom[i].holes_begin();
	 hi!=dom[i].holes_end();hi++) {
      borders.push_back(*hi);
    }
    for (Size j=0;j!=borders.size();j++) {
      for (Polygon::Edge_const_iterator e=borders[j].edges_begin();
	   e!=borders[j].edges_end();e++) {
	if (e->has_on(p1) && e->has_on(p2)) {
	  return 0.0;
	}
      }
    }
  }
  return 1.0;
}

/** \brief Return 1
 *
 * \param f: a (holed) polygon representing a tessellation face. 
 * \param t: the tessellation to be considered.
 * \ingroup features
 *
 * Silly function that can be used by an Energy object for 
 * specifying the number of tessellation face as a feature.*/
double face_number(HPolygon f, TTessel* t){
  return 1;
}

/** \brief Return the squared face area of a tessellation face
 *
 * \param f : tessellation face as a (holed) polygon.
 * \param t : the tessellation to be considered.
 * \ingroup features*/
double face_area_2(HPolygon f, TTessel* t){
  Polygons p = boundaries(f);
  double area = 0.0; 
  for (Size i=0;i!=p.size();i++) {
    area += CGAL::to_double(p[i].area());
  }
  return area;
}
/** \brief Return the perimeter of a tessellation face
 *
 * \param f : tessellation face as a (holed) polygon.
 * \param t : the tessellation to be considered.
 * \ingroup features*/ 
double face_perimeter(HPolygon f, TTessel* t){
  Polygons p = boundaries(f);
  double l = 0.0; 
  for (Size i=0;i!=p.size();i++) {
    for (Size j=0;j!=p[i].size();j++) {
      l += sqrt(CGAL::to_double(squared_distance(p[i][j],
						 p[i][(j+1)%p[i].size()])));
    }
  }
  return l;
}
/** \brief Return the shape parameter of a tessellation face.
 *
 * \param f : tessellation face as a (holed) polygon
 * \param t : the tessellation to be considered
 *
 * The shape parameter is the ratio between the perimeter
 * and 4 times the squared root of the area. */
double face_shape(HPolygon f, TTessel* t){
  return face_perimeter(f,t)/(4*pow(CGAL::to_double(face_area_2(f,t)),0.25)); 
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
 * \param f : T-tessellation face as a (holed) polygon.
 * \param t : the T-tessellation to be considered.
 * \ingroup features
 * 
 * Only vertices of the outer boundary where incident edges form an
 * acute angle contribute to the measure. The additive contribution of
 * such a vertex \f$v\f$ is
 * \f[
 * \frac{\pi}{2}-\phi(v),
 * \f]
 * where \f$\phi(v)\f$ is the (acute) angle between the edges incident
 * to \f$v\f$. The above contribution is close to zero if the incident
 * edges are almost perpendicular. It is close to \f$\pi/2\f$ if the 
 * edges are almost aligned. 
*/
double face_sum_of_angles(HPolygon f, TTessel* t){
  Vector v1,v2;
  Point2 p;
  double tot_angle=0;

  Polygon ob = f.outer_boundary();
  for(int i=0;i!=ob.size();i++) {
    int j=(i-1)%((int)ob.size());
    if (j<0) j=j+ob.size();
    p=ob[i];
    //if (is_point_inside_window(p,t)>.5)
      {
	v1 = Vector(ob[i],ob[j]);
        v2 = Vector(ob[i],ob[(i+1)%((int)ob.size())]);
	double a = CGAL_PI/2-angle_between_vectors(v1,v2);
	if (a>=0)
	  tot_angle += a;
      }
}
  return tot_angle;
}

/** \brief Returns a HPolygon object representing a face of the
 *     tessellation. 
 * \param f : the tessellation face to be converted.
 * \param simplify : if true, only vertices joining non-aligned edges 
 * are specified during the HPolygon construction. Default: true.
 * \return a holed polygon.
 **/
HPolygon face2poly(TTessel::Face_handle f, bool simplify) {
  std::vector<TTessel::Ccb_halfedge_circulator> ccbs;
  ccbs.push_back(f->outer_ccb());
  if (std::distance(f->holes_begin(),f->holes_end())>0) {
    ccbs.insert(ccbs.end(),f->holes_begin(),f->holes_end());
  }
  Polygons poly(ccbs.size());
  for (Size i=0;i<ccbs.size();i++) {
    TTessel::Ccb_halfedge_circulator e = ccbs[i];
    TTessel::Halfedge_handle e_first=e;
    do {
      if (!simplify || e->get_next_hf()==NULL_HALFEDGE_HANDLE || 
	  (e->get_next_hf()!=e->next())) 
	poly[i].push_back(e->target()->point());
      e++;
    } while (e!=e_first);
    
  }
  Polygons::iterator hi(poly.begin());
  hi++;
  HPolygon res(poly[0],hi,poly.end());
  return(res);
}
/** \brief Return the smallest angle in a tessellation face
 *
 * \param f : tessellation face as a (holed) polygon.
 * \param t : the tessellation to be considered.
 * \ingroup features
 *
 * Angles between consecutive edges of the outer boundary are
 * considered. This function can be used by an Energy object for
 * specifying the sum of face smallest angle as a feature. That
 * feature is of interest for penalizing faces with small angles.
*/
double min_angle(HPolygon f, TTessel* t){
  Vector v1,v2;
  Point2 p;
  double min_angle;
 
  Polygon ob = f.outer_boundary();
  min_angle = angle_between_vectors(Vector(ob[ob.size()-1],ob[0]),
				    Vector(ob[0],ob[1]));
  for(int i=0;i!=ob.size();i++) {
    int j=(i-1)%((int)ob.size());
    if (j<0) j=j+ob.size();
    p=ob[i];
    v1 = Vector(ob[i],ob[j]);
    v2 = Vector(ob[i],ob[(i+1)%((int)ob.size())]);
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
  HPolygon poly;
  for (TTessel::Face_iterator f = t->faces_begin();f!=t->faces_end();f++) {
    if (f->data()) { 
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
  HPolygon poly;
  for (TTessel::Face_iterator f = t->faces_begin();f!=t->faces_end();f++) {
    if (f->data()) { 
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
  HPolygon poly;
  for (TTessel::Face_iterator f = t->faces_begin();f!=t->faces_end();f++) {
    if (f->data()) { 
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

/** \brief Return the connected components of a holed polygon boundary
 * 
 * First polygon in the output is the outer boundary oriented 
 * counterclockwise. Subsequent polygons are holes boundaries oriented
 * clockwise.
 */
Polygons boundaries(HPolygon hp) {
  Polygons b;
  Polygon ob(hp.outer_boundary());
  if (ob.area()<0) {
    ob.reverse_orientation();
  }
  b.push_back(ob);
  for (HPolygon::Hole_const_iterator h=hp.holes_begin();h!=hp.holes_end();h++) {
    
    Polygon poly(*h);
    if (poly.area()>0) {
      poly.reverse_orientation();
    }
    b.push_back(poly);
  }
  return b;
}
/** \brief Remove unnecessary vertices on a polygon
 *
 * A vertex is not necessary if it is collinear with its neighbours. 
 */
Polygon simplify(Polygon p) throw(std::domain_error const&) {
  if (p.is_empty()) {
    throw std::domain_error("simplify(Polygon), input polygon is empty");
  }
  Polygon sp;
  
  /* Start by removing all unnecessary vertices */
  Polygon::Vertex_circulator pv2 = p.vertices_circulator();
  pv2--;
  Polygon::Vertex_circulator pv1 = pv2++, pv0 = pv1++, pv_end = pv1;
  pv2++;
  do {
    if(!CGAL::collinear(*pv0,*pv1,*pv2)) {
      sp.push_back(*pv1);
    }
    pv0++; pv1++; pv2++;
  } while (pv1!=pv_end);
  return sp;
}
/** \brief Remove unnecessary vertices on a holed polygon
 */
HPolygon simplify(HPolygon p) {
  Polygons b = boundaries(p);
  Polygons::iterator hb = b.begin();
  hb++;
  Polygons::iterator he = b.end();
  b[0] = simplify(b[0]);
  for (Polygons::iterator h=hb;h!=he;h++) {
    *h = simplify(*h);
  }
  HPolygon sp(b[0],hb,he);
  return sp;
}
/** \brief Compute the intersection of a ray with a face boundary
 * \param r : ray starting inside the face.
 * \param f : LineTes face.
 * \param e : on exit, the edge (on the face boundary) where the ray
 *            exits the face for the first time.
 * \return : the location where the ray exits the face.
 * \pre : the ray source must be inside the face or on its boundary.
 */ 
Point2 ray_exit_face(Rayon &r,LineTes::Face_handle &f,
		     LineTes::Halfedge_handle &e) {
  /* The target halfedge may be on the outer ccb or on a "hole"
     boundary. Approach below: visit all halfedges (on outer ccb or
     along holes), compute intersection, if any keep the closest
     one to the ray source. */
  Point2 res;
  std::vector<LineTes::Halfedge_handle> all_edges;
  // push all halfedges of the outer ccb to vector all_edges
  LineTes::Ccb_halfedge_circulator hc=f->outer_ccb();
  LineTes::Halfedge_handle h0 = hc;
  do {
    all_edges.push_back(hc);
    hc++;
  } while (hc!=h0);
  // push all halfedges on inner boundaries to vector all_edges
  for (LineTes::Hole_iterator ho=f->holes_begin();ho!=f->holes_end();ho++) {
    hc = *ho;
    h0 = hc;
    do {
      all_edges.push_back(hc);
      hc++;
    } while (hc!=h0);
  }
  // find the intersections of the ray and keep the closest one
  // to the ray source
  double lp = std::numeric_limits<double>::max();
  Point2 inter_location;
  bool found(false);
  for (std::vector<LineTes::Halfedge_handle>::iterator ce=all_edges.begin();
       ce!=all_edges.end();ce++) {
    CGAL::Object inter = CGAL::intersection(r,
					    Segment((*ce)->source()->point(),
						    (*ce)->target()->point())); 
    if (CGAL::assign(inter_location,inter)) {
      NT len2 = CGAL::squared_distance(r.source(),inter_location);
      double len = sqrt(CGAL::to_double(len2));
      if (len==0) continue;
      /* It may happen that both an halfedge and its twin bound the face. */
      if (found && *ce==e->twin()) { 
	/* Select the halfedge that has
	   the source of the ray on its left */
	Line ce_line((*ce)->source()->point(),(*ce)->target()->point());
	if (ce_line.oriented_side(r.source())==CGAL::ON_POSITIVE_SIDE) {
	  e = *ce;
	} // otherwise retain the previous minimal halfedge e (twin of ce)
	continue;
      }
      if (len<lp) { // ce is closest than e. Keep it. 
	res = inter_location;
	if (!found) found = true;
	e = *ce;
	lp = len;
      }
    }
  }
  if (!found) {
    std::cerr << "ray_exit_face: intersection of ray with face ";
    std::cerr << "boundaries not found" << std::endl;
  }
  return res;
}
/** Convert a connected component boundary circulator into a polygon edge
 * circulator.
 * \param e : a CCB halfedge circulator to be converted.
 * \param poly : the CCB as a polygon.
 * \return The alter ego of e as a polygon edge circulator. It starts at
 * the same edge than e.
 * \sa ccb2polygon.  
 */
PECirc 
lt2poly_edge_circulator(LineTes::Halfedge_handle &e, Polygon &poly) {
  PECirc e_poly = poly.edges_circulator(), 
    done = e_poly;
  Segment seg(e->source()->point(),e->target()->point());
  do {
    if (*e_poly==seg) {
      return e_poly;
    } else {
      e_poly++;
    }
  } while (e_poly!=done);
  throw std::domain_error("function lt2poly_edge_circulator failed to find"
			  " the given edge in the given polygon");
}
/** \brief Convert a connected component boundary into a polygon
 * \param e : a circulator for the CCB to converted.
 * \return the polygon defined by the CCB.
 */
Polygon ccb2polygon(LineTes::Ccb_halfedge_circulator e) {
  Polygon res;
  LineTes::Ccb_halfedge_circulator done = e;
  do {
    res.push_back(e->source()->point());
    e++;
  } while (e!=done);
  return res;
}
/** \brief Compute a polygon generated when joining two polygon edges by a 
 *  line segment.
 *
 * \param e1 : a polygon edge circulator starting at one of the edges where 
 *             the joining line segment starts.
 * \param e2 : a polygon edge circulator starting at the other edge where 
 *             the joining line segment ends.
 * \param p1 : the point on e1 where the joining line segment
 *             starts.
 * \param p2 : the point on e2 where the joining line segment
 *             ends.
 * \return the generated polygon containing the edge (p1,p2).
 *
 * The edges e1 and e2 may belong to the same polygon. Then, joining e1 and e2 
 * splits the polygon into two polygons. The returned polygon is the one 
 * containing the edge (p1,p2). The other generated polygon can be obtained 
 * by calling the function swapping e1 and e2, p1 and p2.
 *
 * The edges e1 and e2 may also belong to different polygons. Then, joining e1
 * and e2 merges their containing polygons. The returned polygon is the merged
 * polygon. Note that it is not simple as its boundary self-intersects along the
 * line segment (p1,p2).
 */
Polygon polygon_insert_edge(PECirc e1, 
			    PECirc e2, 
			    Point2 p1, Point2 p2) {
  std::vector<Point2> res;
  res.push_back(p1);
  res.push_back(p2);
  PECirc e = e2;
  if (e->target()==p2)
    e++;
  do {
    res.push_back(e->target());
    e++;
  } while (e!=e2 && e!=e1);
  if (e==e1) {
    if (e1->source()==p1)
      res.pop_back();
    return Polygon(res.begin(),res.end());
  }
  // e==e2
  res.push_back(p2);
  res.push_back(p1);
  e = e1;
  if (e1->target()==p1) 
    e++;
  do {
    res.push_back(e->target());
    e++;
  } while (e!=e1);
  if (e1->source()==p1)
    res.pop_back();
  return Polygon(res.begin(),res.end());
}
/** \brief Compute the holed polygon(s) generated when joining two polygon 
 * edges by a line segment.
 *
 * \param hpoly : a holed polygon as a vector of polygons. The first polygon 
 * is the counterclockwise oriented outer boundary. The following polygons are 
 * clockwise oriented inner boundaries bounding holes. Such a representation of
 * a holed polygon can be obtained by using the function boundaries.
 * \param i1 : the index of the border where the inserted line segment starts. 
 * \param i2 : the index of the border where the inserted line segment ends. 
 * \param e1 : a polygon edge circulator starting at the edge where 
 *             the inserted line segment starts.
 * \param e2 : a polygon edge circulator starting at the other edge where 
 *             the inserted line segment ends.
 * \param p1 : the point on e1 where the inserted line segment
 *             starts.
 * \param p2 : the point on e2 where the inserted line segment
 *             ends.
 * \return the holed polygon(s) generated by the insertion of the line segment.
 *
 * The insertion of the line segment can result in a 
 * split  or in a modification of the input holed polygon. The former case 
 * occurs when the inserted linne segment joins two points on the outer 
 * boundary. Then, the generated holed 
 * polygons are returned as a vector of of length 2. In the latter
 * case, the modified holed polygon is returned as a vector of length 1.
 * 
 */
HPolygons hpolygon_insert_edge(Polygons &hpoly, Size i1, Size i2, PECirc e1, 
			       PECirc e2, Point2 p1, Point2 p2) {
  PECirc ob_circ = hpoly[0].edges_circulator();
  HPolygons res;
  int cas = 0;
  /* 4 possible cases
     case 1: the splitting line segment connects two points along 
             the outer boundary.
     case 2: the splitting line segment connects two points along
             the same inner boundary (bounding a hole).
     case 3: the splitting line segment connects a point on the 
             outer boundary to a point on an inner boundary.
     case 4: the splitting line segment connects two points on two
             different inner boundaries.
  */
  // Identify case
  if (hpoly.size()==1) {
    cas = 1;
  } else {
    if (i1==i2) { // e1 and e2 on same ccb
      if (i1==0) cas = 1; // e1 and e2 on outer boundary
      else cas = 2; // e1 and e2 on same inner boundary
    } else { // e1 and e2 on different ccb's
      if (i1==0) { // e1 on outer boundary and e2 on inner boundary
	cas = 3;
      } else { // e1 on inner boundary 
	/* e2 may be on another inner boundary or on the outer
	   boundary */
	if (i2==0) { // e1 on inner boundary, e2 on outer boundary
	  cas = 3;
	} else cas = 4; // e1 and e2 on different inner boundaries
      }
    }
  }
  // Predict generated holed polygon(s)
  Polygon outer1, outer2, ccb_buf;
  Polygons holes1, holes2;
  HPolygon hpoly1, hpoly2;
  std::vector<bool> inside;
  Size j;
  switch(cas) {
  case 1 : {
    outer1 = polygon_insert_edge(e1,e2,p1,p2);
    outer2 = polygon_insert_edge(e2,e1,p2,p1);
    if (hpoly.size()==1) {
      hpoly1 = HPolygon(outer1);
      hpoly2 = HPolygon(outer2);
    } else {
      inside = filter_holes(hpoly.begin()+1,
			    hpoly.end(),outer1);
      for (Size ii=1;ii<hpoly.size();ii++) {
	if (inside[ii-1]) {
	  holes1.push_back(hpoly[ii]);
	}
	else
	  holes2.push_back(hpoly[ii]);
      }
      hpoly1 = HPolygon(outer1,holes1.begin(),holes1.end());
      hpoly2 = HPolygon(outer2,holes2.begin(),holes2.end());
    }
    res.push_back(hpoly1);
    res.push_back(hpoly2);
  } break;
  case 2 : {
    /* the split generates two new faces referenced as mother and
       daughter faces. The daughter face lies along the hole (with 
       inner boundary containing p1 and p2). The mother face is the
       past face from which the daughter is removed. */
    // Outer boundary of the mother face
    outer1 = hpoly[0];
    // Outer boundary of the daughter face
    // new ccb containing (p2,p1)
    Polygon c21 = polygon_insert_edge(e2,e1,p2,p1);
    Polygon enlarged_hole;
    if (CGAL::bounded_side_2(c21.vertices_begin(),c21.vertices_end(),
			     e2->target())!=CGAL::ON_UNBOUNDED_SIDE) {
      // the outer boundary of the daughter face is c12
      outer2 = polygon_insert_edge(e1,e2,p1,p2);
      enlarged_hole = c21;
    } else { 
      /* if target of e2 is outside c21, c21 is the outer boundary of
	 the daughter face */
      outer2 = c21;
      enlarged_hole = polygon_insert_edge(e1,e2,p1,p2);
    }
    // Dispatch holes
    inside = filter_holes(hpoly.begin()+1,hpoly.end(),outer2);
    for (Size ii=1;ii<hpoly.size();ii++) {
      if (ii==i1) { // hole lying along the daughter face
	holes1.push_back(enlarged_hole);
      } else if (inside[ii-1]) {
	holes2.push_back(hpoly[ii]);
      }
      else {
	holes1.push_back(hpoly[ii]);
      }
    }
    hpoly1 = HPolygon(outer1,holes1.begin(),holes1.end());
    hpoly2 = HPolygon(outer2,holes2.begin(),holes2.end());
    res.push_back(hpoly1);
    res.push_back(hpoly2);
  } break;
  case 3 : {
    /* The existing face is modified. Its outer boundary folds 
       inside the face and connects to a hole. */
    // ccb containing e2
    outer1 = polygon_insert_edge(e1,e2,p1,p2);
    if (i1>0) {
      j = i1;
    } else {
      j = i2;
    }
    for (Size ii=1;ii<hpoly.size();ii++) {
      if (ii!=j) {
	holes1.push_back(hpoly[ii]);
      }
    }
    hpoly1 = HPolygon(outer1,holes1.begin(),holes1.end());
    res.push_back(hpoly1);
  } break;
  case 4 : {
    /* Two holes are merged */
    outer1 = hpoly[0];
    for (Size ii=1;ii<hpoly.size();ii++) {
      if (ii!=i1 && ii!=i2) {
	holes1.push_back(hpoly[ii]);
      }
    }
    // merge the two holes and add the result to the list of holes
    ccb_buf = polygon_insert_edge(e1,e2,p1,p2);
    holes1.push_back(ccb_buf);
    hpoly1 = HPolygon(outer1,holes1.begin(),holes1.end());
    res.push_back(hpoly1);
  }
  }
  return res;
}
/** \brief Compute the holed polygon generated by when removing an edge
 * \param hpoly : a holed polygon (represented as a vector of polygons) 
 *                bounded by the oriented edge e to be removed.
 * \param hpoly_twin : a holed polygon bounded by the e_twin the edge
 *                     opposite to e.
 * \param i : the index of the border of hpoly where e lies.
 * \param i_twin : the index of the border of hpoly_twin where e_twin lies.
 * \param e : a polgon edge circulator pointing to the edge to be removed.
 * \param e_twin : a polygon edge circulator pointing to the edge opposite
 *                 to e.
 * \param single polygon : true if hpoly and hpoly_twin are the same holed
 *                         polyhon. False if they are distinct.
 * \return the holed poygon generated by the removal of e (and e_twin). 
 *         When hpoly and hpoly_twin are distinct, they are merged when
 *         shared edge is removed. When hpoly and hpoly_twin coincide,
 *         the initial holed polygon is modified by the removal of e.
 * \note When hpoly and hpoly_twin represent the same holed polygon, it is
 * important to pass e and e_twin as circulators computed on the same object.
 * In particular, e_twin should not be computed on a copy of hpoly.
 */
HPolygon hpolygon_remove_edge(Polygons &hpoly, Polygons &hpoly_twin, Size i, 
			       Size i_twin, PECirc e, PECirc e_twin, 
			       bool single_polygon) {
  int cas = 0;
  /* 4 possible cases
     case 1: the removed edge is separating two tessellation faces. 
     case 2: the removed edge is joining two vertices on the same inner
             boundary of a tessellation face.
     case 3: the removed edge is an isthmus where the outer boundary of a
             tessellation face folds in the face interior and connects to
	     a "hole".
     case 4: the removed edge is an isthmus connecting two "holes". */
  // Identify case
  bool hpoly_equals_hpoly_twin = single_polygon;
  if (!hpoly_equals_hpoly_twin && hpoly.size()==1 && 
      hpoly_twin.size()==1) { 
    cas = 1;
  } else {
    if (hpoly_equals_hpoly_twin && i==i_twin) {
      if (i==0) {
	cas = 3;
      } else {
	cas = 4;
      }
    } else { 
      // hpoly!=hpoly_twin or i!=i_twin i.e. e and its twin on different ccb's
      if (i>0) { // e along an inner boundary
	cas = 2; // only possible case: hpoly!=hpoly_twin i>0 and i_twin==0
      } else { // hpoly!=hpoly_twin or i!=i_twin and i==0
	if (i_twin==0) { // hpoly must be different from hpoly_twin
	  cas = 1;
	} else { // i==0 and i_twin>0, hpoly must be different from hpoly_twin
	  cas = 2;
	}
      }
    }
  }
  // Added faces
  Polygon outer_new, inner_new, buf_poly;
  Polygons holes_new, daughter, mother;
  HPolygon hpoly_new;
  Size j;
  switch(cas) {
  case 1:
    outer_new = polygon_remove_edge(e,e_twin);
    if (hpoly.size()>1) 
      holes_new = Polygons(hpoly.begin()+1,hpoly.end());
    if (hpoly_twin.size()>1)
      holes_new.insert(holes_new.end(),hpoly_twin.begin()+1,hpoly_twin.end());
    break;
  case 2:
    if (i>0) { // hpoly is the mother face
      mother = hpoly;
      daughter = hpoly_twin;
    } else { // hpoly is the daughter face
      mother = hpoly_twin;
      daughter = hpoly;
    }
    outer_new = mother[0];
    if (i>0) {
      j = i;
    } else {
      j = i_twin;
    }
    // push the daughter holes to the list of holes
    if (daughter.size()>1) 
      holes_new = Polygons(daughter.begin()+1,daughter.end());
    // now add the mother holes
    for (Size ii=1;ii<mother.size();ii++) {
      if (ii!=j) {
	holes_new.push_back(mother[ii]);
      } else { // mother hole reduced by merging
	if (i>0) {
	  buf_poly = polygon_remove_edge(e,e_twin);
	} else {
	  buf_poly = polygon_remove_edge(e_twin,e);
	}
	holes_new.push_back(buf_poly);
      }
    }
    break;
  case 3:
    outer_new = polygon_remove_edge(e,e_twin);
    inner_new = polygon_remove_edge(e_twin,e);
    if (outer_new.area()<=0) { // exchange inner_new and outer_new
      buf_poly = outer_new;
      outer_new = inner_new;
      inner_new = buf_poly;
    }
    if (hpoly.size()>1)
      holes_new = Polygons(hpoly.begin()+1,hpoly.end());
    holes_new.push_back(inner_new);
    break;
  case 4:
    outer_new = hpoly[0];
    // Holes : all those of the containing face except the one bounded
    // by e
    if (hpoly.size()>1)
      holes_new = Polygons(hpoly.begin()+1,hpoly.end());
    j = i;
    holes_new.erase(holes_new.begin()+j-1);
    // Merge -> the hole bounded by e divides into two holes
    buf_poly = polygon_remove_edge(e,e_twin);
    holes_new.push_back(buf_poly);
    buf_poly = polygon_remove_edge(e_twin,e);
    holes_new.push_back(buf_poly);
  }
  hpoly_new = HPolygon(outer_new,holes_new.begin(),holes_new.end());
  return hpoly_new;
}
/** \brief Compute polygon(s) generated by the removal of a polygon edge.
 * \param e1 : a polygon edge circulator starting at the edge to be removed.
 * \param e2 : a polygon edge circulator starting at an edge opposite to e1.
 * \return the generated polygon that starts at the end of the removed edge.
 * \pre If e1 and e2 belong to different polygons, both polygons must have the
 * same orientation.
 *
 * There two cases considered. In the first case, the edge to be
 * removed is shared by two adjacent polygons. Removal of that edge
 * merges the initial polygons. The second case involves non simple
 * polygon with a boundary that "visits" the same edge twice. When
 * that edge is removed, the initial polygon is split in two polygons. 
 * In both cases, e1 and e2 are opposite (their starting and ending points
 * are swapped) to each other.
 *
 * In the first case, the merged polygon is returned. In the second case, the 
 * generated polygon "in front of e1" is returned. In order to get the other 
 * generated polygon behind e1, call the function swapping e1 and e2.
*/
Polygon polygon_remove_edge(PECirc &e1,
			    PECirc &e2) {
  Polygon res;
  PECirc circ(e1);
  PECirc done_merge = e1;
  PECirc done_break = e2;
  circ++;
  while (circ!=done_merge && circ!=done_break) {
    res.push_back(circ->target());
    circ++;
  }
  if (circ==done_break)
    return res;
  circ = e2;
  done_merge = circ;
  circ++;
  while (circ!=done_merge) {
   res.push_back(circ->target());
    circ++;
  }
  return res;
}
/** \brief Test whether some halfedges lie on the same connected
 * component boundary as a given halfedge.
 *
 * \param e0 : the reference halfedge.
 * \param es : the halfedges to be tested.
 * \return   : a vector of Boolean, true if the corresponding halfedge
 *             is in the same connected component boundary as e0.
 */
std::vector<bool> is_on_same_ccb(LineTes::Halfedge_handle &e0,
				 std::vector<LineTes::Halfedge_handle> &es) {
  std::vector<bool> res(es.size(),false);
  LineTes::Ccb_halfedge_circulator e = e0->ccb();
  do {
    for (Size i=0;i<es.size();i++) {
      if (es[i]==e) {
	res[i] = true;
      }
    }
    e++;
  } while (e!=e0);
  return res;
}
/** \brief Test whether polygon edges lie on the same connected
 * component boundary as a given edge.
 *
 * \param e0 : the reference edge.
 * \param es : the edges to be tested.
 * \return   : a vector of Boolean, true if the corresponding edge
 *             is in the same connected component boundary as e0.
 */
std::vector<bool> is_on_same_ccb(PECirc e0, std::vector<PECirc> es) {
  std::vector<bool> res(es.size(),false);
  PECirc e = e0;
  do {
    for (Size i=0;i<es.size();i++) {
      if (es[i]==e) {
	res[i] = true;
      }
    }
    e++;
  } while (e!=e0);
  return res;
}
/** \brief Test whether two halfedges lie on the same connected
 * component boundary.
 *
 * \param e0 : a halfedge.
 * \param e1 : another halfedge.
 */
bool is_on_same_ccb(LineTes::Halfedge_handle &e0,
		    LineTes::Halfedge_handle &e1) {
  LineTes::Ccb_halfedge_circulator e = e0->ccb();
  do {
    if (e1==e) {
      return true;
    }
    e++;
  } while (e!=e0);
  return false;
}
/** \brief Test whether edges lie on the same connected
 * component boundary of a holed polygon.
 *
 * \param e0 : an edge.
 * \param e1 : another edge.
 */
bool is_on_same_ccb(PECirc e0, PECirc e1) {
  PECirc e = e0;
  do {
    if (e1==e) {
      return true;
    }
    e++;
  } while (e!=e0);
  return false;
}
/** \brief Test whether a face has one or more holes
 * \param f : the face to be considered.
 * \return true if there is a least one hole in f, false if
 * f has no no hole.
 */
bool has_holes(LineTes::Face_handle &f) {
  return std::distance(f->holes_begin(),f->holes_end())>0;
}
/** \brief Test whether a holed polygon has one or more holes
 * \param p : the holed polygon to be considered.
 * \return true if there is a least one hole in p, false if
 * p has no no hole.
 */
bool has_holes(HPolygon &p) {
  return std::distance(p.holes_begin(),p.holes_end())>0;
}
/** \brief Test whether the first points of holes boundaries are 
 * inside a reference polygon.
 *
 * \param begin : iterator to the first hole to be tested.
 * \param end : iterator past the last hole to be tested.
 * \param ref : reference polygon.
 * \return : a vector of Boolean, true for each tested hole
 *           with its first point contained in the reference 
 *           polygon.
 * \note{Can be used a an inclusion test when the tested holes
 *       are assumed to be either inside or outside 
 *       the reference polygon.}
 */
std::vector<bool> filter_holes(HPolygon::Hole_const_iterator begin,
			       HPolygon::Hole_const_iterator end,
			       Polygon &poly) {
  std::vector<bool> res;
  for (HPolygon::Hole_const_iterator hi=begin;hi!=end;hi++) {
    //if (poly.bounded_side((*hi)[0])==CGAL::ON_BOUNDED_SIDE)
    if (CGAL::bounded_side_2(poly.vertices_begin(),poly.vertices_end(),(*hi)[0])
	==CGAL::ON_BOUNDED_SIDE)
      res.push_back(true);
    else
      res.push_back(false);
  }
  return res;
}
/** \brief Test whether the first points of polygons are 
 * inside a reference polygon.
 *
 * \param begin : iterator to the first polygon to be tested.
 * \param end : iterator past the last polygon to be tested.
 * \param ref : reference polygon.
 * \return : a vector of Boolean, true for each tested polygon
 *           with its first point contained in the reference 
 *           polygon.
 * \note{Can be used a an inclusion test when the tested polygons
 *       are assumed to be either inside or outside 
 *       the reference polygon.}
 */
 std::vector<bool> filter_holes(Polygons::const_iterator begin,
			       Polygons::const_iterator end,
			       Polygon &poly) {
  std::vector<bool> res;
  for (Polygons::const_iterator hi=begin;hi!=end;hi++) {
    if (CGAL::bounded_side_2(poly.vertices_begin(),poly.vertices_end(),(*hi)[0])
	==CGAL::ON_BOUNDED_SIDE)
      res.push_back(true);
    else
      res.push_back(false);
  }
  return res;
}
/** \brief Compute the index of the hole containing a given halfedge in
 * its boundary
 * \param e : halfedge.
 * \param begin : first hole boundary to be searched.
 * \param end : last hole boundary to be searched.
 * \return : the position of the hole whose boundary contains e. For
 * instance, if 0, the first hole was found to countain e in its boundary.
 * \exception std::domain_error if e was not found in any hole.
 */
Size hole_index(LineTes::Halfedge_handle e,LineTes::Hole_iterator begin,
		LineTes::Hole_iterator end) throw(std::domain_error const&)
{
  std::vector<TTessel::Halfedge_handle> vec_holes(begin,end);
  std::vector<bool> test_holes = is_on_same_ccb(e,vec_holes);
  std::vector<bool>::iterator hole_pos = std::find(test_holes.begin(),
						   test_holes.end(),
						   true);
  if (hole_pos==test_holes.end()) {
    std::ostringstream msg;
    msg << "hole_index failed to find halfedge "
	<< "(" << e->source()->point().x() << ","
	<< e->source()->point().y() << ")-("
	<< e->target()->point().x() << "," << e->target()->point().y()
	<< ") in any of the provided holes.";
    throw std::domain_error(msg.str());
  }
  Size j = std::distance(test_holes.begin(),hole_pos);
  return j;
}
/** \brief Find among polygons the polygon that contains a given edge
 * \param e : edge.
 * \param begin : first polygon to be searched.
 * \param end : last polygon to be searched.
 * \return : the position of the polygon whose boundary contains e. For
 * instance, if 0, the first polygon was found to countain e in its boundary.
 * \exception std::domain_error if e was not found in any polygon.
 */
Size polygon_index(PECirc e,Polygons::const_iterator begin,
		Polygons::const_iterator end) throw(std::domain_error const&)
{
  std::vector<PECirc> vec_edges;
  for (Polygons::const_iterator p=begin;p!=end;p++) {
    vec_edges.push_back(p->edges_circulator());
  }
  std::vector<bool> test_polys = is_on_same_ccb(e,vec_edges);
  std::vector<bool>::iterator poly_pos = std::find(test_polys.begin(),
						   test_polys.end(),
						   true);
  if (poly_pos==test_polys.end()) {
    std::ostringstream msg;
    msg << "poly_index failed to find edge "
	<< "(" << e->source().x() << ","
	<< e->source().y() << ")-("
	<< e->target().x() << "," << e->target().y()
	<< ") in any of the provided holes.";
    throw std::domain_error(msg.str());
  }
  Size j = std::distance(test_polys.begin(),poly_pos);
  return j;
}
/** \brief Search a segment among the edges of a series of polygons
 * \param polys : the series of polygons to be searched.
 * \param seg : the segment to be found.
 * \param e : on exit, an edge circulator starting at the edge equal to
 *            the given segment.
 * \return the index of the polygon where the segment was found.
 * \note if s is part of the boundaries of several polygons, the search 
 * stops at the first match.
 */
Size find_edge_in_polygons(Polygons &polys, Segment seg, PECirc &e) 
  throw (std::domain_error const&) {
  bool found = false;
  Size index;
  for (Size i=0;i<polys.size();i++) {
    e = polys[i].edges_circulator();
    PECirc done = e;
    do {
      if (e->source()==seg.source() && e->target()==seg.target()) {
	found = true;
	index = i; 
	break;
      } else {
	e++;
      }
    } while (e!=done);
    if (found) break;
  }
  if (!found) {
    std::ostringstream msg;
    msg << "find_edge_in_polygons: the segment ("
	<< seg.source().x() << "," << seg.source().y() << ")-("
	<< seg.target().x() << "," << seg.target().y() << ") is "
	<< "not an edge of the given polygons";
    throw std::domain_error(msg.str());
  } else {
    return index;
  }
}
/** \brief Search a segment among the edges of holed polygon
 * \param poly : the holed polygon to be searched.
 * \param seg : the segment to be found.
 * \param e : on exit, an edge circulator starting at the edge equal to
 *            the given segment.
 * \return the index of the border where the segment was found. Zero if the
 * segment is on the outer boundary. A positive value if the segment is on
 * an inner boundary.
 */
Size find_edge_in_hpolygon(HPolygon &poly, Segment seg, PECirc &e) {
  Polygons borders = boundaries(poly);
  return find_edge_in_polygons(borders, seg, e);
}
