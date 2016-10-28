#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SMFPredictions 
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "check_lite.h"

Point2 upper_left_vertex(const Polygon &poly) {
  return *(std::min_element(poly.vertices_begin(),poly.vertices_end()));
}
bool llt(const Polygon &a,const Polygon &b) {
  Point2 ul_a(upper_left_vertex(a)), ul_b(upper_left_vertex(b));
  return ul_a<ul_b;
}
bool operator==(const HPolygon &a,const HPolygon &b) {
  if(a.outer_boundary()!=b.outer_boundary()) {
    return false;
  }
  Polygons holes_a(a.holes_begin(),a.holes_end()),
    holes_b(b.holes_begin(),b.holes_end());
  if (holes_a.size()!=holes_b.size()) return false;
  std::sort(holes_a.begin(),holes_a.end(),llt);
  std::sort(holes_b.begin(),holes_b.end(),llt);
  bool buf = std::mismatch(holes_a.begin(),holes_a.end(),
		       holes_b.begin()).first==holes_a.end();
  return std::mismatch(holes_a.begin(),holes_a.end(),
		       holes_b.begin()).first==holes_a.end();
}
HPolygons get_faces(TTessel &tes) {
  HPolygons faces;
  for (TTessel::Face_iterator fi=tes.faces_begin();fi!=tes.faces_end();fi++) {
    if (fi->data()) {
      HPolygon theface;
      theface = face2poly(fi);
      faces.push_back(theface);
    }
  }
  return faces;
}
HPolygons set_diff(HPolygons &a,HPolygons &b) {
  HPolygons res;
  for (HPolygons::const_iterator p=a.begin();p!=a.end();p++) {
    HPolygons::const_iterator found;
    found = std::find(b.begin(),b.end(),*p);
    if (found==b.end()) {
      res.push_back(*p);
    }
  }
  return res;
}
std::ostream& operator<<(std::ostream &os,const HPolygon &p) {
  Polygons bds = boundaries(p);
  for (Polygons::const_iterator bi=bds.begin();bi!=bds.end();bi++) {
    for (unsigned j=0;j<bi->size();j++) {
      os << CGAL::to_double(bi->vertex(j).x());
      os << "," << CGAL::to_double(bi->vertex(j).y());
      os << std::endl;
    }
    os << "NA,NA" << std::endl;
  }
  return os;
}
std::ostream& operator<<(std::ostream &os,const HPolygons &ps) {
  for (HPolygons::const_iterator pi=ps.begin();pi!=ps.end();pi++) {
    HPolygon poly(*pi);
    os << poly;
  }
  return os;
}
std::istream& operator>>(std::istream &is, HPolygons &ps) {
  /* Read holed polygons from a CVS file */
  char sep = ',';
  std::string line;
  Polygon poly, outer;
  Polygons holes;
  HPolygon hpoly;
  bool got_first_hpoly(false);
  double number;
  // read the file line by line until EOF
  while (std::getline(is,line)) {
    std::istringstream iss(line);
    std::string cell;
    double number;
    // Get first cell on the current line
    std::getline(iss,cell,sep);
    if (cell.compare("NA")==0) { //line "NA,NA" -> new hole or new holed polygon
      if (poly.orientation()==CGAL::COUNTERCLOCKWISE) { // new holed polygon
	if (got_first_hpoly) {
	  /* This is not the first holed polygon we are reading, store the 
	     previously read one into ps */
	  hpoly = HPolygon(outer,holes.begin(),holes.end());
	  ps.push_back(hpoly);
	} else {
	  got_first_hpoly = true;
	}
	// push poly into outer
	outer = poly;
	// clear the series of holes
	holes.clear();
      } else { // feed holes with poly into 
	holes.push_back(poly);
      }
      poly.clear();
    } else { // pair of coordinates
      number = atof(cell.c_str()); // string -> double
      NT x(number);
      std::getline(iss,cell,sep); // read y-coordinate
      number = atof(cell.c_str());
      NT y(number);
      poly.push_back(Point2(x,y)); // push new point to poly
    }
  }
  // push last read data to the list of holed polygons
  hpoly = HPolygon(outer,holes.begin(),holes.end());
  ps.push_back(hpoly);
  return is;
}
TTessel::Halfedge_handle find_halfedge(TTessel &tes,Point2 p1,Point2 p2) {
  TTessel::Halfedge_handle e;
  for (TTessel::Halfedge_iterator hi=tes.halfedges_begin();
       hi!=tes.halfedges_end();hi++) {
    if (check_point_closeness(hi->source()->point(),p1) && 
	check_point_closeness(hi->target()->point(),p2)) {
      e = hi;
      return e;
    }
   }
  return NULL_HALFEDGE_HANDLE;
}
HPolygons holed_polygons() {
  Point2 outerVertices[6] = {
    Point2(0,0),
    Point2(6,3),
    Point2(14,0),
    Point2(9,10),
    Point2(12,15),
    Point2(0,15)};
  Polygon outer(outerVertices,outerVertices+6);
  Polygon hole;
  hole.push_back(Point2(3,4));
  hole.push_back(Point2(2,5));
  hole.push_back(Point2(3,7));
  hole.push_back(Point2(4,7));
  Polygons holes;
  holes.push_back(hole);
  hole.clear();
  hole.push_back(Point2(7,4));
  hole.push_back(Point2(7,5));
  hole.push_back(Point2(8,7));
  hole.push_back(Point2(8,5));
  hole.push_back(Point2(10,6));
  hole.push_back(Point2(10,4));
  holes.push_back(hole);
  hole.clear();
  hole.push_back(Point2(2,8 ));
  hole.push_back(Point2(2,13));
  hole.push_back(Point2(6,13));
  hole.push_back(Point2(4,11));
  hole.push_back(Point2(7,10));
  hole.push_back(Point2(7,9 ));
  hole.push_back(Point2(3,10));
  holes.push_back(hole);
  HPolygon dom_el(outer,holes.begin(),holes.end());
  HPolygons res;
  res.push_back(dom_el);
  outer.clear();
  outer.push_back(Point2(14,3 ));
  outer.push_back(Point2(21,3 ));
  outer.push_back(Point2(21,15));
  outer.push_back(Point2(14,15));
  outer.push_back(Point2(11,10));
  outer.push_back(Point2(14,5 ));
  holes.clear();
  hole.clear();
  hole.push_back(Point2(19,4));
  hole.push_back(Point2(19,5));
  hole.push_back(Point2(20,5));
  hole.push_back(Point2(20,4));
  holes.push_back(hole);
  hole.clear();
  hole.push_back(Point2(15,6));
  hole.push_back(Point2(17,8));
  hole.push_back(Point2(18,6));
  holes.push_back(hole);
  hole.clear();
  hole.push_back(Point2(13,9 ));
  hole.push_back(Point2(13,10));
  hole.push_back(Point2(14,10));
  hole.push_back(Point2(14,9 ));
  holes.push_back(hole);
  hole.clear();
  hole.push_back(Point2(17,10));
  hole.push_back(Point2(15,12));
  hole.push_back(Point2(17,14));
  hole.push_back(Point2(17,12));
  hole.push_back(Point2(19,12));
  holes.push_back(hole);
  dom_el = HPolygon(outer,holes.begin(),holes.end());
  res.push_back(dom_el);
  return res;
}
HPolygons holed_polygons_for_flips() {
  Point2 outerVertices[4] = {
    Point2(0,0),
    Point2(14,0),
    Point2(14,11),
    Point2(0,11)};
  Polygon outer(outerVertices,outerVertices+4);
  Polygon hole_C;
  hole_C.push_back(Point2(2,2));
  hole_C.push_back(Point2(2,9));
  hole_C.push_back(Point2(7,9));
  hole_C.push_back(Point2(7,8));
  hole_C.push_back(Point2(3,8));
  hole_C.push_back(Point2(3,3));
  hole_C.push_back(Point2(7,3));
  hole_C.push_back(Point2(7,2));
  Polygons holes;
  holes.push_back(hole_C);
  Polygon hole_O;
  hole_O.push_back(Point2(7,5));
  hole_O.push_back(Point2(7,7));
  hole_O.push_back(Point2(9,7));
  hole_O.push_back(Point2(9,5));
  holes.push_back(hole_O);
  HPolygon hpoly(outer,holes.begin(),holes.end());
  HPolygons res;
  res.push_back(hpoly);
  return res;
}
bool check_point_closeness(Point2 a, Point2 b, double tolerance) {
  boost::test_tools::percent_tolerance_t<double> tol(tolerance);
  boost::test_tools::close_at_tolerance<double> tester(tol);
  double ax = CGAL::to_double(a.x()), ay = CGAL::to_double(a.y());
  double bx = CGAL::to_double(b.x()), by = CGAL::to_double(b.y());
  if (tester(ax,bx) && tester(ay,by))
    return true;
  else
    return false;
}
/* Comparison of two polygons up to a cyclic permutation of their
   vertices. The polygons are not assumed to be simple. Still, for
   registration of both polygons, it is assumed that there exists in
   the first polygon a vertex that is not repeated. */
bool check_polygon_closeness(Polygon a,Polygon b,double tolerance) {
  if (a.size()!=b.size()) {
    return false;
  }
  Polygon::Vertex_circulator vc_b = b.vertices_circulator(),
    vc_a = a.vertices_circulator(), vc_a0, done, search_a0_done;;
  // Search for a vertex in polygon a that is not repeated
  vc_a0 = vc_a;
  search_a0_done = vc_a0;
  bool found_a0 = false;
  do { // loop on vc_a0 for searching a unique vertex on polygon a
    vc_a = vc_a0;
    done = vc_a;
    while (++vc_a!=done) { // search vc_a in polygon a using vc_a
      if (check_point_closeness(*vc_a,*vc_a0,tolerance)) break;
    }
    if (vc_a==done) { // vc_a0 is not repeated
      found_a0 = true;
      break;
    } else { // look for the next vertex
      vc_a0++;
    }
  } while (vc_a0!=search_a0_done);
  if (!found_a0) {
    throw std::invalid_argument("no unique vertex on first polygon");
  }
  vc_a = vc_a0;
  done = vc_b;
  found_a0 = false;
  do {
    if (check_point_closeness(*vc_b,*vc_a,tolerance)) {
      found_a0 = true;
      break;
    }
    vc_b++;
  } while(vc_b!=done);
  if (!found_a0) return false;
  done = vc_a;
  bool res = true;
  do {
    res = res & check_point_closeness(*vc_b,*vc_a,tolerance);
    if (!res) return false;
    vc_a++;
    vc_b++;
  } while (vc_a!=done);
  return true;
}    
bool check_holed_polygon_closeness(HPolygon a, HPolygon b) {
  if(!check_polygon_closeness(a.outer_boundary(),b.outer_boundary())) {
    return false;
  }
  Polygons holes_a(a.holes_begin(),a.holes_end()),
    holes_b(b.holes_begin(),b.holes_end());
  if (holes_a.size()!=holes_b.size()) return false;
  /* Hole sorting below is not very robust. Suppose there are in holed
     polygon a two holes, say h1 and h2, with equal upper left vertex,
     same pair of holes h1 and h2 in b. After sorting, h2 may follow
     h1 in a and vice-versa in b. Because of the inversion, a and b
     will be considered as different while they might not be.

     In practice, this may not happen since upper left vertices of
     different holes should not be close. */
  std::sort(holes_a.begin(),holes_a.end(),llt);
  std::sort(holes_b.begin(),holes_b.end(),llt);
  return std::mismatch(holes_a.begin(),holes_a.end(),holes_b.begin(),
		       check_polygon_closeness_default).first==holes_a.end();
}
bool match(HPolygon &x,HPolygons &Y) {
  for (Size i=0;i<Y.size();i++) {
    if (check_holed_polygon_closeness(x,Y[i])) {
      Y.erase(Y.begin()+i);
      return true;
    }
  }
  return false;
}
bool is_point_outside(Point2 p, HPolygon hpoly) {
  bool beyond_outer_boundary = 
    hpoly.outer_boundary().bounded_side(p)==CGAL::ON_UNBOUNDED_SIDE;
  bool inside_hole = false;
  for (HPolygon::Hole_const_iterator hi=hpoly.holes_begin();
       hi!=hpoly.holes_end();hi++) {
    inside_hole = inside_hole | hi->bounded_side(p)==CGAL::ON_BOUNDED_SIDE;
    if (inside_hole) break;
  }
  return beyond_outer_boundary | inside_hole;
}
bool is_point_outside(Point2 p, HPolygons hpolys) {
  bool res = true;
  for (HPolygons::iterator hi=hpolys.begin();hi!=hpolys.end();hi++) {
    res = res & is_point_outside(p,*hi);
    if (!res) break;
  }
  return res;
}
template <class ModifClass>
bool prediction_is_right(TTessel &tes,
			 ModifClass &modif,
			 HPolygons &deleted_faces_unpredicted, 
			 HPolygons &deleted_faces_unrealized, 
			 HPolygons &added_faces_unpredicted, 
			 HPolygons &added_faces_unrealized,
			 std::ostringstream &info) {
  bool is_right = true;
  deleted_faces_unrealized.clear();
  deleted_faces_unpredicted.clear();
  added_faces_unrealized.clear();
  added_faces_unpredicted.clear();
  ModList predict;
  HPolygons past_polygons = get_faces(tes);
  try {
    predict = modif.modified_elements();
  } catch(std::logic_error const& exception) {
    is_right =  false;
    info << "Error occured when getting predictions. "
	 << exception.what();
    return false;
  }
  deleted_faces_unrealized = predict.del_faces;
  added_faces_unrealized = predict.add_faces;
  tes.update(modif);
  HPolygons new_polygons = get_faces(tes);
  HPolygons removed_polygons = set_diff(past_polygons,new_polygons);
  HPolygons added_polygons = set_diff(new_polygons,past_polygons);
  // Compare removed faces with prediction
  for (HPolygons::const_iterator p=removed_polygons.begin();
       p!=removed_polygons.end();p++) {
    HPolygon curr_poly(*p);
    bool found = match(curr_poly,deleted_faces_unrealized);
    if (!found) {
      is_right = false;
      HPolygon unpredicted(*p);
      deleted_faces_unpredicted.push_back(unpredicted);
    }
  }
  if (deleted_faces_unrealized.size()>0) {
    is_right = false;
  }
  // Compare new faces with prediction
  for (HPolygons::const_iterator p=added_polygons.begin();
       p!=added_polygons.end();p++) {
    HPolygon curr_poly(*p);
    bool found = match(curr_poly,added_faces_unrealized);
    if (!found) {
      is_right = false;
      HPolygon unpredicted(*p);
      added_faces_unrealized.push_back(unpredicted);
    }
  }
  if (added_faces_unpredicted.size()>0) {
    is_right = false;
  }
  return is_right;
}
void process_wrong_predictions(HPolygons deleted_faces_unpredicted,
			       HPolygons deleted_faces_unrealized,
			       HPolygons added_faces_unpredicted,
			       HPolygons added_faces_unrealized,
			       const char old_face_unpredicted_path[],
			       const char old_face_unrealized_path[],
			       const char new_face_unpredicted_path[],
			       const char new_face_unrealized_path[],
			       std::ostringstream &info) {
      if (deleted_faces_unpredicted.size()>0) {
	info << deleted_faces_unpredicted.size() 
		<< " unpredicted removed face(s)";
	std::ofstream output_file;
	output_file.open(old_face_unpredicted_path,std::ios::out);
    	output_file << "x,y" << std::endl;
    	output_file << deleted_faces_unpredicted;
    	output_file.close();

      }
      if (deleted_faces_unrealized.size()>0) {
	info << deleted_faces_unrealized.size() 
		<< " face(s) predicted as being removed was/were not";
	std::ofstream output_file;
	output_file.open(old_face_unrealized_path,std::ios::out);
    	output_file << "x,y" << std::endl;
    	output_file << deleted_faces_unrealized;
    	output_file.close();
      }
      if (added_faces_unpredicted.size()>0) {
	info << added_faces_unpredicted.size() 
		<< " unpredicted added face(s)";
	std::ofstream output_file;
	output_file.open(new_face_unpredicted_path,std::ios::out);
    	output_file << "x,y" << std::endl;
    	output_file << added_faces_unpredicted;
    	output_file.close();
      }
      if (added_faces_unrealized.size()>0) {
	info << added_faces_unrealized.size() 
		<< " face(s) predicted as being added was/were not";
	std::ofstream output_file;
	output_file.open(new_face_unrealized_path,std::ios::out);
    	output_file << "x,y" << std::endl;
    	output_file << added_faces_unrealized;
    	output_file.close();
      }
}

void check_predictions_4_random_smf(unsigned seed,unsigned n,unsigned m,
		       HPolygons &dom) {
  rnd = new CGAL::Random(seed);
  TTessel tesl;
  tesl.insert_window(dom);
  SMFChain::ModType modif_type;
  for (int i=0;i<n+m;i++) {
    //std::clog << "Iteration " << i << std::endl;//debug
    std::ofstream tessel_file;
    tessel_file.open("ttessellation.txt",std::ios::out);
    tesl.write(tessel_file);
    tessel_file.close();
    TTessel::Split s;
    TTessel::Merge m;
    TTessel::Flip f;
    std::ostringstream details(std::ostringstream::out);
    details << "iteration " << i << ", ";
    if (i<n) { // only splits
      modif_type = SMFChain::SPLIT;
      s = tesl.propose_split();
      details <<  "split. Splitting segment: "
      	      << "(" << s.get_p1() << ")-(" << s.get_p2() << ")";
    } else { // splits or merges
      switch(i%3) {
      case 0:
	modif_type = SMFChain::SPLIT;
	s = tesl.propose_split();
	details <<  "split. Splitting segment: "
		<< "(" << s.get_p1() << ")-(" << s.get_p2() << "). ";
	break;
      case 1:
	if (tesl.number_of_non_blocking_segments()==0)
	  continue;
	modif_type = SMFChain::MERGE;
	m = tesl.propose_merge();
	details << "merge. Removed segment: "
		<< "(" << m.get_e()->source()->point() << ")-(" 
		<< m.get_e()->target()->point() << "). ";
	break;
      default:
	if (tesl.number_of_blocking_segments()==0)
	  continue;
	modif_type = SMFChain::FLIP;
	f = tesl.propose_flip();
	details << "flip. Removed edge: "
		<< "(" << f.get_e1()->source()->point() << ")-("
		<< f.get_e1()->target()->point() << "). Inserted edge: ("
		<< f.get_e1()->source()->point() << ")-("
		<< f.get_p2() << "). ";
      }
    }
    HPolygons deleted_faces_unpredicted, deleted_faces_unrealized, 
      added_faces_unpredicted, added_faces_unrealized;
    bool ok;
    switch(modif_type) {
    case SMFChain::SPLIT : 
      ok = prediction_is_right<TTessel::Split>(tesl, s, 
					       deleted_faces_unpredicted, 
					       deleted_faces_unrealized, 
					       added_faces_unpredicted, 
					       added_faces_unrealized,
					       details);
      break;
    case SMFChain::MERGE : 
      ok = prediction_is_right<TTessel::Merge>(tesl, m, 
					       deleted_faces_unpredicted, 
					       deleted_faces_unrealized, 
					       added_faces_unpredicted, 
					       added_faces_unrealized,
					       details);
      break;
    case SMFChain::FLIP :
      ok = prediction_is_right<TTessel::Flip>(tesl,f,
					       deleted_faces_unpredicted, 
					       deleted_faces_unrealized, 
					       added_faces_unpredicted, 
					       added_faces_unrealized,
					       details);
    }
    if (!ok) {
      process_wrong_predictions(deleted_faces_unpredicted, 
				deleted_faces_unrealized, 
				added_faces_unpredicted, 
				added_faces_unrealized,
				"smf_old_face_unpredicted.csv",
				"smf_old_face_unrealized.csv",
				"smf_new_face_unpredicted.csv",
				"smf_new_face_unrealized.csv",details);
      BOOST_FAIL(details.str());
    }
    if (!tesl.is_a_T_tessellation(true)) {
      details << "Test halted because the tessellation is not a T-tessellation";
      BOOST_FAIL(details.str());
    }
  }
  delete rnd;
}

BOOST_AUTO_TEST_CASE(clip_segment_by_holed_polygon) {
  HPolygon dom = holed_polygons()[0];
  Segment seg(Point2(-1,8),Point2(3,12));
  Segment cseg = clip_segment_by_polygon(seg,dom);
  Segment expected_cseg(Point2(0,9),Point2(2,11));
  BOOST_CHECK_MESSAGE(cseg==expected_cseg || cseg==expected_cseg.opposite(),
		      "The clipped segment is not correct. It is ("
		      << cseg.source() << ")-( " << cseg.target() << ")"
		      << " instead of ("
		      << expected_cseg.source() << ")-(" 
		      << expected_cseg.target() << ")");
  // test case where the segment is outside the holed polygon
  seg = Segment(Point2(-1,-1),Point2(-2,-2));
  BOOST_CHECK_THROW(cseg = clip_segment_by_polygon(seg,dom),
		    std::domain_error);
}

BOOST_AUTO_TEST_CASE(clip_segment_by_holed_polygons) {
  HPolygons dom = holed_polygons();
  Segment seg(Point2(19/2,9/2),Point2(39/2,9/2));
  Segment cseg = clip_segment_by_polygon(seg,dom);
  Segment expected_cseg(Point2(14,9/2),Point2(19,9/2));
  BOOST_CHECK_MESSAGE(cseg==expected_cseg || cseg==expected_cseg.opposite(),
		      "The clipped segment is not correct. It is ("
		      << cseg.source() << ")-( " << cseg.target() << ")"
		      << " instead of ("
		      << expected_cseg.source() << ")-(" 
		      << expected_cseg.target() << ")");
  // test case where the segment is outside both holed polygons
  seg = Segment(Point2(-1,-1),Point2(-2,-2));
  BOOST_CHECK_THROW(cseg = clip_segment_by_polygon(seg,dom),
		    std::domain_error);
  // test case where the segment hits only one holed polygon
  seg = Segment(Point2(19/2,9/2),Point2(13,9/2));
  BOOST_CHECK_NO_THROW(cseg = clip_segment_by_polygon(seg,dom));
}

BOOST_AUTO_TEST_CASE(insert_segment) {
  HPolygons dom = holed_polygons();
  rnd = new CGAL::Random(0);
  TTessel tesl;
  tesl.insert_window(dom);
  /* The line segment to be inserted extends beyond the domain outer
     boundary at one end and inside a hole at the other end. Should be
     clipped at both ends. */
  Point2 seg_p1(3,12), seg_p2(3,16);
  Segment seg(seg_p1,seg_p2);
  Segment seg_test(seg);
  seg_test = seg;
  // Segment seg(Point2(3,12),Point2(3,16));
  tesl.insert_segment(seg);
  TTessel::Halfedge_handle e = find_halfedge(tesl,Point2(3,13),Point2(3,15));
  BOOST_CHECK_MESSAGE(e!=NULL_HALFEDGE_HANDLE,
		      "Segment should have been clipped to edge "
		      << "(" << Point2(3,13) << ")-(" << Point2(3,15)
		      << "). But that edge was not found");
  TTessel::Halfedge_handle e_next = e->get_next_hf();
  BOOST_CHECK_MESSAGE(e_next==NULL_HALFEDGE_HANDLE,
		      "the upper part of the inserted segment extends to the "
		      << "unbounded side of the domain outer boundary");
  TTessel::Halfedge_handle e_prev = e->get_prev_hf();
  BOOST_CHECK_MESSAGE(e_prev==NULL_HALFEDGE_HANDLE,
		      "the lower part of the insert segment extends into "
		      << "a domain hole");
  /* The segment to be inserted has one end inside the domain and the
     other end beyond the domain outer boundary. Should be clipped on
     one side. */
  seg = Segment(Point2(1,2),Point2(5,2));
  tesl.insert_segment(seg);
  e = find_halfedge(tesl,Point2(1,2),Point2(4,2));
  BOOST_CHECK_MESSAGE(e!=NULL_HALFEDGE_HANDLE,
		      "Segment should have been clipped to edge "
		      << "(" << Point2(1,2) << ")-(" << Point2(4,2)
		      << "). But that edge was not found");
  e_next = e->get_next_hf();
  BOOST_CHECK_MESSAGE(e_next==NULL_HALFEDGE_HANDLE,
		      "the inserted segment extends too much towards right");
  e_prev = e->get_prev_hf();
  BOOST_CHECK_MESSAGE(e_prev==NULL_HALFEDGE_HANDLE,
		      "the inserted segment extends too much towards left");
}
  
BOOST_AUTO_TEST_CASE(remove_ivertices) {
  HPolygons dom = holed_polygons();
  rnd = new CGAL::Random(0);
  TTessel tesl;
  tesl.insert_window(dom);
  /* Insert a segment with both free ends, not too far from the domain
     boundary. Both ends are I-vertices. They should be removed by
     extending the segment to the domain boundary on both sides. */
  Segment seg(Point2(5,3),Point2(5,9));
  tesl.insert_segment(seg);
  tesl.remove_ivertices();
  TTessel::Halfedge_handle e = find_halfedge(tesl,Point2(5,2.5),Point2(5,9.5));
  BOOST_CHECK_MESSAGE(e!=NULL_HALFEDGE_HANDLE,
		      "Segment should have been extended to edge "
		      << "(" << Point2(2,2.5) << ")-(" << Point2(5,9.5)
		      << "). But that edge was not found");
  TTessel::Halfedge_handle e_next = e->get_next_hf();
  BOOST_CHECK_MESSAGE(e_next==NULL_HALFEDGE_HANDLE,
		      "the segment extends too much upside");
  TTessel::Halfedge_handle e_prev = e->get_prev_hf();
  BOOST_CHECK_MESSAGE(e_prev==NULL_HALFEDGE_HANDLE,
		      "the segments extends too much downside");
  /* Case where a segment must be shortened at one end */
  tesl.clear();
  tesl.insert_window(dom);
  e = find_halfedge(tesl,Point2(0,15),Point2(0,0));
  TTessel::Split split(e,0.5,CGAL_PI/2);
  tesl.update(split);
  seg = Segment(Point2(8.5,15),Point2(8.5,7));
  tesl.insert_segment(seg);
  e = find_halfedge(tesl,Point2(8.5,7),Point2(8.5,7.5));
  if (e==NULL_HALFEDGE_HANDLE) {
    std::cerr << "Insertion of segment (" << seg.source() << ")-("
	      << seg.target() << ") should have generated ";
    std::cerr << "halfedge ("
	      << Point2(8.5,7) << ")-(" << Point2(8.5,7.5) << "). "
	      << "But that halfege was not found." << std::endl;
  }
  tesl.remove_ivertices();
  e = find_halfedge(tesl,Point2(8.5,7.5),Point2(8.5,7));
  BOOST_CHECK_MESSAGE(e==NULL_HALFEDGE_HANDLE,
		      "Halfedge (" <<  Point2(8.5,7.5) << ")-(" 
		      << Point2(8.5,7) << ") has not been removed");
}

BOOST_AUTO_TEST_CASE(split_constructor) {
  HPolygons dom = holed_polygons();
  rnd = new CGAL::Random(0);
  TTessel tesl;
  tesl.insert_window(dom);
  // First test: the split segment joins two points on the outer boundary
  // find halfedge (0,15)-(0,0)
  TTessel::Halfedge_handle e1 = find_halfedge(tesl,Point2(0,15),Point2(0,0));
  // split starting from (0,1) horizontally
  TTessel::Split s1(e1,14.0/15,CGAL_PI/2);
  BOOST_CHECK_MESSAGE(s1.get_p1()==Point2(0,1),
		      "split segment joining two points on the outer boundary:"
		      << "it starts from (" << s1.get_p1() << ")  instead of " 
		      << "(" << Point2(0,1) << ")");
  // find halfedge (0,0)-(6,3)
  TTessel::Halfedge_handle e2 = find_halfedge(tesl,Point2(0,0),Point2(6,3));
  BOOST_CHECK_MESSAGE(s1.get_e2()==e2,
		      "split segment joining two points on the outer boundary:"
		      << " it ends on edge (" << s1.get_e2()->source()->point()
		      << ")-(" << s1.get_e2()->target()->point() << ") "
		      << "instead of edge (" << e2->source()->point()
		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(s1.get_p2()==Point2(2,1),
		      "split segment joining two points on the outer boundary:"
		      << "it ends at (" << s1.get_p2() << ")  instead of " 
		      << "(" << Point2(2,1) << ")");
  // Second test: the split segment joins two points on the same inner boundary
  e1 = find_halfedge(tesl,Point2(6,13),Point2(4,11));
  TTessel::Split s2(e1,1.0/2,CGAL_PI/4);
  BOOST_CHECK_MESSAGE(check_point_closeness(s2.get_p1(),Point2(5,12)),
		      "split segment joining points on the same inner boundary:"
		      << "it starts from (" << s2.get_p1() << ")  instead of " 
		      << "(" << Point2(5,12) << ")");
  e2 = find_halfedge(tesl,Point2(4,11),Point2(7,10));
  BOOST_CHECK_MESSAGE(s2.get_e2()==e2,
		      "split segment joining two points on the same inner"
		      << " boundary: it ends on edge ("
		      << s2.get_e2()->source()->point()
		      << ")-(" << s2.get_e2()->target()->point() << ") "
		      << "instead of edge (" << e2->source()->point()
		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(check_point_closeness(s2.get_p2(),Point2(5,11-1.0/3)),
		      "split segment joining two points on the same inner "
		      << "boundary: it ends at (" << s2.get_p2() 
		      << ")  instead of " << "(" << Point2(5,11-1.0/3) << ")");
  //Third test: the split segment joins the outer boundary to an inner one
  e1 = find_halfedge(tesl,Point2(12,15),Point2(0,15));
  TTessel::Split s3(e1,8.0/12,CGAL_PI/2);
  BOOST_CHECK_MESSAGE(check_point_closeness(s3.get_p1(),Point2(4,15)),
		      "split segment joining the outer boundary to an inner "
		      << "boundary:"
		      << "it starts from (" << s3.get_p1() << ")  instead of " 
		      << "(" << Point2(4,15) << ")");
  e2 = find_halfedge(tesl,Point2(2,13),Point2(6,13));
  BOOST_CHECK_MESSAGE(s3.get_e2()==e2,
		      "split segment joining the outer boundary to an inner "
		      << " boundary: it ends on edge ("
		      << s3.get_e2()->source()->point()
		      << ")-(" << s3.get_e2()->target()->point() << ") "
		      << "instead of edge (" << e2->source()->point()
		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(check_point_closeness(s3.get_p2(),Point2(4,13)),
		      "split segment joining the outer boundary to an inner "
		      << "boundary: it ends at (" << s3.get_p2() 
		      << ")  instead of " << "(" << Point2(4,13) << ")");
  // Fourth and last test: the split segment joins two different inner 
  // boundaries
  e1 = find_halfedge(tesl,Point2(4,7),Point2(3,4));
  TTessel::Split s4(e1,1.0/3,CGAL_PI/2+atan(1.0/3));
  BOOST_CHECK_MESSAGE(check_point_closeness(s4.get_p1(),Point2(3+2.0/3,6)),
	      "split segment joining two different inner boundaries:"
		      << " it starts from (" << s4.get_p1() << ")  instead of " 
		      << "(" << Point2(3+2.0/3,6) << ")");
  e2 = find_halfedge(tesl,Point2(7,5),Point2(8,7));
  BOOST_CHECK_MESSAGE(s4.get_e2()==e2,
		      "split segment joining two different inner boundaries: "
		      << "it ends on edge ("
		      << s4.get_e2()->source()->point()
		      << ")-(" << s4.get_e2()->target()->point() << ") "
		      << "instead of edge (" << e2->source()->point()
		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(check_point_closeness(s4.get_p2(),Point2(7.5,6)),
		      "split segment joining two different inner boundaries: "
		      << "it ends at (" << s4.get_p2() 
		      << ")  instead of " << "(" << Point2(7.5,6) << ")");
}

BOOST_AUTO_TEST_CASE(flip_constructor) {
  HPolygons dom = holed_polygons_for_flips();
  rnd = new CGAL::Random(0);
  TTessel tesl;
  tesl.insert_window(dom);
  // FIRST TEST
  TTessel::Halfedge_handle s1_e1 = find_halfedge(tesl,Point2(0,11),
						 Point2(0,0));
  TTessel::Split s1(s1_e1,10./11,CGAL_PI/2);
  tesl.update(s1);
  TTessel::Halfedge_handle s2_e1 = find_halfedge(tesl,Point2(0,0),
						 Point2(14,0));
  TTessel::Split s2(s2_e1,11./14,CGAL_PI/2);
  tesl.update(s2);
  TTessel::Halfedge_handle e = find_halfedge(tesl,Point2(11,1),Point2(14,1));
  TTessel::Flip f(e);
  TTessel::Halfedge_handle e2 = find_halfedge(tesl,Point2(14,11),Point2(0,11));
  BOOST_CHECK_MESSAGE(f.get_e2()==e2,
		      "flip test 1,"
		      << " the split edge is (" << f.get_e2()->source()->point()
		      << ")-(" << f.get_e2()->target()->point() << ") "
		      << "instead of edge (" << e2->source()->point()
		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(check_point_closeness(f.get_p2(),Point2(11,11)),
		      "flip test 1, new vertex is "
		      << " (" << f.get_p2() 
		      << ")  instead of " << "(" << Point2(11,11) << ")");
  // TEST 2
  tesl.clear(false);
  s1_e1 = find_halfedge(tesl,Point2(14,0), Point2(14,11));
  TTessel::Split s1a(s1_e1,6./11,CGAL_PI/2);
  tesl.update(s1a);
  s2_e1 = find_halfedge(tesl,Point2(0,0), Point2(14,0));
  TTessel::Split s2a(s2_e1,11./14,CGAL_PI/2);
  tesl.update(s2a);
  e = find_halfedge(tesl,Point2(11,6),Point2(9,6));
  TTessel::Flip fa(e);
  e2 = find_halfedge(tesl,Point2(14,11),Point2(0,11));
  BOOST_CHECK_MESSAGE(fa.get_e2()==e2,
		      "flip test 2,"
		      << " the split edge is (" 
		      << fa.get_e2()->source()->point()
		      << ")-(" << fa.get_e2()->target()->point() << ") "
		      << "instead of edge (" << e2->source()->point()
		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(check_point_closeness(fa.get_p2(),Point2(11,11)),
		      "flip test 2, new vertex is "
		      << " (" << fa.get_p2() 
		      << ")  instead of " << "(" << Point2(11,11) << ")");
  // TEST 3
  tesl.clear(false);
  s1_e1 = find_halfedge(tesl,Point2(3,3), Point2(7,3));
  TTessel::Split s1b(s1_e1,1./2,CGAL_PI/2);
  tesl.update(s1b);
  s2_e1 = find_halfedge(tesl,Point2(14,0), Point2(14,11));
  TTessel::Split s2b(s2_e1,4./11,CGAL_PI/2);
  tesl.update(s2b);
  e = find_halfedge(tesl,Point2(5,4),Point2(5,8));
  TTessel::Flip fb(e);
  e2 = find_halfedge(tesl,Point2(3,8),Point2(3,3));
  BOOST_CHECK_MESSAGE(fb.get_e2()==e2,
		      "flip test 3,"
		      << " the split edge is (" 
		      << fb.get_e2()->source()->point()
		      << ")-(" << fb.get_e2()->target()->point() << ") "
		      << "instead of edge (" << e2->source()->point()
		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(check_point_closeness(fb.get_p2(),Point2(3,4)),
		      "flip test 3, new vertex is "
		      << " (" << fb.get_p2() 
		      << ")  instead of " << "(" << Point2(3,4) << ")");
  // TEST 4
  tesl.clear(false);
  s1_e1 = find_halfedge(tesl,Point2(3,3), Point2(7,3));
  TTessel::Split s1c(s1_e1,1./2,CGAL_PI/2);
  tesl.update(s1c);
  s2_e1 = find_halfedge(tesl,Point2(7,5), Point2(7,7));
  TTessel::Split s2c(s2_e1,1.0/2,CGAL_PI/2);
  tesl.update(s2c);
  e = find_halfedge(tesl,Point2(5,6),Point2(5,8));
  TTessel::Flip fc(e);
  e2 = find_halfedge(tesl,Point2(3,8),Point2(3,3));
  BOOST_CHECK_MESSAGE(fc.get_e2()==e2,
		      "flip test 4,"
		      << " the split edge is (" 
		      << fc.get_e2()->source()->point()
		      << ")-(" << fc.get_e2()->target()->point() << ") "
		      << "instead of edge (" << e2->source()->point()
		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(check_point_closeness(fc.get_p2(),Point2(3,6)),
		      "flip test 4, new vertex is "
		      << " (" << fc.get_p2() 
		      << ")  instead of " << "(" << Point2(3,6) << ")");
  // TEST 5
  tesl.clear(false);
  s1_e1 = find_halfedge(tesl,Point2(3,3), Point2(7,3));
  TTessel::Split s1d(s1_e1,1./4,CGAL_PI/2);
  tesl.update(s1d);
  s2_e1 = find_halfedge(tesl,Point2(3,8), Point2(3,3));
  TTessel::Split s2d(s2_e1,2.0/5,3*CGAL_PI/4);
  tesl.update(s2d);
  e = find_halfedge(tesl,Point2(4,7),Point2(4,8));
  TTessel::Flip fd(e);
  e2 = find_halfedge(tesl,Point2(7,8),Point2(4,8));
  BOOST_CHECK_MESSAGE(fd.get_e2()==e2,
  		      "flip test 5,"
  		      << " the split edge is (" 
  		      << fd.get_e2()->source()->point()
  		      << ")-(" << fd.get_e2()->target()->point() << ") "
  		      << "instead of edge (" << e2->source()->point()
  		      << ")-( " << e2->target()->point() << ") ");
  BOOST_CHECK_MESSAGE(check_point_closeness(fd.get_p2(),Point2(5,8)),
  		      "flip test 5, new vertex is "
  		      << " (" << fd.get_p2() 
  		      << ")  instead of " << "(" << Point2(5,8) << ")");
}

BOOST_AUTO_TEST_CASE(flip_modified_elements) {
  HPolygons dom = holed_polygons_for_flips();
  rnd = new CGAL::Random(0);
  TTessel tesl;
  tesl.insert_window(dom);
  // FIRST TEST
  TTessel::Halfedge_handle s1_e1 = find_halfedge(tesl,Point2(0,11),
						 Point2(0,0));
  TTessel::Split s1(s1_e1,10./11,CGAL_PI/2);
  tesl.update(s1);
  TTessel::Halfedge_handle s2_e1 = find_halfedge(tesl,Point2(0,0),
						 Point2(14,0));
  TTessel::Split s2(s2_e1,11./14,CGAL_PI/2);
  tesl.update(s2);
  TTessel::Halfedge_handle e = find_halfedge(tesl,Point2(11,1),Point2(14,1));
  TTessel::Flip f(e);
  HPolygons deleted_faces_unpredicted, deleted_faces_unrealized, 
    added_faces_unpredicted, added_faces_unrealized;
  std::ostringstream details(std::ostringstream::out);
  bool ok = prediction_is_right<TTessel::Flip>(tesl,f,			
					       deleted_faces_unpredicted, 
					       deleted_faces_unrealized, 
					       added_faces_unpredicted, 
					       added_faces_unrealized,
					       details);
  if (!ok) {
    process_wrong_predictions(deleted_faces_unpredicted, 
			      deleted_faces_unrealized, 
			      added_faces_unpredicted, 
			      added_faces_unrealized,
			      "flip_old_face_unpredicted.csv",
			      "flip_old_face_unrealized.csv",
			      "flip_new_face_unpredicted.csv",
			      "flip_new_face_unrealized.csv",details);
    BOOST_FAIL(details);
    
  }
  // TEST 2
  tesl.clear(false);
  s1_e1 = find_halfedge(tesl,Point2(14,0), Point2(14,11));
  TTessel::Split s1a(s1_e1,6./11,CGAL_PI/2);
  tesl.update(s1a);
  s2_e1 = find_halfedge(tesl,Point2(0,0), Point2(14,0));
  TTessel::Split s2a(s2_e1,11./14,CGAL_PI/2);
  tesl.update(s2a);
  e = find_halfedge(tesl,Point2(11,6),Point2(9,6));
  TTessel::Flip fa(e);
  details.str("");
  details.clear();
  ok = prediction_is_right<TTessel::Flip>(tesl,fa,			
					  deleted_faces_unpredicted, 
					  deleted_faces_unrealized, 
					  added_faces_unpredicted, 
					  added_faces_unrealized,
					  details);
  if (!ok) {
    process_wrong_predictions(deleted_faces_unpredicted, 
			      deleted_faces_unrealized, 
			      added_faces_unpredicted, 
			      added_faces_unrealized,
			      "flip_old_face_unpredicted.csv",
			      "flip_old_face_unrealized.csv",
			      "flip_new_face_unpredicted.csv",
			      "flip_new_face_unrealized.csv",details);
    BOOST_FAIL(details.str());
  }
  // TEST 3
  tesl.clear(false);
  s1_e1 = find_halfedge(tesl,Point2(3,3), Point2(7,3));
  TTessel::Split s1b(s1_e1,1./2,CGAL_PI/2);
  tesl.update(s1b);
  s2_e1 = find_halfedge(tesl,Point2(14,0), Point2(14,11));
  TTessel::Split s2b(s2_e1,4./11,CGAL_PI/2);
  tesl.update(s2b);
  e = find_halfedge(tesl,Point2(5,4),Point2(5,8));
  TTessel::Flip fb(e);
  details.str("");
  details.clear();
  ok = prediction_is_right<TTessel::Flip>(tesl,fb,			
					  deleted_faces_unpredicted, 
					  deleted_faces_unrealized, 
					  added_faces_unpredicted, 
					  added_faces_unrealized,
					  details);
  if (!ok) {
    process_wrong_predictions(deleted_faces_unpredicted, 
			      deleted_faces_unrealized, 
			      added_faces_unpredicted, 
			      added_faces_unrealized,
			      "flip_old_face_unpredicted.csv",
			      "flip_old_face_unrealized.csv",
			      "flip_new_face_unpredicted.csv",
			      "flip_new_face_unrealized.csv",details);
    BOOST_FAIL(details.str());
  }
  // TEST 4
  tesl.clear(false);
  s1_e1 = find_halfedge(tesl,Point2(3,3), Point2(7,3));
  TTessel::Split s1c(s1_e1,1./2,CGAL_PI/2);
  tesl.update(s1c);
  s2_e1 = find_halfedge(tesl,Point2(7,5), Point2(7,7));
  TTessel::Split s2c(s2_e1,1.0/2,CGAL_PI/2);
  tesl.update(s2c);
  e = find_halfedge(tesl,Point2(5,6),Point2(5,8));
  TTessel::Flip fc(e);
  details.str("");
  details.clear();
  ok = prediction_is_right<TTessel::Flip>(tesl,fc,			
					  deleted_faces_unpredicted, 
					  deleted_faces_unrealized, 
					  added_faces_unpredicted, 
					  added_faces_unrealized,
					  details);
  if (!ok) {
    process_wrong_predictions(deleted_faces_unpredicted, 
			      deleted_faces_unrealized, 
			      added_faces_unpredicted, 
			      added_faces_unrealized,
			      "flip_old_face_unpredicted.csv",
			      "flip_old_face_unrealized.csv",
			      "flip_new_face_unpredicted.csv",
			      "flip_new_face_unrealized.csv",details);
    BOOST_FAIL(details.str());
  }
  // TEST 5
  tesl.clear(false);
  s1_e1 = find_halfedge(tesl,Point2(3,3), Point2(7,3));
  TTessel::Split s1d(s1_e1,1./4,CGAL_PI/2);
  tesl.update(s1d);
  s2_e1 = find_halfedge(tesl,Point2(3,8), Point2(3,3));
  TTessel::Split s2d(s2_e1,2.0/5,3*CGAL_PI/4);
  tesl.update(s2d);
  e = find_halfedge(tesl,Point2(4,7),Point2(4,8));
  TTessel::Flip fd(e);
  details.str("");
  details.clear();
  ok = prediction_is_right<TTessel::Flip>(tesl,fd,			
					  deleted_faces_unpredicted, 
					  deleted_faces_unrealized, 
					  added_faces_unpredicted, 
					  added_faces_unrealized,
					  details);
  if (!ok) {
    process_wrong_predictions(deleted_faces_unpredicted, 
			      deleted_faces_unrealized, 
			      added_faces_unpredicted, 
			      added_faces_unrealized,
			      "flip_old_face_unpredicted.csv",
			      "flip_old_face_unrealized.csv",
			      "flip_new_face_unpredicted.csv",
			      "flip_new_face_unrealized.csv",details);
    BOOST_FAIL(details.str());
  }
  delete rnd;
}

BOOST_AUTO_TEST_CASE(squared_domain) {
  Polygon square;
  square.push_back(Point2(0,0));
  square.push_back(Point2(1,0));
  square.push_back(Point2(1,1));
  square.push_back(Point2(0,1));
  HPolygons dom;
  dom.push_back(HPolygon(square));
  check_predictions_4_random_smf(5,100,50,dom);
}

BOOST_AUTO_TEST_CASE(non_convex_holed_domain) {
  HPolygons dom = holed_polygons();
  /* It is important to keep the number of internal segments
     low. Otherwise there are too few faces with holes. With the
     parameters below, each case is observed at least about 40
     times. */
  check_predictions_4_random_smf(5,5,5000,dom);
}

BOOST_AUTO_TEST_CASE(squared_diameter_points) {
  typedef CGAL::Cartesian_d<double> Kd;
  typedef Kd::Point_d Point;
  int n = 100;
  CGAL::Random_points_in_ball_d<Point> gen (2, 5.0);
  Points sample;
  for (int k=0;k<n;k++) {
    Point p = *gen++;
    Point2 p2(p.x(),p.y());
    sample.push_back(p2);
  }
  NT diam2_rc = squared_diameter(sample.begin(),sample.end());
  NT diam2_bf = 0;
  for (int k=0;k<sample.size();k++) {
    for (int l=0;l<k;l++) {
      NT d = CGAL::squared_distance(sample[k],sample[l]);
      if (d>diam2_bf) diam2_bf = d;
    }
  }
  BOOST_CHECK_MESSAGE(diam2_rc==diam2_bf,"computed squared diameter"
		      << " is " << diam2_rc << " instead of " << diam2_bf);
}
		      
BOOST_AUTO_TEST_CASE(hausdorff_segments) {
  std::vector<Segment> P, Q;
  P.push_back(Segment(Point2(0,0),Point2(3,0)));
  P.push_back(Segment(Point2(3,0),Point2(0,3)));
  P.push_back(Segment(Point2(0,3),Point2(2,3)));
  P.push_back(Segment(Point2(2,3),Point2(0,1)));
  P.push_back(Segment(Point2(0,1),Point2(0,0)));
  Q.push_back(Segment(Point2(0,0),Point2(3,0)));
  Q.push_back(Segment(Point2(0,0),Point2(3,0)));
  Q.push_back(Segment(Point2(3,0),Point2(3,3)));
  Q.push_back(Segment(Point2(3,3),Point2(0,0)));
  Q.push_back(Segment(Point2(1,2),Point2(3,2)));
  NT shd = squared_Hausdorff_distance(P,Q);
  BOOST_CHECK_MESSAGE(shd==0,"temporary dummy test for debugging");
}
