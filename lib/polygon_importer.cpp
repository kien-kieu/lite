
// Line Tessellation (LiTe) library
// |||Development version
// Authors: Katarzyna Adamczyk and Kiên Kiêu.
// |||Copyright INRA 2006-yyyy.
// Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
// License: GPL v3.

#include "ttessel.h"

double infinity_double = std::numeric_limits<double>::infinity();

/** \brief Empty constructor */
PolygonImporter::PolygonImporter() {}
/** \brief Read input polygons
 *
 * Feed data member input_polygons from input stream. Expected format:
 * for each polygon, a line with numbers of vertices (separated by a white
 * space) followed by a series of xy-coordinates (one by line, separated by
 * a white space). The first line provides the number of vertices of the
 * polygon borders. The first border must be the outer boundaries. The
 * subsequent borders are inner boundaries. Vertices along the outer boundary
 * must be orderd anti-clockwise. Vertices along the inner boundaries must be
 * ordered clockwise.
 */
void PolygonImporter::read_polygons(std::istream& input) {
  std::string line;
  while (std::getline(input,line)) {
    std::istringstream size_iss(line);
    std::vector<Size> sizes;
    Size buf;
    while (size_iss>>buf) sizes.push_back(buf);
    Polygons borders;
    for (std::vector<Size>::iterator the_size=sizes.begin();
         the_size!=sizes.end();the_size++) {
      Points vertices;
      for (Size i=0;i<*the_size;i++) {
        std::getline(input,line);
        std::istringstream coord_iss(line);
        NT x, y;
        coord_iss >> x;
        coord_iss >> y;
        Point2 vertex(x,y);
        vertices.push_back(vertex);
      }
      Polygon border(vertices.begin(),vertices.end());
      borders.push_back(border);
    }
    HPolygon hpoly(borders[0],borders.begin()+1,borders.end());
    input_polygons.push_back(hpoly);
  }
}

/** \brief Return the polygon sides */ 
PolygonImporter::Sides PolygonImporter::get_polygon_sides() {
  Sides res;
  Size poly_index = 1;
  for (HPolygons::iterator hi=input_polygons.begin();hi!=input_polygons.end();
       hi++) {
    Polygon outer, inner;
    outer = hi->outer_boundary();
    for (Polygon::Edge_const_iterator ei=outer.edges_begin();
         ei!=outer.edges_end();ei++) {
      Side side;
      side.first = Segment(ei->source(),ei->target());
      side.second = poly_index;
      res.push_back(side);
    }
    for (HPolygon::Hole_const_iterator hoi=hi->holes_begin();
         hoi!=hi->holes_end();hoi++) {
      inner = *hoi;
      for (Polygon::Edge_const_iterator ei=inner.edges_begin();
           ei!=inner.edges_end();ei++) {
        Side side;
        side.first = Segment(ei->source(),ei->target());
        side.second = poly_index;
        res.push_back(side);
      }
    }
    poly_index++;
  }
  return res;
}

/** \brief Read clustering data for polygon sides from an input stream
 *
 * Data must be formatted as follows. For each side cluster:
 * + first line, 4 Cartesian coordinates of the representative segment ends 
 *   separated by spaces, followed by the number of sides grouped together,
 * + one line per side with the curvilinear abscissae of the side projection
 *   on the segment, followed by the index of the polygon the side comes from.
 * The curvilinear abscissae are relative to the segment: 0 locates at the
 * segment start, 1 at the segment end.
 */
void PolygonImporter::read_side_clusters(std::istream& input) {
  std::string line;
  while (std::getline(input,line)) {
    SideCluster sclust;
    std::vector<NT> seg_coords;
    unsigned int nb_sides;
    std::istringstream iss;
    iss.str(line);
    for (unsigned int i=0;i<4;i++) {
      CGAL::Gmpq buf;
      iss >> buf;
      seg_coords.push_back(buf);
    }
    sclust.representant = Segment(Point2(seg_coords[0],seg_coords[1]),
                                  Point2(seg_coords[2],seg_coords[3]));
    iss >> nb_sides;
    for (unsigned int i=0;i<nb_sides;i++) {
      double buf; // relative abscissae of side ends
      unsigned int pn; // polygon index
      std::getline(input,line);
      iss.clear();
      iss.str(line);
      iss >> buf;
      sclust.side_start.push_back(buf);
      iss >> buf;
      sclust.side_end.push_back(buf);
      iss >> pn;
      sclust.side_polygon.push_back(pn);
    }
    side_clusters.push_back(sclust);
  }
}

/** \brief Build an arrangement from the representative segments
 *
 * \param expand : optional extra length to be added to the segments before
 * their insertion into the arrangement. Default to zero.
 *
 * Build an arrangement with history (attribute arr) by insertion of the
 * segments listed in attribute side_clusters. Segments can be expanded at
 * both ends by an extra length.  
 */
void PolygonImporter::insert_segments(double expand) {
  arr.clear();
  for (unsigned int i=0;i<side_clusters.size();i++) {
    Segment seg = side_clusters[i].representant;
    if (expand!=0) {
      Vector vec = seg.to_vector();
      vec = (expand/sqrt(CGAL::to_double(vec.squared_length())))*vec;
      seg = Segment(seg.source()-vec,seg.target()+vec);
    }
    Arr_curve ac(seg.source(),seg.target());
    side_clusters[i].curve_ref = CGAL::insert(arr,ac);
  }
}

/** \brief Count the number of I-vertices in the arrangement.
 */
Size PolygonImporter::number_of_I_vertices() {
  Size count = 0;
  for (HistArrangement::Vertex_iterator v=arr.vertices_begin();
       v!=arr.vertices_end();v++) {
    if (v->degree()==1)
      count++;
  }
  return count;
}
/** \brief Remove I-vertices in the arrangement
 * \param within : edges ending by an I-vertex and with length larger or equal
 * to within will not be removed. Default to infinity (i.e. all edges ending 
 * an I-vertex will be removed).
 * \note Running this method does not necessarily yield an arrangement free
 * of I-vertices. Indeed, a removed edge ends by an I-vertex and a L-vertex, the
 * remaining vertex switch from L to I type.
 */
Size PolygonImporter::remove_I_vertices(double within) {
  Size count = 0;
  for (HistArrangement::Edge_iterator ei=arr.edges_begin();
       ei!=arr.edges_end();ei++) {
    if (ei->source()->degree()==1 || ei->target()->degree()==1) {
      NT len2 = CGAL::squared_distance(ei->curve().source(),
                                       ei->curve().target());
      if (sqrt(CGAL::to_double(len2))>=within)
        continue;
      arr.remove_edge(ei);
      count++;
    }
  }
  return count;
}

/** \brief Returns a HPolygon object representing an arrangement face
 * \param f : the arrangement face to be converted.
 * \param simplify : if true, only vertices joining non-aligned edges 
 * are specified during the HPolygon construction. Default: true.
 * \return a holed polygon.
 **/
HPolygon aface2poly(PolygonImporter::HistArrangement::Face_handle f, 
                    bool poly_simplify) {
  std::vector<PolygonImporter::HistArrangement::Ccb_halfedge_circulator> ccbs;
  if (f->has_outer_ccb())
    ccbs.push_back(f->outer_ccb());
  if (std::distance(f->holes_begin(),f->holes_end())>0) {
    ccbs.insert(ccbs.end(),f->holes_begin(),f->holes_end());
  }
  Polygons poly(ccbs.size());
  for (Size i=0;i<ccbs.size();i++) {
    PolygonImporter::HistArrangement::Ccb_halfedge_circulator e = ccbs[i];
    PolygonImporter::HistArrangement::Halfedge_handle e_first=e;
    do {
        poly[i].push_back(e->target()->point());
      e++;
    } while (e!=e_first);
    if (poly_simplify)
      poly[i] = simplify(poly[i]);
  }
  if (f->has_outer_ccb()) {
    Polygons::iterator hi(poly.begin());
    hi++;
    HPolygon res(poly[0],hi,poly.end());
    return(res);
  } else {
    Polygon empty;
    HPolygon res(empty,poly.begin(),poly.end());
    return(res);
  }
}
/** \brief Get arrangement faces as holed polygons*/
HPolygons PolygonImporter::get_faces(bool simplify) {
  HPolygons hpolys;
  for (HistArrangement::Face_iterator fi=arr.faces_begin();fi!=arr.faces_end();
       fi++) {
    hpolys.push_back(aface2poly(fi,simplify));
  }
  return hpolys;
}

/** \brief Curvilinear abscissa of the projection of a point onto a line
 * \param p : point to be projected.
 * \param s : segment onto the target line. The segment start is considered
 * as the origin (abscissa 0), the segment end has abscissa 1.
 */
NT curvilinear_coordinate(Point2 p,Segment s) {
  if (s.source().x()!=s.target().x()) {
    return (p.x()-s.source().x())/(s.target().x()-s.source().x());
  } else {
    return (p.y()-s.source().y())/(s.target().y()-s.source().y());
  }
}
/** \brief Length of the complement of a union of intervals within two bounds
 * \param start : first bound.
 * \param end : second bound.
 * \param intervals : union of intervals as a vector of pairs. The pair 
 * elements are the interval ends.
 */
double interval_free_length(double start,double end,
                            std::vector<std::pair<double,double> > intervals) {
  std::vector<std::pair<double, unsigned int> > bounds;
  for (unsigned int i=0;i!=intervals.size();i++) {
    if (intervals[i].first>end || intervals[i].second<start) continue;
    bounds.push_back(std::make_pair(fmax(intervals[i].first,start),i));
    bounds.push_back(std::make_pair(fmin(intervals[i].second,end),i));
  }
  if (bounds.empty()) return end-start;
  std::sort(bounds.begin(),bounds.end());
  std::set<unsigned int> covering;
  double res = 0.0;
  std::vector<std::pair<double,unsigned int> >::iterator bi = bounds.begin();
  double last_out = start;
  for (;bi!=bounds.end();bi++) {
    if (covering.empty()) {
      res += bi->first-last_out;
    }
    if (covering.erase(bi->second)==0) {
      covering.insert(bi->second);
    }
    if (covering.empty()) {
      last_out = bi->first;
    }
  }
  bi--;
  res += end-bi->first;
  return res;
}
/** \brief Find the side cluster represented by an arrangement segment
 * \param ch : the arrangement segment as a handle.
 * \return an iterator to the corresponding side cluster.
 */
PolygonImporter::SideClusters::iterator 
PolygonImporter::fetch_side_cluster(Curve_handle ch) {
  for (SideClusters::iterator sci=side_clusters.begin();
       sci!=side_clusters.end();sci++) {
    if (sci->curve_ref==ch)
      return sci;
  }
}

/** \brief Compute polygon votes at halfedge level
 * \param e : halfedge of attribute arr. It must come from the insertion of a
 * representative segment stored in side_clusters.
 * \return votes of sides for polygons.
 *
 * Sides allowed to vote are those represented by the segment whose insertion 
 * has generated the halfedge e. Sides opposite to e are excluded from the vote.
 * A side votes for the polygon it comes from and its vote is weighted by the
 * length of its projection on the line supporting the segment. Thus the total 
 * vote for a given polygon is a sum of lengths. A null vote is also computed
 * (vote for fictitious polygon 0): it is equal to the length of the halfedge
 * free of side projection.
 *
 * The total sum of votes is at least equal to the halfedge length. It may
 * exceed the halfedge length as soon as some side projections overlap.
 */
PolygonImporter::PolygonVotes 
PolygonImporter::polygon_sides(HistArrangement::Halfedge_handle e) {
  PolygonVotes res;
  HistArrangement::Originating_curve_iterator ci = arr.originating_curves_begin(e);
  Curve_handle ch = ci;
  SideClusters::iterator sci = fetch_side_cluster(ch);
  double e_start = CGAL::to_double(curvilinear_coordinate(e->source()->point(),
                                                          sci->representant));
  double e_end = CGAL::to_double(curvilinear_coordinate(e->target()->point(),
                                                          sci->representant));
  double e_right, e_left;
  if (e_start<e_end) {
    e_right = e_start;
    e_left = e_end;
  } else {
    e_right = e_end;
    e_left = e_start;
  }
  std::vector<std::pair<double,double> > side_ends;
  for (unsigned int i=0;i!=sci->side_start.size();i++) {
    if ((e_end-e_start)*(sci->side_end[i]-sci->side_start[i])<=0) {
      continue;
    }
    double side_right, side_left;
    if (sci->side_start[i]<sci->side_end[i]) {
      side_right = sci->side_start[i];
      side_left = sci->side_end[i];
    } else {
      side_right = sci->side_end[i];
      side_left = sci->side_start[i];
    }
    double inter_right = fmax(e_right,side_right);
    double inter_left = fmin(e_left,side_left);
    if (inter_left<=inter_right) {
      continue;
    }
    if (res.find(sci->side_polygon[i])==res.end()) {
      res[sci->side_polygon[i]] = inter_left-inter_right;
    } else {
      res[sci->side_polygon[i]] += inter_left-inter_right;
    }
    side_ends.push_back(std::make_pair(side_right,side_left));
  }
  double ifl = interval_free_length(e_right,e_left,side_ends);
  if (ifl>0)  
    res[0] = ifl;
  double seg_length = sqrt(CGAL::to_double(sci->representant.squared_length()));
  for (PolygonVotes::iterator pvi=res.begin();pvi!=res.end();pvi++) {
    pvi->second *= seg_length;
  } 
  return res;
}

/** \brief Compute polygon votes at face level
 * \param f : face of attribute arr.
 * \return votes of sides for polygons.
 *
 * Votes are collected for each halfedge bounding the face (see
 * polygon_sides(HistArrangement::Halfedge_handle)). Results are summed
 * and normalized by the halfedge perimeter.
 */
PolygonImporter::PolygonVotes 
PolygonImporter::polygon_sides(HistArrangement::Face_handle f) {
  double perimeter = 0.0;
  PolygonVotes res;
  std::vector<HistArrangement::Ccb_halfedge_circulator> ccbs;
  if (f->has_outer_ccb())
    ccbs.push_back(f->outer_ccb());
  for (HistArrangement::Hole_iterator hi=f->holes_begin();hi!=f->holes_end();
       hi++)
    ccbs.push_back(*hi);
  for (std::vector<HistArrangement::Ccb_halfedge_circulator>::iterator
         ccbi=ccbs.begin();ccbi!=ccbs.end();ccbi++) {
    HistArrangement::Ccb_halfedge_circulator e = *ccbi, done = *ccbi;
    do {
      PolygonVotes e_votes = polygon_sides(e);
      for (PolygonVotes::iterator pvi=e_votes.begin();pvi!=e_votes.end();
           pvi++) {
        if (res.find(pvi->first)==res.end()) {
          res[pvi->first] = pvi->second;
        } else {
          res[pvi->first] += pvi->second;
        }
      }
      NT len2 =  CGAL::squared_distance(e->source()->point(),
                                        e->target()->point());
      perimeter += sqrt(CGAL::to_double(len2));
      e++;
    } while (e!=done);
  }
  for (PolygonVotes::iterator pvi=res.begin();pvi!=res.end();pvi++) {
    pvi->second /= perimeter;
  }
  return res;
}

/** \brief Compute which polygon wins a vote*/
PolygonImporter::PolygonVote vote_winner(PolygonImporter::PolygonVotes pvs) {
  PolygonImporter::PolygonVote res;
  res.first = pvs.size();
  res.second = 0.0;
  for (PolygonImporter::PolygonVotes::iterator pvi=pvs.begin();
       pvi!=pvs.end();pvi++) {
    if (pvi->second>res.second) {
      res.first = pvi->first;
      res.second = pvi->second;
    }
  }
  return res;
}
/** \brief Compute the polygon mostly associated with a face
 * \param f : face of arr.
 * \return the index of the polygon mostly associated with the face together
 * with its score.
 *
 * Computations are based on polygon votes as returned by
 * polygon_sides(HistArrangement::Face_handle). Note that polygon
 * indices start from 1 and that a vote mostly for polygon 0 means
 * that no input polygon was found to match the face. The face may be
 * outside the tessellated domain (unbounded external face or a domain
 * hole).
 */
PolygonImporter::PolygonVote PolygonImporter::elected_polygon(HistArrangement::Face_handle f) {
  PolygonVotes pvs = polygon_sides(f);
  PolygonVote res = vote_winner(pvs);
  return res;
}

std::vector<PolygonImporter::PolygonVote>
PolygonImporter::elected_polygons(HistArrangement::Ccb_halfedge_circulator e) {
  if (e==NULL) {
    throw std::domain_error("empty circulator provided to "
                            "PolygonImporter::elected_polygons");
  }
  std::vector<PolygonVote> res;
  HistArrangement::Ccb_halfedge_circulator done = e;
  do {
    PolygonVote poly_edge = vote_winner(polygon_sides(e));
    NT len2 =  CGAL::squared_distance(e->source()->point(),
                                      e->target()->point());
    poly_edge.second /= sqrt(CGAL::to_double(len2));
    res.push_back(poly_edge);
    e++;
  } while (e!=done);
  return res;
}

/** \brief Area of a holed polygon */
NT area(HPolygon hp) {
  Polygons p = boundaries(hp);
  NT res = 0; 
  for (Polygons::iterator pi=p.begin();pi!=p.end();pi++) {
    res += CGAL::polygon_area_2(pi->vertices_begin(),pi->vertices_end(),
                                Polygon::Traits());
  }
  return res;
}
/** \brief Sample points along a line segment
 */
Points digitize(Segment seg, NT step_length) {
  Points res;
  if (seg.source()==seg.target()) {
    res.push_back(seg.source());
    return res;
  }
  Vector vec(seg.source(),seg.target());
  NT len = sqrt(CGAL::to_double(vec.squared_length()));
  vec = vec /len;
  for (NT delta=0;delta<=len;delta+=step_length) {
    res.push_back(seg.source()+delta*vec);
  }
  return res;
}
/** \brief Sample points along line segments */
Points digitize(std::vector<Segment> segs,NT step_length) {
  Points res;
  for (std::vector<Segment>::iterator si=segs.begin();si!=segs.end();si++) {
    Points buf = digitize(*si,step_length);
    res.insert(res.end(),buf.begin(),buf.end());
  }
  return res;
}
  
/** \brief Sample points along the boundaries of a holed polygon
 * \param hpoly : the input holed polygon.
 * \param step_length: distance between consecutive sampling points.
 * \return set of sampling points.
 */
Points digitize(HPolygon hpoly,NT step_length) {
  std::vector<Segment> segs;
  if (!hpoly.is_unbounded()) {
    for (Polygon::Edge_const_iterator ei=hpoly.outer_boundary().edges_begin();
         ei!=hpoly.outer_boundary().edges_end();ei++) {
      segs.push_back(*ei);
    }
  }
  for (HPolygon::Hole_const_iterator hi=hpoly.holes_begin();
       hi!=hpoly.holes_end();hi++) {
    for (Polygon::Edge_const_iterator ei=hi->edges_begin();
         ei!=hi->edges_end();ei++) {
      segs.push_back(*ei);
    }
  }
  Points res = digitize(segs,step_length);
  return res;
}
/** \brief Compute the squared Hausdorff distance between two point sets
 */
NT squared_Hausdorff_distance(Points& a,Points& b) {
  std::vector<NT> da(a.size(),-1), db(b.size(),-1);
  for (Size i=0;i<a.size();i++) {
    for (Size j=0;j<b.size();j++) {
      NT dist = CGAL::squared_distance(a[i],b[j]);
      if (dist<da[i] || da[i]<0) da[i] = dist;
      if (dist<db[j] || db[j]<0) db[j] = dist;
    }
  }
  std::vector<NT> dab(da);
  dab.insert(dab.end(),db.begin(),db.end());
  std::vector<NT>::iterator max_it = std::max_element(dab.begin(),dab.end());
  return *max_it;
}
NT squared_Hausdorff_distance_vd(Points& a,Points& b) {
  typedef CGAL::Delaunay_triangulation_2<Kernel> DT;
  typedef CGAL::Delaunay_triangulation_adaptation_traits_2<DT> AT;
  typedef CGAL::Identity_policy_2<DT,AT> AP;
  typedef CGAL::Voronoi_diagram_2<DT,AT,AP> VD;
 
  VD vda(a.begin(),a.end());
  VD vdb(b.begin(),b.end());
  NT res = 0;
  for (Points::iterator pi=a.begin();pi!=a.end();pi++) {
    VD::Locate_result lr = vdb.locate(*pi);
    if (VD::Face_handle* f = boost::get<VD::Face_handle>(&lr)) {
      Point2 nearest = (*f)->dual()->point();
      NT dist2 = CGAL::squared_distance(*pi,nearest);
      if (dist2>res)
        res = dist2;
    } else {
      throw std::domain_error("squared_Hausdorff_distance_v,"
                              " case not supported");
    }
  }
  for (Points::iterator pi=b.begin();pi!=b.end();pi++) {
    VD::Locate_result lr = vda.locate(*pi);
    if (VD::Face_handle* f = boost::get<VD::Face_handle>(&lr)) {
      Point2 nearest = (*f)->dual()->point();
      NT dist2 = CGAL::squared_distance(*pi,nearest);
      if (dist2>res)
        res = dist2;
    } else {
      throw std::domain_error("squared_Hausdorff_distance_v,"
                              " case not supported");
    }
  }
  return res;
}
NT squared_Hausdorff_distance_dt(Points& a, Points& b) {
  typedef CGAL::Delaunay_triangulation_2<Kernel> DT;
  DT dta, dtb;
  dta.insert(a.begin(),a.end());
  dtb.insert(b.begin(),b.end());
  NT res = 0;
  for (Points::iterator pi=a.begin();pi!=a.end();pi++) {
    DT::Vertex_handle v_nearest = dtb.nearest_vertex(*pi);
    NT dist2 = CGAL::squared_distance(*pi,v_nearest->point());
    if (dist2>res)
      res = dist2;
  }
  for (Points::iterator pi=b.begin();pi!=b.end();pi++) {
    DT::Vertex_handle v_nearest = dta.nearest_vertex(*pi);
    NT dist2 = CGAL::squared_distance(*pi,v_nearest->point());
    if (dist2>res)
      res = dist2;
  }
  return res;
}
NT squared_Hausdorff_distance_dtbis(Points& a,Points& b) {
  typedef CGAL::Point_set_2<Kernel> PS;
  PS psa(a.begin(),a.end()), psb(b.begin(),b.end());
  NT res = 0;
  for (Points::iterator pi=a.begin();pi!=a.end();pi++) {
    PS::Vertex_handle v_nearest = psb.nearest_neighbor(*pi);
    NT dist2 = CGAL::squared_distance(*pi,v_nearest->point());
    if (dist2>res)
      res = dist2;
  }
  for (Points::iterator pi=b.begin();pi!=b.end();pi++) {
    PS::Vertex_handle v_nearest = psa.nearest_neighbor(*pi);
    NT dist2 = CGAL::squared_distance(*pi,v_nearest->point());
    if (dist2>res)
      res = dist2;
  }
  return res;
}
NT squared_Hausdorff_distance_os(PolygonImporter::OS::Tree* a,
                                 PolygonImporter::OS::Tree* b,
                                 NT eps) {
  NT res = 0;
  for (PolygonImporter::OS::Tree::iterator pi=a->begin();pi!=a->end();pi++) {
    PolygonImporter::OS search(*b,*pi,1,eps);
    NT dist2 = search.begin()->second;
    if (dist2>res)
      res = dist2;
  }
  for (PolygonImporter::OS::Tree::iterator pi=b->begin();pi!=b->end();pi++) {
    PolygonImporter::OS search(*a,*pi,1,eps);
    NT dist2 = search.begin()->second;
    if (dist2>res)
      res = dist2;
  }
  return res;
}
NT squared_Hausdorff_distance(std::vector<Segment>& P,std::vector<Segment>& Q) {
  typedef CGAL::Segment_Delaunay_graph_traits_2<Kernel> SDGTraits;
  typedef CGAL::Segment_Delaunay_graph_2<SDGTraits> SDG;
  typedef CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG> SDGAdaptTraits;
  typedef CGAL::Segment_Delaunay_graph_degeneracy_removal_policy_2<SDG> SDGPol;
  typedef CGAL::Voronoi_diagram_2<SDG,SDGAdaptTraits,SDGPol> SVD;
  SVD VorP;
  SDGAdaptTraits::Point_2 a(0,0), b(1,0);
  SDGAdaptTraits::Site_2 toto;
  toto.construct_site_2(a,b);
  VorP.insert(toto);
  return 0;
}
/** \brief Compare a face to an input polygon
 * \param f : face of arr.
 * \param an index of an input polygon (starting from 1).
 * \param eps : step length to be used for digitization.
 * border.
 * \return the squared discrete Hausdorff distance between the face and the
 * polygon borders.
 */
NT PolygonImporter::compare(Size fidx,Size pidx,NT eps) {
  NT res = squared_Hausdorff_distance_os(polygon_trees[pidx],face_trees[fidx],
                                         eps);
  return res;
}

void PolygonImporter::preprocess_polygons_4compare(NT eps) {
  polygon_trees.clear();
  for (HPolygons::iterator hpi = input_polygons.begin();
       hpi != input_polygons.end();hpi++) {
    Points sample = digitize(*hpi,eps);
    OS::Tree* tree = new OS::Tree(sample.begin(),sample.end());
    polygon_trees.push_back(tree);
  }
}
void PolygonImporter::preprocess_faces_4compare(NT eps) {
  face_trees.clear();
  for (HistArrangement::Face_iterator fi = arr.faces_begin();
       fi != arr.faces_end(); fi++) {
    HPolygon f_poly = aface2poly(fi,true);
    Points sample = digitize(f_poly,eps);
    OS::Tree* tree = new OS::Tree(sample.begin(),sample.end());
    face_trees.push_back(tree);
  }
}
/** \brief Assess globally the fit polygon/face
 */
FMatrix PolygonImporter::goodness_of_fit(NT eps) {
  std::vector<FVector> res;
  NT dmax = 0;
  std::vector<Size> missing;
  for (Size i=0;i<input_polygons.size();i++)
    missing.push_back(i);
  std::vector<Size>::iterator missing_begin = missing.begin(), 
    missing_end = missing.end();
  Size face_counter = 0;
  for (HistArrangement::Face_iterator fi=arr.faces_begin();
       fi!=arr.faces_end();fi++) {
    Size face_index = std::distance(arr.faces_begin(),fi);
    PolygonVotes pv = polygon_sides(fi);
    for (PolygonVotes::iterator pvi=pv.begin();pvi!=pv.end();pvi++) {
      if (pvi->first>0) {
        NT dist2 = compare(face_index,pvi->first-1,eps);
        FVector col(4);
        col[0] = (double) std::distance(arr.faces_begin(),fi);
        col[1] = (double) pvi->first-1;
        col[2] = pvi->second;
        col[3] = CGAL::to_double(dist2);
        res.push_back(col);
        missing_end = std::remove(missing_begin,missing_end,pvi->first-1);
      }
    }
    face_counter++;
  }
  if (std::distance(missing_begin,missing_end)>0) {
    std::clog << "Warning: there are " << std::distance(missing_begin,missing_end) << " orphan polygons" << std::endl;
  }
  FMatrix matres(res.begin(),res.end());
  FMatrix mtransp = LA::transpose(matres);
  return mtransp;
}
