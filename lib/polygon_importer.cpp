#include "ttessel.h"
double infinity_double = std::numeric_limits<double>::infinity();
NT infinity_NT = std::numeric_limits<NT>::infinity();
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

void PolygonImporter::insert_segments(double expand) {
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

Size PolygonImporter::number_of_I_vertices() {
  Size count = 0;
  for (HistArrangement::Vertex_iterator v=arr.vertices_begin();
       v!=arr.vertices_end();v++) {
    if (v->degree()==1)
      count++;
  }
  return count;
}
Size PolygonImporter::remove_I_vertices(double within) {
  Size count = 0;
  for (HistArrangement::Vertex_iterator v=arr.vertices_begin();
       v!=arr.vertices_end();v++) {
    if (v->degree()==1) {
      HistArrangement::Halfedge_around_vertex_circulator 
        e = v->incident_halfedges();
      NT len2 = CGAL::squared_distance(e->curve().source(),
                                       e->curve().target());
      if (sqrt(CGAL::to_double(len2))>=within)
        continue;
      arr.remove_edge(e);
      count++;
    }
  }
  return count;
}

NT curvilinear_coordinate(Point2 p,Segment s) {
  if (s.source().x()!=s.target().x()) {
    return (p.x()-s.source().x())/(s.target().x()-s.source().x());
  } else {
    return (p.y()-s.source().y())/(s.target().y()-s.source().y());
  }
}

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
    
PolygonImporter::SideClusters::iterator 
PolygonImporter::fetch_side_cluster(Curve_handle ch) {
  for (SideClusters::iterator sci=side_clusters.begin();
       sci!=side_clusters.end();sci++) {
    if (sci->curve_ref==ch)
      return sci;
  }
}

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

PolygonImporter::PolygonVote PolygonImporter::elected_polygon(HistArrangement::Face_handle f) {
  PolygonVote res;
  PolygonVotes pv = polygon_sides(f);
  res.first = pv.size();
  res.second = 0.0;
  for (PolygonImporter::PolygonVotes::iterator pvi=pv.begin();
       pvi!=pv.end();pvi++) {
    if (pvi->second>res.second) {
      res.first = pvi->first;
      res.second = pvi->second;
    }
  }
  return res;
}

NT area(HPolygon hp) {
  Polygons p = boundaries(hp);
  NT res = 0; 
  for (Polygons::iterator pi=p.begin();pi!=p.end();pi++) {
    res += CGAL::polygon_area_2(pi->vertices_begin(),pi->vertices_end(),
                                Polygon::Traits());
  }
  return res;
}
Points digitize(Points polyline,NT step_length) {
  Points res;
  for (Points::iterator src=polyline.begin();true;src++) {
    Points::iterator end = src+1;
    if (end==polyline.end()) break;
    Vector vec(*src,*end);
    NT len = sqrt(CGAL::to_double(vec.squared_length()));
    vec = vec /len;
    for (NT delta=0;delta<=len;delta+=step_length) {
      res.push_back(*src+delta*vec);
    }
  }
  return res;
}
NT squared_Hausdorff_distance(Points& a,Points& b) {
  std::vector<NT> da(a.size()), db(b.size());
  for (Size i=0;i<a.size();i++) {
    for (Size j=0;j<b.size();j++) {
      NT dist = CGAL::squared_distance(a[i],b[j]);
      if (dist<da[i]) da[i] = dist;
      if (dist<db[j]) db[j] = dist;
    }
  }
  std::vector<NT> dab(da);
  dab.insert(dab.end(),db.begin(),db.end());
  return *std::max_element(dab.begin(),dab.end());
}
NT PolygonImporter::compare(HistArrangement::Face_handle f,unsigned int pidx,
                            NT eps) {
  Points face_sample;
  Points poly_sample;
  HistArrangement::Ccb_halfedge_circulator ccb_begin, ccb_end;
  if (f->has_outer_ccb()) {
    Points vertices;
    HistArrangement::Ccb_halfedge_circulator e, done;
    e = f->outer_ccb();
    done = e;
    do {
      vertices.push_back(e->source()->point());
      e++;
    } while (e!=done);
    vertices.push_back(e->source()->point());
    Points buf = digitize(vertices,eps);
    face_sample.insert(face_sample.end(),buf.begin(),buf.end());
  }
  for (HistArrangement::Hole_iterator hi=f->holes_begin();hi!=f->holes_end();
       hi++) {
    Points vertices;
    HistArrangement::Ccb_halfedge_circulator e, done;
    e = *hi;
    done = e;
    do {
      vertices.push_back(e->source()->point());
      e++;
    } while (e!=done);
    vertices.push_back(e->source()->point());
    Points buf = digitize(vertices,eps);
    face_sample.insert(face_sample.end(),buf.begin(),buf.end());
  }
  HPolygon ref_hp = input_polygons[pidx-1];
  if (!ref_hp.is_unbounded()) {
    Points vertices(ref_hp.outer_boundary().vertices_begin(),
                    ref_hp.outer_boundary().vertices_end());
    vertices.push_back(vertices[0]);
    Points buf = digitize(vertices,eps);
    poly_sample.insert(poly_sample.end(),buf.begin(),buf.end());
  }
  for (HPolygon::Hole_const_iterator hi=ref_hp.holes_begin();
       hi!=ref_hp.holes_end();hi++) {
    Points vertices(hi->vertices_begin(),hi->vertices_end());
    vertices.push_back(vertices[0]);
    Points buf = digitize(vertices,eps);
    poly_sample.insert(poly_sample.end(),buf.begin(),buf.end());
  }
  NT res = squared_Hausdorff_distance(face_sample,poly_sample);
  return res;
}
