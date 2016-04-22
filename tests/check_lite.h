#include <fstream>
#include <vector>
#include <algorithm>
#include <exception>
#include "ttessel.h"

Point2 upper_left_vertex(const Polygon&);
bool llt(const Polygon&,const Polygon&);
bool operator==(const HPolygon&,const HPolygon&);
HPolygons get_faces(TTessel&);
HPolygons set_diff(HPolygons&, HPolygons&);
std::ostream& operator<<(std::ostream&,const HPolygon&);
std::ostream& operator<<(std::ostream&,const HPolygons&);
std::istream& operator>>(std::istream&,const HPolygons&);
TTessel::Halfedge_handle find_halfedge(TTessel&,Point2,Point2);
HPolygons holed_polygons();
HPolygons holed_polygons_for_flips();
bool check_point_closeness(Point2, Point2, double tolerance=1e-6);
bool check_polygon_closeness(Polygon, Polygon, double tolerance=1e-6);
inline bool check_polygon_closeness_default(Polygon a, Polygon b) {
  return check_polygon_closeness(a,b);};
bool check_holed_polygon_closeness(HPolygon, HPolygon);
bool match(HPolygon&,HPolygons&);
bool is_point_outside(Point2,HPolygon);
bool is_point_outside(Point2,HPolygons);
template <class ModifType>
bool prediction_is_right(TTessel&, ModifType&,HPolygons&,
			 HPolygons&,HPolygons&,HPolygons&,std::ostringstream&);
void process_wrong_predictions(HPolygons,HPolygons,HPolygons,HPolygons,
			       const char[],const char[],const char[],
			       const char[],std::ostringstream&);
void check_predictions_4_random_smf(unsigned,unsigned,unsigned,HPolygons&);
