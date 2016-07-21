/* Line Tessellation (LiTe) library
   |||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
   |||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

/** \file */

/******************************************************************************/
/*             LINE BASED TESSELLATION CLASS - DECLARATIONS                   */
/******************************************************************************/

/******************************************************************************/
/*                               INCLUDES                                     */
/******************************************************************************/

// Define shorter names to please linker (g++)
//#include "short_names.h"
//#include <cotime>
#include <iostream>
#include <exception>
#include <boost/unordered_map.hpp>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
/* Calcul exact pour +,-,*,/. On peut sauvegarder la tessellation dans
un fichier texte sans perte de précision. Par contre, la racine carrée
n'est implémentée.*/ 
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/double.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Random.h>
// Two following headers required for clip_segment_by_convex_polygon
#include <CGAL/Extended_cartesian.h>
#include <CGAL/Nef_polyhedron_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Arrangement_with_history_2.h>

/******************************************************************************/
/*          TYPE DEFINITIONS AND FORWARD CLASS DECLARATIONS                   */
/******************************************************************************/
/** \typedef NT 
 * \brief The CGAL numeric type used for exact computations in LiTe*/
typedef CGAL::Lazy_exact_nt<CGAL::Gmpq>                         NT;
/** \typedef Kernel
 * \brief The CGAL type of coordinates used in LiTe*/
typedef CGAL::Cartesian<NT>                                     Kernel;
/** \typedef Traits
 * \brief The CGAL arrangement traits class used for representing 
 * tessellations*/
typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits;
/** \typedef Vector
 * \brief The CGAL class used for representing vectors*/
typedef CGAL::Vector_2<Kernel>                                  Vector;
/** \typedef Point2
 * \brief CGAL class representing a point from a tessellation*/
typedef Traits::Point_2                                         Point2;
/** \typedef Segment
 * \brief CGAL class representing a line segment from a tessellation*/
typedef Traits::Segment_2                                       Segment;
/** \typedef Curve
 * \brief CGAL class representing a curve (line segment) from 
 * a tessellation*/
typedef Traits::X_monotone_curve_2                              Curve;
/** typedef \brief CGAL class representing a rectangle with horizontal
 * and vertical sides*/
typedef CGAL::Iso_rectangle_2<Kernel>                           Rectangle;
/** \typedef Line
 * \brief CGAL class representing an infinite line*/
typedef CGAL::Line_2<Kernel>                                    Line;
/** \typedef Rayon
 * \brief CGAL class representing a half-line*/
typedef CGAL::Ray_2<Kernel>                                     Rayon;
/** \typedef Direction
 * \brief CGAL class representing a direction in the plane*/
typedef CGAL::Direction_2<Kernel>                               Direction;
/** \typedef Transformation
 * \brief CGAL class representing an affine transformation*/
typedef CGAL::Aff_transformation_2<Kernel>                      Transformation;
/** \typedef Polygon
 * \brief CGAL class representing a polygon*/
typedef CGAL::Polygon_2<Kernel>                                 Polygon;
/** \typedef PECirc
 * \brief CGAL Polygon Edge circulator*/
typedef Polygon::Edge_const_circulator                          PECirc;
/** \typedef Polygons
 * \brief List of Polygons*/
typedef std::vector<Polygon>                                    Polygons;
/** \typedef HPolygon
 * \brief CGAL class representing a polygon with holes*/
typedef CGAL::Polygon_with_holes_2<Kernel>                      HPolygon;
/** \typedef Points
 * \brief Class representing a series of points*/
typedef std::vector<Point2>                                     Points;
/** \typedef LA
 * \brief Class used for basic linear algebra
 * 
 * This class is used in particular for some computations on 
 * matrices.*/
typedef CGAL::Linear_algebraCd<double>                          LA;
/** \typedef FVector
 * \brief Class used for representing a vector (in the linear algebra
 * sense)*/
typedef LA::Vector                                              FVector;
/** \typedef FMatrix
 * \brief Class used for representing a matrix */
typedef LA::Matrix                                              FMatrix;

class LineTes_halfedge; // forward declaration

/** \typedef Dcel
 * \brief Class used as a template parameter for the Arrangement class
 * 
 * Use of non-standard halfedges (class LineTes_halfedge) and 
 * non-standard faces*/
typedef CGAL::Arr_face_extended_dcel<Traits,bool,
  CGAL::Arr_vertex_base<Point2>,
  LineTes_halfedge,CGAL::Arr_face_base>                          Dcel;
/** \typedef Arrangement
 * \brief Class representing a tessellation*/
typedef CGAL::Arrangement_2<Traits,Dcel>                        Arrangement;
/** \typedef Size
 * \brief Number type used for indexing tessellation features*/
typedef Dcel::Size                                              Size;

/** \brief LiTe random generator
 *
 * LiTe implements several stochastic algorihtms where this 
 * generator is used.*/
extern CGAL::Random *rnd; // variable globale

/** \brief Infinite value for doubles
 *
 * Any double is smaller than that value.
 */
extern double infinity_double;

// definition that prevent type expansion in the documentation generated by Rcpp
/** \brief List of polygons with holes
 *
 * Just a vector of HPolygon objects.
 */
class HPolygons : public std::vector<HPolygon> {};

/******************************************************************************/
/*                  DÉFINITION DE LA CLASSE LineTes                           */
/******************************************************************************/

/** \brief Line tessellation of a bounded polygonal domain
 *
 * A line tessellation is defined as a planar tessellation whose
 * network of edges is a subset of a union of lines. The lines are
 * supposed to be in general position. For instance, line
 * configurations with more than two lines crossing at the same point
 * are excluded. 
 *
 * Line tessellations can be considered as patterns of line
 * segments. A segment is defined as a maximal subset of edges that
 * are aligned and connected.  In a line tessellation, two distinct
 * segments cannot be supported by the same line.  The %LineTes class
 * inherits from the Arrangement_2 class provided by CGAL.  Compared
 * to its parent class, the %LineTes class is of interest for
 * tessellations with segments that consist of more than one
 * edge. Segments of a line tessellation are represented as Seg
 * objects.
*/

class LineTes : public Arrangement {

 public:
  
  /** \brief Segment in a line tessellation
   *
   * In a line tessellation represented by a LineTes object, a segment
   * is defined as a maximal subset of edges that are aligned and
   * connected. A %Segment object is mainly a Halfedge_handle referring
   * to an edge lying on the segment.
   */
  class Seg {
  public:
    /** \name Initialization */
    /** \{ */
    Seg();
    Seg(Halfedge_handle);
    /** \} */
    /** \name Access */
    /** \{ */
    Halfedge_handle         halfedges_start();
    Halfedge_handle         halfedges_end();
    Point2                  pointSource();
    Point2                  pointTarget();
    Points                  list_of_points();
    /** \brief Return the halfedge handle identifying the segment*/
    inline Halfedge_handle  get_halfedge_handle(){return e;}
    void                    set_halfedge_handle(Halfedge_handle);
    /** \} */
    /** \name Features */
    /** \{ */
    Size                    number_of_edges();
    bool                    number_of_edges_is_greater_than(unsigned int);
    /** \} */
    /** \name Input/output */
    /** \{ */
    void                    print(bool endsOnly=false);
    /** \} */
  private:
    Halfedge_handle e;
    Size            or_idx;
  };

  /** \typedef Seg_handle
   * \brief  Pointer to a Seg object*/
  typedef Seg* Seg_handle;
 
  /** \brief List of Seg objects
   *
   * This class is intended mainly for representing all segments in a
   * LineTes object. For representing only a subset of segments, use
   * thee Seg_sublist class instead.
   */
  class Seg_list : public std::vector<Seg_handle> {
  public:
    /** \brief Add a segment at the end of the segment list */
    inline void add(Seg_handle s)     { push_back(s);}
    /** \brief Delete a Seg object and remove it from the segment list. 
     *
     * \param s : handle of the segment to be suppressed.
     *
     * Yield a warning and return false if the input segment is not
     * found in the list. */
    inline bool suppress(Seg_handle s){
      iterator ou_s = std::find(begin(),end(),s);
      if (ou_s==end()) {
	std::clog << "LineTes::Seg_list::suppress : item to suppress not found"
		  << std::endl;
	return false;
      }
      erase(ou_s);
      delete s;
      return true;
    }
  };
  /** \brief Sublist of Seg objects
   *
   * To be used for storing e.g. a sublist of segments of a given type.
   */
 class Seg_sublist : public std::vector<Seg_handle> {
 public:
   /** \brief Add a segment at the end of the segment list */
   inline void add(Seg_handle s)      { push_back(s);}
   /** \brief Remove a segment from the segment list. 
    *
    * \param s : handle of the segment to be suppressed.
    *
    * Yield a warning and return false if the input segment is not
    * found in the list. */
   inline bool suppress(Seg_handle s ){
     iterator ou_s = std::find(begin(),end(),s);
     if (ou_s==end()) {
       std::clog << "LineTes::Seg_sublist::suppress : item to suppress not found"
		 << std::endl;
       return false;
     }
     erase(ou_s);
     return true;
   }
 };
 /** \typedef Seg_list_iterator
  * \brief Iterator for accessing segments in a list*/  
 typedef Seg_list::iterator      Seg_list_iterator;
 /** \typedef Seg_sublist_iterator
  * \brief Iterator for accessing segments in a sublist of segments*/
 typedef Seg_sublist::iterator   Seg_sublist_iterator;

 /**\name Initialization  */ 
 /**\{*/
  LineTes();
  virtual void                  insert_window(Rectangle);
  virtual void                  insert_window(Polygon&);
  virtual void                  insert_window(HPolygons&);
  /**\}*/
  /**\name Modification */
  /**\{ */
  Halfedge_handle               split_from_edge(Halfedge_handle, Point2, 
					             Halfedge_handle, Point2);
  Halfedge_handle               split_from_vertex(Vertex_handle, Halfedge_handle,
					          Point2);
  Face_handle                   suppress_edge(Halfedge_handle);
  void                          clear(bool remove_window=true);
  Seg_handle                    insert_segment(Segment);
  void                          remove_ivertices(Size imax=0,
						 bool verbose=false);
  void                          remove_lvertices(Size imax=0,
						 bool verbose=false);
  /**\}*/
  /** \name Access */
  /** \{ */
  /** \brief Return an iterator pointing to the beginning of the list of
   * tessellation segments
   */
  inline Seg_list_iterator      segments_begin(){return all_segments.begin();}
  /** \brief Return an iterator pointing to the end of the list of
   * tessellation segments
   */
  inline Seg_list_iterator      segments_end(){return all_segments.end();}
  /** \brief Return the domain that is tessellated 
   */
  inline HPolygons              get_window(){return window;} 
  /** \} */
  /** \name Features */
  /** \{ */
  /** \brief Return the total number of segments
   *
   * Boundary and internal segments are included
   */
  inline Size                   number_of_segments(){return all_segments.size();}
  Size                          number_of_internal_segments();
  Size                          number_of_window_edges();
  double                        get_window_perimeter();
  /** \} */
  /** \name Input/output */
  /** \{ */
  void                          print_all_segments(bool endsOnly=false);
  void                          write(std::ostream&);
  void                          read(std::istream&);
  /** \} */
  /** \name Checking and testing */
  /** \{*/
  int                           is_on_boundary(Halfedge_handle);
  /** \brief Test whether a segment lies on the domain boundary
   *
   * \param s : a handle to the segment to be tested
   */
  inline bool                   is_on_boundary(Seg_handle s){return is_on_boundary(s->get_halfedge_handle())>0;}
  bool                          is_valid(bool verbose=false);
  bool                          check_all_segments(bool verbose=false);
  bool                          check_segment(Seg_handle s, bool verbose=false);
  bool                          check_all_halfedges(bool verbose=false);
  bool                          check_halfedge(LineTes::Halfedge_handle e,
					       bool verbose=false);
  bool                          halfedge_exists(LineTes::Halfedge_handle);
  bool                          check_halfedge_neighbours(LineTes::Halfedge_handle e,
						       bool verbose=false);
  bool                          check_halfedge_prev_hf(LineTes::Halfedge_handle e,
						       bool verbose=false);
  bool                          check_halfedge_next_hf(LineTes::Halfedge_handle e,
						       bool verbose=false);
  bool                          check_halfedge_segment(LineTes::Halfedge_handle e,
						       bool verbose=false);
  bool                          check_halfedge_dir(LineTes::Halfedge_handle e,
						       bool verbose=false);
  bool                          is_a_T_tessellation(bool verbose=false,
						    std::ostream& out=std::clog);
  /** \}*/
private:
  HPolygons                     window;
  Seg_list                      all_segments;
  Halfedge_handle               split_edge(Halfedge_handle,X_monotone_curve_2,
				           X_monotone_curve_2);
  Halfedge_handle               merge_edge(Halfedge_handle, Halfedge_handle,
			                   X_monotone_curve_2); 
  Halfedge_handle               add_edge(X_monotone_curve_2,Vertex_handle, 
				         Vertex_handle); 
  Face_handle                   remove_edge(Halfedge_handle);
};

/* Variables globales declarées avant leur utilisation */

/** \brief Null value for a LineTes halfedge handle*/
extern LineTes::Halfedge_handle NULL_HALFEDGE_HANDLE;
/** \brief Null value for a LineTes face handle*/
extern LineTes::Face_handle NULL_FACE_HANDLE;
/** \brief Null value for a LineTes segment handle*/
extern LineTes::Seg_handle NULL_SEG_HANDLE;

/******************************************************************************/
/*            DEFINITION OF THE CLASS LineTes_halfedge                        */
/******************************************************************************/

/** \brief Halfedge class used by the LineTes class
 *
 * The LineTes_halfedge class derives from the CGAL class
 * Arr_halfedge_base. Several attributes related to segments have been
 * added.
 *
 * A segment represented by a LineTes::Seg object is identified by one
 * of its halfedge. It is therefore oriented and its orientation is
 * determined by the orientation of its identifying halfedge.
 */
class LineTes_halfedge : public CGAL::Arr_halfedge_base<Curve> {
public:
  LineTes_halfedge();
  /** \brief Return previous halfedge along the segment
   */
  inline LineTes::Halfedge_handle get_prev_hf() {return prev_hf;}
  /** \brief Return next halfedge along the segment
   */
  inline LineTes::Halfedge_handle get_next_hf() {return next_hf;}
  /** \brief Test halfedge direction against the segment direction
   * \return true if the halfedge and the segment have the same direction,
   * false if they are opposite.
   */
  inline bool                     get_dir() {return dir;}
  /** \brief Return the halfedge length
   */
  inline double                   get_length() {return length;}
  /** \brief Return a handle of the segment containing the halfedge
   */
  inline LineTes::Seg_handle      get_segment() {return segment_ptr;}
  /** \brief Return a handle of the segment containing the halfedge
   * \note Redondant method. Is it really used?
   */
  inline LineTes::Seg_handle      segment(){return segment_ptr;}
  void                            set_length();
  /** \brief Set the length attribute of an halfedge
   */
  inline void                            set_length(double l) {length=l;}
  void                            set_prev_hf(LineTes::Halfedge_handle);
  void                            set_next_hf(LineTes::Halfedge_handle);
  void                            set_dir();
  void                            set_dir(bool d);
  /** \brief Define the segment containing the halfedge
   */
  inline void                     set_segment(LineTes::Seg_handle s){
    segment_ptr = s;}
  /** \brief Test whether a halfedge starting a segment is also ending the
   * segment
   * \pre the halfedge handle should not be null.
   */
  inline bool              is_one_edge_segment(){
    return next_hf==NULL_HALFEDGE_HANDLE; 
  }
 private:
  double                          length;
  LineTes::Halfedge_handle        next_hf;
  LineTes::Halfedge_handle        prev_hf;
  bool dir;
  LineTes::Seg_handle             segment_ptr;
};
/******************************************************************************/
/*                DEFINITIONS RELATED TO I/L-VERTEX PROCESSING                */
/******************************************************************************/
/** \brief Data about removal of an I-vertex
 *
 * An I-vertex in a line tessellation may be removed either by removing its
 * incident edge or by lengthening it until it meets another tessellation
 * edge. 
 */
struct IVertex {
  LineTes::Vertex_handle v; ///< handle of the I-vertex
  LineTes::Halfedge_handle e; ///< incident halfedge
  LineTes::Halfedge_handle esplit; ///< halfedge split by lengthening
  Point2 p; ///< location of new vertex added by lengthening
  double changed_length; ///< smallest length variation
  bool shorten; ///< true if smallest length variation due to shortening
};
/** \brief Comparison of I-vertices based on length variation 
 * \param v1 : data related to the first I-vertex.
 * \param v2 : data related to the second I-vertex.
 * \return true if first I-vertex removal induces a length change smaller than
 * the second I-vertex removal.
 */
inline bool compare_ivertices(IVertex v1,IVertex v2) {
  return v1.changed_length<v2.changed_length;
}


/******************************************************************************/
/*                   DEFINITION OF STRUCTURE ModList                          */
/******************************************************************************/

/** \brief Representation of a modification to be applied to a line tessellation
 *
 * The modification is represented by lists of deleted and new vertices, edges,
 * faces and segments.
 */
struct ModList {
  std::vector<Point2>               del_vertices; ///< list of deleted vertices
  std::vector<Point2>               add_vertices; ///< list of new vertices
  std::vector<Segment>              del_edges;    ///< list of deleted edges
  std::vector<Segment>              add_edges;    ///< list of new edges
  HPolygons                         del_faces;    ///< list of deleted faces
  HPolygons                         add_faces;    ///< list of new faces
  std::vector<std::vector<Point2> > del_segs;     ///< list of deleted segments
  std::vector<std::vector<Point2> > add_segs;     ///< list of new segments
};

/******************************************************************************/
/*                   DEFINITION OF CLASS PolygonImporter                      */
/******************************************************************************/

/** \brief Generation of a LineTes object from polygons with errors
 *
 * This class gather data and methods useful for trying to build a
 * line tessellation from a list of polygons which are supposed to be
 * tessellation faces. Polygon vertex coordinates are subject to
 * numerical errors. In particular, input polygons may have non simple
 * (even relatively) borders, their interiors may overlap, they may
 * not be space filling.
 *
 * Another problem handled here is the determination of the
 * tessellated domain when it is not defined explicitely and
 * separately. Methods for inferring the domain from input polygons
 * are provided.
 */
class PolygonImporter {
public:
  /** \typedef Side
   * \brief Side of a holed polygon
   *
   * The side is oriented. It is represented as a line segment together
   * with the polygon index (starting from 1. */
  typedef std::pair<Segment,Size> Side;
  /** \typedef Sides
   * \brief Vector of Side */
  typedef std::vector<Side> Sides;
  /** \typedef HistArrangement
   * \brief Arrangement of segments with history */
  typedef CGAL::Arrangement_with_history_2<Traits> HistArrangement;
  /** \typedef Arr_curve
   * \brief Segment of the arrangement with history */
  typedef Traits::Curve_2 Arr_curve;
  /** \typedef Curve_handle
   * \brief Handle to a segment of the arrangement with history */
  typedef HistArrangement::Curve_handle Curve_handle;
  /** \typedef SideCluster
   * \brief Group of polygons sides supposed to be part of a tessellation 
   * segment
   *
   * Each polygon side is represented by its projection onto the segment 
   * and the index of the polygon it belongs to (side_polygon). The side 
   * projection is represented by the curvilinear abscissae of its both 
   * ends, side_start and side_end. A curvilinear abscissa of 0 
   * (respectively 1) corresponds to the segment start (respectively end).
   */
  struct SideCluster {
    Segment representant; ///< tessellation segment
    std::vector<double> side_start; ///< side start abscissa on the segment 
    std::vector<double> side_end; ///< side end abscissa on the segment
    std::vector<double> side_index; ///< side index 
    std::vector<unsigned int> side_polygon; ///< polygon index of side
    Curve_handle curve_ref; ///< handle to the segment inserted into the arrangement
  };
  /** \typedef SideClusters
   * \brief Set of side groups */
  typedef std::vector<SideCluster> SideClusters;
  /** \typedef PolygonVote
   * \brief A polygon vote
   *
   * A pair where the first element is a (polygon) index and the second
   * element is a vote. */
  typedef std::pair<unsigned int,double> PolygonVote;
  /** \typedef PolygonVotes
   * \brief A set of polygon votes
   *
   *  Hash map where the polygon indices are used as hash keys. */
  typedef boost::unordered_map<unsigned int,double> PolygonVotes;
  // data members
  /** \brief Input polygons
   *
   * Input polygons represent the faces of the line tessellation to be
   * reconstructed. They may be defined approximately. */
  HPolygons input_polygons;
  /** \brief Results of the polygon side clustering
   *
   * Polygon sides supposed to be parts of the same tessellation segment are
   * grouped together by a clustering algorithm. For the time being, the 
   * clustering algorithm is not implemented and an external tool such as R
   * must be used. The clustering results may be read using the 
   * read_side_clusters method.
   */
  SideClusters side_clusters;
  /** \brief The reconstructed tessellation as an arrangement
   *
   * Arrangement with history is used in order to be able to associate 
   * halfedges to segments.
   */
  HistArrangement arr;
  // methods
  PolygonImporter();
  /** \brief Set input polygons
   * 
   * Input polygons are holed polgygons supposed to represent
   * approximately the tessellation faces. */
  inline void set_polygons(HPolygons hpolys) {input_polygons = hpolys;};
  /** \brief Get input polygons */
  inline HPolygons get_polygons() {return input_polygons;};
  void read_polygons(std::istream&);
  void read_side_clusters(std::istream&);
  Sides get_polygon_sides();
  /** \brief Set side clustering results*/
  inline void set_side_clusters(SideClusters sc) {side_clusters = sc;}
  /** \brief Get stored side clustering results*/
  inline SideClusters get_side_clusters() {return side_clusters;};
  void insert_segments(double expand=0.0);
  Size number_of_I_vertices(); 
  Size remove_I_vertices(double within=infinity_double);
  SideClusters::iterator fetch_side_cluster(Curve_handle);
  HPolygons get_faces();
  PolygonVotes polygon_sides(HistArrangement::Halfedge_handle);
  PolygonVotes polygon_sides(HistArrangement::Face_handle);
  PolygonVote elected_polygon(HistArrangement::Face_handle);
  NT compare(HistArrangement::Face_handle,unsigned int,NT);    
  NT goodness_of_fit(NT eps,double weight_missing);
};

/******************************************************************************/
/*                   DEFINITION OF CLASS TTessel                              */
/******************************************************************************/

/** \brief Dynamic T-tessellation
 *
 * Possible updates are splits, merges and flips.
 */
class TTessel : public LineTes {

public:
  /** \brief Abstract class representing a modification of a T-tessellation
   *
   * Concrete derived classes must implement the computation of
   * elements of the T-tessellation that are modified by the
   * modification.
   */ 
  class Modification {
  public:
    virtual ModList modified_elements() ;
  };
  /** \brief A split that can be applied to a T-tessellation
   *
   * A split is defined as the insertion of a line segment in
   * a tessellation face. The inserted line segment must not
   * intersect the existing face boundary except at its ends.
   * When the inserted line segment connects two points on the
   * face outer boundary, the face is divided into two new smaller
   * faces. 
   * \image html split_square_domain.svg "A standard split"
   *
   * A split is represented by both ends \f$(p_1,p_2)\f$ of the inserted line
   * segment and handles \f$(e_1,e_2)\f$ to the halfedges where they lie. 
   * These data can be accessed 
   * through Split::get_e1, Split::get_p1, Split::get_e2 and 
   * Split::get_p2 methods.
   *
   * The Split(Halfedge_handle,double,double) constructor takes as input
   * arguments the handle of the halfedge where the splitting line segment
   * starts, the relative position (between 0 and 1) of the
   * starting point along that halfedge and the angle betweeen the halfedge
   * and the splitting line segment.
   *
   * The way a T-tessellation is modified by a split can be predicted using
   * Split::modified_elements method. A split has several effects on 
   * tessellation elements:
   * - Two edges bounding the split face are split.
   * - Two vertices along these two edges are created.
   * - The face where the line segment is inserted is
   *   either modified or 
   *   split into two faces.
   * - A new non-blocking segment is created.
   *
   * The insertion of the line segment may not result in a face splitting 
   * only when the face has at least one hole. The examples below show 
   * cases that can occur with a holed face. The standard case where the 
   * inserted line segment joins two points of the outer boundary is not
   * shown.
   * \image html split_holed_domain.svg "Splits of a holed face"
   * Above, the boundaries of the holed face are shown in black and the
   * inserted line segment in red. Left: the face is modified by a merging
   * of an inner boundary to the outer boundary. Middle: the face is modified
   * by a merging of two inner boundaries. Right: the face is split in two
   * faces, one being included into the other one.
   * \sa TTessel::update(Split), TTessel::Modification::modified_elements
   */
  class Split : public Modification {
  public:
    /** \name Initialization */
    /** \{ */
    Split();
    Split(Halfedge_handle, double=.5,double=CGAL_PI/2);
    /** \} */

    /** \name Access */
    /** \{ */
    /** \brief Return a handle to one of the split edges
     *
     * The other split edge can be accessed using Split::get_e2.*/
    inline Halfedge_handle get_e1(){return e1;}
    /** \brief Get the location of the new vertex on one of the new edge
     *
     * The new vertex lies on the new edge returned by Split::get_e1*/
    inline Point2          get_p1(){return p1;}
    /** \brief Return a handle to one of the split edges
     *
     * The other split edge can be accessed using Split::get_e1.*/
    inline Halfedge_handle get_e2(){return e2;}
    /** \brief Get the location of the new vertex on one of the new edge
     *
     * The new vertex lies on the new edge returned by Split::get_e2*/
    inline Point2          get_p2(){return p2;}
    /** \} */

    /** \name Computations */
    /** \{ */
    virtual                ModList modified_elements();
    /** \} */

    /** \name Checking */
    /** \{ */
    bool                   is_valid();
    /** \} */

  private:
    Halfedge_handle         e1; // Split edge
    Halfedge_handle         e2; // Split edge
    Point2                  p1; // New point on e1
    Point2                  p2; // New point on e2
  };
  /** \brief A merge that can be applied to a T-tessellation
   *
   * A merge is the removal of a non-blocking segment. When the
   * segment separates two distinct faces, they are merged.
   * \image html merge_square_domain.svg "A standard merge"
   * In the example above, the merge shown is the only possible
   * one since there is only one non-blocking segment.
   *
   * A merge is represented by the (half)edge to be removed.
   *
   * A merge has several effects on tessellation elements:
   * - Two vertices are removed.
   * - Two pairs of edges are merged.
   * - Either two faces are replaced by a single one or a face is modified.
   * - A non-blocking segment is removed.
   * - The status (blocking or not) of segments incident to 
   *   the removed segment may change.
   *
   * The removal of a non-blocking segment may not result into a merging of
   * two faces. This happens in particular when the segment is associated
   * with two halfedges that belong to the same (non simple) boundary. 
   * Examples can be found in the documentation of class Split.
   */
  class Merge : public Modification { 
  public:
    /** \name Initialization */
    /** \{ */
    Merge();
    Merge(Halfedge_handle);
    /** \} */

    /** \name Access */
    /** \{ */
    /** \brief Return a handle to the suppressed edge*/
    inline Halfedge_handle get_e(){return e;}
    /** \} */

    /** \name Computations */
    /** \{ */
    virtual                ModList modified_elements();
    /** \} */

    /** \name Checking */
    /** \{ */
    bool                   is_valid();
    /** \} */
  private:
    Halfedge_handle e;
    inline Halfedge_handle e1(){return e->twin()->next();}
    inline Halfedge_handle e2(){return e->next();}
    inline Point2           p1(){return e->source()->point();}
    inline Point2           p2(){return e->target()->point();}

  };
  /** \brief A flip that can be applied to a T-tessellation
   *
   * A flip is a shortening of a segment combined with a lengthening of
   * another segment. A standard flip is shown below.
   * \image html flip_square_domain.svg "A standard flip"
   * The shortened segment must be blocking (made of more than one edge)
   * and the removed halfedge must be at one of its ends. 
   *
   * A flip is represented by its removed halfedge. That halfedge is
   * to be provided to Flip(Halfedge_handle) constructor.
   *
   * Method Flip::modified_elements predicts how a flip modifies a 
   * T-tessellation.
   */
  class Flip : public Modification {
  public:
    /** \name Initialization */
    /** \{ */
    Flip();
    Flip(Halfedge_handle);
    /** \} */

    /** \name Access */
    /** \{ */
    /** \brief Return a handle to the suppressed edge*/
    inline Halfedge_handle get_e1(){return e1;}
    /** \brief Return a handle to the split halfedge
     *
     * The split halfedge is the halfedge where the extending segment
     * ends. */
    inline Halfedge_handle get_e2(){return e2;}
    /** \brief Return the location of the new vertex
     *
     * The new vertex lies at the end of the extending segment.*/
    inline Point2          get_p2(){return p2;}
    /** \brief Return true when the halfedge flips towards the righthand-side*/
    inline bool to_right(){return right;}
    /** \} */

    /** \name Computations */
    /** \{ */
    virtual                ModList modified_elements();
    /** \} */
  private:
    Halfedge_handle        e1; // Suppressed edge
    Halfedge_handle        e2; // Split edge
    Point2                 p2; // New vertex
    bool                   right; // True if flip to right
  };

  /** \typedef Split_list
   * \brief A list of splits*/
  typedef std::vector<Split>   Split_list;
  /** \typedef Merge_list
   * \brief A list of merges*/
  typedef std::vector<Merge>   Merge_list;
  /** \typedef Flip_list
   * \brief A list of flips*/
  typedef std::vector<Flip>    Flip_list;
  /** \typedef Split_list_iterator
   * \brief Iterator for accessing splits in a list*/
  typedef Split_list::iterator Split_list_iterator;
  /** \typedef Merge_list_iterator
   * \brief Iterator for accessing merges in a list*/
  typedef Merge_list::iterator Merge_list_iterator;
  /** \typedef Flip_list_iterator
   * \brief Iterator for accessing flips in a list*/
  typedef Flip_list::iterator  Flip_list_iterator;

  /** \name Initialization */
  /** \{ */
  TTessel();
  TTessel(LineTes&);
  virtual void         insert_window(Rectangle);
  virtual void         insert_window(Polygon&);
  virtual void         insert_window(HPolygons&);
  /** \} */
  /** \name Modification */
  /** \{ */
  Halfedge_handle      update(Split);
  Face_handle          update(Merge);
  Halfedge_handle      update(Flip);
  void                 clear(bool remove_window=true);
  Split                propose_split();
  Merge                propose_merge();
  Flip                 propose_flip();
  Split_list           split_sample(Size);
  Split_list           poisson_splits(double);
  Merge_list           all_merges();
  Flip_list            all_flips();
  /** \} */
  /** \name Access */
  /** \{ */
  Seg_sublist_iterator blocking_segments_begin();
  Seg_sublist_iterator blocking_segments_end();
  HPolygons            all_faces();
  Seg_sublist_iterator non_blocking_segments_begin();
  Seg_sublist_iterator non_blocking_segments_end();
  /** \} */
  /** \name Features */
  /** \{ */
  Size                 number_of_blocking_segments();
  Size                 number_of_non_blocking_segments();
  /** \brief Sum of all halfedge lengths
   *
   * The sum of halfedge lengths is equal to the sum of cell perimeters.
   */
  inline double        get_total_internal_length(){return int_length;}
  /** \} */
  /** \name Checking and testing*/
  /** \{ */
  bool                 is_valid(bool verbose=false);
  /** \} */
  /** \name Input/output */
  /** \{ */
  void                 print_blocking_segments();
  void                 print_non_blocking_segments();
  void                 printRCALI(std::ostream&);
  /** \} */
private:
  Seg_sublist          blocking_segments;
  Seg_sublist          non_blocking_segments;
  double               int_length;
  Halfedge_handle      random_halfedge();
  inline double        total_internal_length(){
                       /* Return the total length of halfedges bounding an internal face. */
                        return int_length;
                       }
  Halfedge_handle      length_weighted_random_halfedge();
  Halfedge_handle      alt_length_weighted_random_halfedge();
  std::vector<Halfedge_handle> alt_length_weighted_random_halfedge(unsigned int);
};

/******************************************************************************/
/*                   DEFINITION OF STRUCTURE CatItems                         */
/******************************************************************************/

/** \struct CatItems
 * \brief Storage structure for categorized items
 *
 * Item categories: vertices, edges, faces and segs.
 */
template <typename T>
struct CatItems {
  std::vector<T> vertices; ///< a vector of T objects
  std::vector<T> edges;    ///< a vector of T objects
  std::vector<T> faces;    ///< a vector of T objects
  std::vector<T> segs;     ///< a vector of T objects
  const CatItems<T>& operator+=(const CatItems<T> &);
  const CatItems<T>& operator-=(const CatItems<T> &);
  bool IsEmpty();
};
template <typename T>
const CatItems<T> operator+(const CatItems<T> &, const CatItems<T> &);
template <typename T>
const CatItems<T> operator*(const double &, const CatItems<T> &);
template <typename T>
const CatItems<T> operator-(const CatItems<T> &, const CatItems<T> &);
template <typename T>
const CatItems<T> operator*(const CatItems<T> &, const CatItems<T> &);

/** \brief Increment operator for CatItems objects
 */
template<typename T>
const CatItems<T>& CatItems<T>::operator+=(const CatItems<T> &incr) {
  *this = (*this)+incr;
  return *this;
}

/** \brief Decrement operator for CatItems objects
 */
template<typename T>
const CatItems<T>& CatItems<T>::operator-=(const CatItems<T> &decr) {
  *this += -1.0*decr;
  return *this;
}

/** \brief Addition operator for CatItems objects
 */
template<typename T>
const CatItems<T> operator+(const CatItems<T> &x1, 
				const CatItems<T> &x2) {
  CatItems<T> res = CatItems<T>(x1);
  for(unsigned int i=0; i!= x2.vertices.size(); i++) {
    res.vertices[i] += x2.vertices[i];
  }
  for(unsigned int i=0; i!= x2.edges.size(); i++) {
    res.edges[i] += x2.edges[i];
  }
  for(unsigned int i=0; i!= x2.faces.size(); i++) {
    res.faces[i] += x2.faces[i];
  }
  for(unsigned int i=0; i!= x2.segs.size(); i++) {
    res.segs[i] += x2.segs[i];
  }
  return res;
}

/** \brief Product of a CatItems object with a scalar
 */
template<typename T>
const CatItems<T> operator*(const double &scalar, const CatItems<T> &x) {
  CatItems<T> res = CatItems<T>(x);
  for(unsigned int i=0; i!= x.vertices.size(); i++) {
    res.vertices[i] = scalar*res.vertices[i];
  }
  for(unsigned int i=0; i!= x.edges.size(); i++) {
    res.edges[i] = scalar*res.edges[i];
  }
  for(unsigned int i=0; i!= x.faces.size(); i++) {
    res.faces[i] = scalar*res.faces[i];
  }
  for(unsigned int i=0; i!= x.segs.size(); i++) {
    res.segs[i] = scalar*res.segs[i];
  }
  return res;
}

/** \brief Subtraction of two CatItems objects
 */
template<typename T>
const CatItems<T> operator-(const CatItems<T>& x1, const CatItems<T>& x2) {
  CatItems<T> res = x1+(-1.0*x2);
  return res;
}

/** \brief Term by term product of CatItems objects
 */
template<typename T>
const CatItems<T> operator*(const CatItems<T> &x1, 
			      const CatItems<T> &x2) {
  CatItems<T> res = CatItems<T>(x1);
  for(unsigned int i=0; i!= x2.vertices.size(); i++) {
    res.vertices[i] *= x2.vertices[i];
  }
  for(unsigned int i=0; i!= x2.edges.size(); i++) {
    res.edges[i] *= x2.edges[i];
  }
  for(unsigned int i=0; i!= x2.faces.size(); i++) {
    res.faces[i] *= x2.faces[i];
  }
  for(unsigned int i=0; i!= x2.segs.size(); i++) {
    res.segs[i] *= x2.segs[i];
  }
  return res;
}

/** \typedef CatVector
 * \brief Storage structure for vector of doubles associated with categories
 * vertices, edges, faces or segs.*/
typedef CatItems<double> CatVector;
void fill(CatVector&,FVector&);
// Don't know how to specify standard user-defined conversion for CatVector :-(
std::vector<double> asVectorOfDoubles(const CatVector&);
FVector asFVector(const CatVector&);
std::ostream& operator<<(std::ostream &, const CatVector &);
double sum(const CatVector &);
/** \typedef CatMatrix
 * \brief Storage structure for matrices (vector of vectors) of doubles
 * associated with categories vertices, edges, faces or segs.*/
typedef CatItems<CatVector> CatMatrix;
FMatrix asFMatrix(const CatMatrix&);
CatMatrix outer(const CatVector &, const CatVector &);

/******************************************************************************/
/*                   DEFINITION OF STRUCTURE Features                         */
/******************************************************************************/
/** \brief Set of vertex, edge, face and segment features*/
struct Features {
  std::vector<double (*)(Point2,TTessel*)>              vertices;///< vertex features (vector of functions)
  std::vector<double (*)(Segment,TTessel*)>             edges;///< edge features
  std::vector<double (*)(HPolygon,TTessel*)>            faces;///< face featuress
  std::vector<double (*)(std::vector<Point2>,TTessel*)> segs;///< segment features
};



/******************************************************************************/
/*                   DEFINITION OF CLASS Energy                               */
/******************************************************************************/

/** \brief Energy of a Gibbsian random T-tessellation
 *
 * The energy must have the following form:
 * \f[
 * \sum_i \theta_i \phi_i(T)
 * \f]
 * where the \f$\theta_i\f$'s are real-valued parameters and the 
 * \f$\phi_i\f$'s are statistics of the type
 * \f[
 * \phi_i(T) = \sum_x f_i(x)
 * \f]
 * where the sum runs over all vertices, or all edges, or all faces, or all 
 * segments of the T-tessellation.
 *
 * Methods for defining the form of the energy and the values of the
 * \f$\theta_i\f$'s are provided. In order to compute the energy of a
 * given T-tessellation, the latter must be attached to the %Energy
 * object through the Energy::set_ttessel method. If the T-tessellation is
 * modified using splits, merges or flips, the energy value may be
 * updated using the Energy::variation and Energy::add_value methods.
 */
class Energy{

public:
  /** \name Initialization */
  /** \{ */
  Energy();
  void                     add_features_vertices(double (*)(Point2,TTessel*));
  void                     add_features_edges(double (*)(Segment,TTessel*));
  void                     add_features_faces(double (*)(HPolygon,TTessel*));
  void                     add_features_segs(double (*)(std::vector<Point2>,TTessel*) );
  void                     add_theta_vertices(double);
  void                     add_theta_edges(double);
  void                     add_theta_faces(double);
  void                     add_theta_segs(double);
  /** \brief Set the theta parameter */
  inline void              set_theta(CatVector par){theta = par;}
  void                     set_ttessel(TTessel*);
  /** \} */

  /** \name Access */
  /** \{ */
  /** \brief Return the pointer to the associated TTessel object */
  inline TTessel*          get_ttessel(){return ttes;}
  /** \brief Return the current value of the energy for the associated
      T-tessellation */
  inline double            get_value(){return value;}
  /** \brief Return the theta parameter */
  inline CatVector         get_theta(){return theta;}
  /** \} */

  /** \name Modification */
  /** \{ */
  /** \brief Increment the current value of the energy */
  inline void              add_value(double delta){value +=delta;}
  void                     del_theta_vertices();
  void                     del_theta_edges();
  void                     del_theta_faces();
  void                     del_theta_segs();
  void                     del_features_vertices();
  void                     del_features_edges();
  void                     del_features_faces();
  void                     del_features_segs();
  /** \} */

  /** \name Computation */
  /** \{ */
  double                   variation(TTessel::Modification&);
  CatVector                statistic_variation(TTessel::Modification&);
  /** \} */
private:
  TTessel*                 ttes;
  CatVector                theta;
  Features                 features;
  double                   value;
};


/******************************************************************************/
/*                   DEFINITION OF STRUCTURE ModifCounts                      */
/******************************************************************************/
/** \brief Data summarizing changes occuring through SMF iterations*/
struct ModifCounts {
  int  proposed_S; ///< number of proposed splits
  int  accepted_S; ///< number of accepted splits
  int  proposed_M; ///< number of proposed merges
  int  accepted_M; ///< number of accepted merges
  int  proposed_F; ///< number of proposed flips
  int  accepted_F; ///< number of accepted flips
};



/******************************************************************************/
/*                   DEFINITION OF CLASS SMFChain                             */
/******************************************************************************/

/** \brief Simulation of a Gibbsian T-tessellation
 *
 * Updates are splits, merges and flips. Transition rules according a
 * Metropolis-Hastings-Green construction.
 *
 * The probability distribution to be simulated must be instanciated
 * as an Energy object. It is associated to the %SMFChain object
 * through method SMFChain::set_energy. The probabilities of the three
 * types of modifications are set up through method
 * SMFChain::set_smf_prob. The Markov chain advances by invoking its
 * SMFChain::step method.
 */
class SMFChain {
public:
  /** \enum ModType
   * \brief Types of modifications applicable to a T-tessellation*/
  enum ModType {
    SPLIT = 1,///< split (of a face)
    MERGE ,///< merge (of two faces separated by a non-blocking segment)
    FLIP ///< a segment is shortened, another one is lengthened
  };

  /** \name Initialization */
  /** \{ */
  /** \brief Default constructor
   *
   * Does nothing special. */
  SMFChain() {};
  SMFChain(Energy *,double, double);
  /** \brief Set the energy function of the model to be simulated
   */
  inline void    set_energy(Energy *e) {engy = e;};
  void           set_smf_prob(double,double);
  /** \} */

  /** \name Access */
  /** \{ */
  /** \brief Get the energy function of the model to be simulated
   */
  inline Energy* get_energy() {return engy;};
  /** \} */

  /** \name Computation */
  /** \{ */
  ModType        propose_modif_type();
  double         Hasting_ratio(TTessel::Split&,double*);
  double         Hasting_ratio(TTessel::Merge&,double*);
  double         Hasting_ratio(TTessel::Flip&,double*);
  /** \} */

  /** \name Modification */
  /** \{ */
  ModifCounts    step(unsigned long int=1);
  /** \} */

private:
  Energy  *engy;
  double  p_split;
  double  p_merge;
  double  p_flip;
  double  p_split_merge;
};

/******************************************************************************/
/*                   DEFINITION OF CLASS PseudoLikDiscrete                    */
/******************************************************************************/
/** \brief Discrete approximation of the log-pseudolikelihood of a Gibbsian T-tessellation
 *
 * PseudoLikDiscrete is a technical class that can be used when implementing
 * inference algorithms. An example of user-end class based on 
 * %PseudoLikDiscrete is PLInferenceNOIS.
 *
 * A pseudolikelihood method for Gibbsian T-tessellation was
 * introduced in Research Report "Pseudolikelihood inference for
 * Gibbsian T-tessellations... and point processes" available on
 * [HAL](http://hal.archives-ouvertes.fr/hal-01234119) and
 * [arXiv](http: //arxiv.org/abs/1512.08407). Consider a parametric
 * family of energy functions \f$E_\theta\f$. For each \f$\theta\f$,
 * each T-tessellation \f$T\f$ and each local modification (a split,
 * a merge or a flip) \f$o\f$ of
 * \f$T\f$ , let
 * \f[
 * \Delta E_\theta(o,T) = E_\theta(oT)-E_\theta(T).
 * \f]
 * Note that energy variations may be predicted based on the 
 * Energy::variation method.
 * The discrete approximation 
 * of log-pseudolikelihood
 * is the sum of two terms. The first term is associated to merges and
 * splits:
 * \f[
 * -\sum_{m\in\mathbf{M}_T} \Delta E_\theta(m,T) -
 * \frac{u(T)}{\pi|S|}
 * \sum_{s\in S} \exp(-\Delta E_\theta(s,T)) ds
 * \f]
 * where \f$\mathbf{M}_T\f$ is the set of merges that can be applied to 
 * the tessellation \f$T\f$ and \f$S\f$ is a sample of splits
 * that can be applied to \f$T\f$, \f$u(T)\f$ is the perimeter sum 
 * of faces of \f$T\f$ and \f$|S|\f$ is the size of \f$T\f$. Splits in 
 * \f$S\f$ are called dummy splits.   
 *
 * The second term is associated to flips:
 * \f[
 * -\sum_{f\in\mathbf{F}_T} \Delta E_\theta(f,T) -
 * \sum_{\mathbf{F}_T} \exp(-\Delta E_\theta(f,T))
 * \f]
 * where \f$\mathbf{F}_T\f$ is the set of flips that can be applied to 
 * the tessellation \f$T\f$.
 *
 * A PseudoLikDiscrete object must be made aware of the parametric family
 * of energy functions. This can be done using one of the constructors or
 * PseudoLikDiscrete::SetEnergy method. The value of \f$\theta\f$ assigned
 * in the Energy object does not matter. The observed T-tessellation 
 * \f$T\f$ is the tessellation associated with the Energy object.
 *
 * The second ingredient of the discrete approximation of 
 * log-pseudolikelihood is the sample \f$S\f$ of dummy splits. At start, it
 * is empty. Splits drawn independently
 * according to the uniform distribution can be added to the current
 * sample of dummy splits.
 *
 * The value of the approximation can be obtained by 
 * PseudoLikDiscrete::GetValue method. Other methods useful for
 * optimization are PseudoLikDiscrete::GetGradient and 
 * PseudoLikDiscrete::GetHessian.
 *
 * If the energy function is written as
 * \f[
 * E_\theta(T) = \sum_i\theta_i\phi_i(T),
 * \f]
 * then
 * \f[
 * \sum_{m\in\mathbf{M}_T}\Delta E_\theta(m,T) =
 * \sum_i \theta_i \sum_{m\in\mathbf{M}_T} (\phi_i(mT)-\phi_(T)).
 * \f]
 * Therefore when assessing the value above for various \f$\theta\f$'s, it
 * saves computation effort to precompute the vector of
 * \f[
 * \sum_{m\in\mathbf{M}_T} (\phi_i(mT)-\phi_(T)).
 * \f]
 * This computation is performed when invoking the constructor with an Energy
 * object as input argument. And the result is stored as a private member
 * of the %PseudoLikDiscrete object. Concerning the term depending on splits
 * \f[
 * \sum_{s\in S}\exp(-\Delta E_\theta(s,T)),
 * \f]
 * one needs to precompute and store a larger vector containing the
 * \f[
 * \phi_i(sT)-\phi_i(T)
 * \f]
 * for all \f$s\in S\f$ and all \f$i\f$'s. Moreover this vector is not computed
 * a build time but each time the vector of dummy splits is incremented.
 * Similarly, the vector of 
 * \f[
 * \phi_i(fT)-\phi_i(T)
 * \f]
 * for all flips \f$f\in\mathbf{F}_T\f$ and all \f$i\f$'s and its sum over 
 * the flips are computed and stored at build time.
 */
class PseudoLikDiscrete {
 public:
  /** \name Initialization */
  /** \{ */
  /** \brief Default constructor */
  PseudoLikDiscrete() {};
  PseudoLikDiscrete(Energy*);
  void SetEnergy(Energy*);
  /** \} */

  /** \name Modification */
  /** \{ */
  void AddSplits(Size);
  void AddGivenSplit(TTessel::Split); // only for debugging?
  void ClearSplits();
  /** \} */

  /** \name Access */
  /** \{ */
  /** \brief Get the Energy data member
   */
  inline Energy* GetEnergy() {return engy;};
  /** \} */

  /** \name Computations */
  /** \{ */
  double GetValue(CatVector);
  CatVector GetGradient(CatVector);
  CatMatrix GetHessian(CatVector);
  /** \} */

 private:
  // data members
  Energy *engy;
  CatVector sum_stat_merges;
  CatVector sum_stat_flips;
  std::vector<CatVector> stat_splits;
  std::vector<CatVector> stat_flips;
  // method members
};
  
/******************************************************************************/
/*                   DEFINITION OF CLASS PLInferenceNOIS                      */
/******************************************************************************/

/** \brief Pseudolikelihood inference of a Gibsian T-tessellation with alternate updates
 *
 * Numerical maximization of the log-pseudolikelihood. The acronym NOIS
 * stands for Newton Optimization and Increasing Splitting. Each step consists
 * of a 
 * double update. First, dummy splits are added to the discrete approximation
 * of the pseudolikelihood (see class PseudoLikDiscrete) in order to improve 
 * its accuracy. The number of added 
 * dummy splits is equal to the number of merges that can be applied to the
 * observed T-tessellation. Second, the current
 * parameter estimate is updated according to Newton's method:
 * \f[
 * \theta \leftarrow \theta+\lambda H^{-1}G,
 * \f]
 * where \f$G\f$ and \f$H\f$ are the gradient and the Hessian of the
 * approximated log-pseudolikelihood. The parameter \f$\lambda\f$
 * referred to as the stepsize parameter is a tuning parameter of 
 * Newton's method. Iterations of the maximization algorithm
 * are performed using the PLInferenceNOIS::Step method.
 *
 * The stopping criterion implemented in %PLInferenceNOIS is as follows: stop
 * if the number of steps exceeds a predefined maximal number of if the
 * approximated log-pseudolikelihood (say lp) has been reduced by a factor less
 * than tol*(abs(lp)+tol) during the last step. The algorithm can be run until 
 * the stopping criterion is met using PLInferenceNOIS::Run method.
 */
class PLInferenceNOIS : public PseudoLikDiscrete{
 public:
  /** \name Initialization */
  /** \{ */
  PLInferenceNOIS(Energy*,bool=false,double=1.0);
  /** \brief Set step size used in Newton's method */
  inline void SetStepSize(double lambda) { stepsize = lambda;};
  /** \} */

  /** \name Access */
  /** \{ */
  /** \brief Return step size used in Newton's method */
  inline double GetStepSize() {return stepsize;};
  CatVector GetEstimate();
  /** \brief Return all intermediate parameter estimates
   */
  inline std::vector<CatVector> GetEstimates() {return estimates;};
  /** \brief Return all intermediate log-pseudolikelihood values
   */
  inline std::vector<double>    GetValues() {return values;};
  /** \} */

  /** \name Computations */
  /** \{ */
  void Step(unsigned int=1);
  unsigned int Run(double=0.05, unsigned int=100);
  /** \} */

 private:
  double                 stepsize; /** Stepsize to be used for Newton's method.*/
  CatVector              previous_theta; /** Estimate at previous step. Data
					  *  used for computing the step size.*/
  double                 previous_lpl; /** Log-pseudolikelihood at previous step.
                                        *  Data used for the stopping criterion. 
					*/
  bool                   store_path; /** If true, intermediate estimates and
				      * values of the log-pseudolikelihood are
				      * stored.*/
  std::vector<CatVector> estimates; /** For storing intermediate estimates.*/
  std::vector<double>    values; /** For storing intermediate values of the
				  * log-pseudolikelihood.*/
};

/******************************************************************************/
/*                       FUNCTION DECLARATIONS                                */
/******************************************************************************/

bool                     are_aligned(Point2 p, Point2 q, Point2 r, 
				     bool verbose=false);
bool                     is_inside(Point2, HPolygon&);
bool                     is_inside(Point2, HPolygons&);
Segment                  clip_segment_by_convex_polygon(Segment, Polygon);
Segment clip_segment_by_polygon(Segment, Polygon) 
     throw(std::domain_error const&);
Segment clip_segment_by_polygon(Segment, HPolygon) 
     throw(std::domain_error const&);
Segment clip_segment_by_polygon(Segment, HPolygons) 
     throw(std::domain_error const&);
double                   precompute_lengthening(Arrangement::Halfedge_handle,
						Arrangement::Halfedge_handle*,
						Point2*);
bool                     exist_halfedge(LineTes&,LineTes::Halfedge_handle);
LineTes::Halfedge_handle compute_prev_hf(LineTes::Halfedge_handle);
LineTes::Halfedge_handle compute_next_hf(LineTes::Halfedge_handle);
void                     set_junction(LineTes::Halfedge_handle,LineTes::Halfedge_handle);
LineTes::Halfedge_handle find_halfedge(LineTes&,Point2,Point2);
bool                     is_a_T_vertex(LineTes::Vertex_handle, 
				       bool verbose=false);
unsigned long int        number_of_internal_vertices(TTessel&);
// conflict overloading/default argument -> aface2poly instead of face2poly
HPolygon              aface2poly(PolygonImporter::HistArrangement::Face_handle,
	   	                 bool poly_simplify=true);
HPolygon                 face2poly(TTessel::Face_handle, bool simplify=true);
double                   angle_between_vectors(Vector v1,Vector v2);
double                   sum_of_faces_squared_areas(TTessel*);
double                   sum_of_min_angles(TTessel*);
double                   sum_of_angles_obt(TTessel*);
double                   sum_of_segment_squared_sizes(TTessel* t);
Polygons                 boundaries(HPolygon);
Polygon                  simplify(Polygon) 
     throw(std::domain_error const&);;
HPolygon                 simplify(HPolygon);
Point2                   ray_exit_face(Rayon&,LineTes::Face_handle&,
				       LineTes::Halfedge_handle&);
Polygon                  polygon_insert_edge(PECirc,
					     PECirc,
					     Point2,Point2);
HPolygons hpolygon_insert_edge(Polygons &hpoly, Size, Size, 
			       PECirc, PECirc, Point2, Point2);
HPolygon hpolygon_remove_edge(Polygons&, Polygons&, Size, Size, 
			      PECirc, PECirc,bool);
Polygon                  polygon_remove_edge(PECirc&,
					     PECirc&);

bool                     has_holes(LineTes::Face_handle&);
bool                     has_holes(HPolygon&);
std::vector<bool>        filter_holes(HPolygon::Hole_const_iterator,
				      HPolygon::Hole_const_iterator,
				      Polygon&);
std::vector<bool>        filter_holes(Polygons::const_iterator,
				      Polygons::const_iterator,
				      Polygon&);
Size find_edge_in_polygons(Polygons&, Segment, PECirc&);
Size find_edge_in_hpolygon(HPolygon&, Segment, PECirc&);
NT curvilinear_coordinate(Point2,Segment);
double interval_free_length(double,double,
                            std::vector<std::pair<double,double> >);
NT area(HPolygon);
Points digitize(Points,NT);
NT squared_Hausdorff_distance(Points&,Points&);
/** \defgroup features Features of T-tessellations
 *
 * Functions that can be used by Energy objects for defining a Gibbs model
 * of T-tessellation. %Features are measured on vertices, edges, faces 
 * or segments. @{*/
double                   is_point_inside_window(Point2, LineTes*);
double                   is_point_inside_window(Point2, TTessel*);
double                   edge_length(Segment,TTessel*);
double                   seg_number(std::vector<Point2>,TTessel*);
double                   is_segment_internal(std::vector<Point2>,TTessel*);
/** \brief Return -1
 *
 * \param s : a tessellation segment as a vector of its 2 ends.
 * \param t : the tessellation to be considered.
 *
 * Silly function that can be used by an Energy object for specifying
 * minus the number of segment as a tessellation feature.*/
inline double            minus_is_segment_internal(std::vector<Point2> s,
						   TTessel* t){
  return -is_segment_internal(s,t);}
double                   face_number(HPolygon,TTessel*);
double                   face_area_2(HPolygon,TTessel*);
double                   face_perimeter(HPolygon, TTessel*);
double                   face_shape(HPolygon,TTessel*);
double                   face_sum_of_angles(HPolygon,TTessel*);
double                   min_angle(HPolygon, TTessel*);
double                   segment_size_2(std::vector<Point2>,TTessel*);
/** @}*/
