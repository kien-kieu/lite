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
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
/* Calcul exact pour +,-,*,/. On peut sauvegarder la tessellation dans
un fichier texte sans perte de précision. Par contre, la racine carrée
n'est implémentée.*/ 
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Linear_algebraCd.h>
#include <CGAL/double.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_default_dcel.h>
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
 * Use of non-standard halfedges (class LineTes_halfedge)*/
typedef CGAL::Arr_dcel_base<CGAL::Arr_vertex_base<Point2>,
		            LineTes_halfedge,
		            CGAL::Arr_face_base>                Dcel;
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

// definition that prevent type expansion in the documentation generated by Rcpp
class Polygons : public std::vector<Polygon> {};

/******************************************************************************/
/*                  DÉFINITION DE LA CLASSE LineTes                           */
/******************************************************************************/

/** \brief Polygonal tessellation of a bounded rectangular domain
 *
 * This class inherits from the Arrangement_2 class provided by
 * CGAL. Its main feature concerns the so-called segments. A segment
 * is defined as a maximal subset of edges that are aligned and
 * connected. The %LineTes class is of interest for tessellations with
 * segments that consist of more than one edge.
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
    Seg();
    Seg(Halfedge_handle);
    Halfedge_handle         halfedges_start();
    Halfedge_handle         halfedges_end();
    Size                    number_of_edges();
    bool                    number_of_edges_is_greater_than(unsigned int);
    Point2                  pointSource();
    Point2                  pointTarget();
    Points                  list_of_points();
    /** \brief Return the halfedge handle identifying the segment*/
    inline Halfedge_handle  get_halfedge_handle(){return e;}
    void                    set_halfedge_handle(Halfedge_handle);
    void                    print(bool endsOnly=false);
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
   
  LineTes();
  
  virtual void                  insert_window(Rectangle);
  virtual void                  insert_window(Polygon&);
  Halfedge_handle               split_from_edge(Halfedge_handle, Point2, 
					             Halfedge_handle, Point2);
  Halfedge_handle               split_from_vertex(Vertex_handle, Halfedge_handle,
					          Point2);
  Face_handle                   suppress_edge(Halfedge_handle);
  void                          clear();
  /** \brief Return an iterator pointing to the beginning of the list of
   * tessellation segments
   */
  inline Seg_list_iterator      segments_begin(){return all_segments.begin();}
  /** \brief Return an iterator pointing to the end of the list of
   * tessellation segments
   */
  inline Seg_list_iterator      segments_end(){return all_segments.end();}
  /** \brief Return the total number of segments
   *
   * Boundary and internal segments are included
   */
  inline Size                   number_of_segments(){return all_segments.size();}
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
  void                          print_all_segments(bool endsOnly=false);
  /** \brief Return the domain that is tessellated 
   */
  inline Polygon                get_window(){return window;} 
  Size                          number_of_window_edges();
  double                        get_window_perimeter();
  int                           is_on_boundary(Halfedge_handle);
  /** \brief Test whether a segment lies on the domain boundary
   *
   * \param s : a handle to the segment to be tested
   */
  inline bool                   is_on_boundary(Seg_handle s){return is_on_boundary(s->get_halfedge_handle())>0;}
  Seg_handle                    insert_segment(Segment);
  void                          remove_ivertices(Size imax=0,
						 bool verbose=false);
  void                          remove_lvertices(Size imax=0,
						 bool verbose=false);
  bool                          is_a_T_tessellation(bool verbose=false,
						    std::ostream& out=std::clog);
  void                          write(std::ostream&);
  void                          read(std::istream&);

private:
  Polygon                       window;
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
/*            DÉFINITION DE LA CLASSE LineTes_halfedge                        */
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
  std::vector<Polygon>              del_faces;    ///< list of deleted faces
  std::vector<Polygon>              add_faces;    ///< list of new faces
  std::vector<std::vector<Point2> > del_segs;     ///< list of deleted segments
  std::vector<std::vector<Point2> > add_segs;     ///< list of new segments
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
   * A split is the division of a face by a line segment.
   *
   * A split is represented by a halfedge bounding the face to be
   * split, a relative position (between 0 and 1) along that halfedge
   * and an angle greater than zero and smaller than pi. The splitting
   * segment starts at the point lying on the halfedge at the given
   * relative position and makes the given angle with the halfedge.
   *
   * A split has several effects on tessellation elements:
   * - Two edges bounding the split face are split.
   * - Two vertices along these two edges are created.
   * - The split face is replaced by two faces.
   * - A new non-blocking segment is created.
   */
  class Split : public Modification {
  public:
    Split();
    Split(Halfedge_handle, double=.5,double=CGAL_PI/2);
    virtual                ModList modified_elements();
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
    bool                   is_valid();
  private:
    Halfedge_handle         e1; // Split edge
    Halfedge_handle         e2; // Split edge
    Point2                  p1; // New point on e1
    Point2                  p2; // New point on e2
  };
  /** \brief A merge that can be applied to a T-tessellation
   *
   * A merge is the removal of a non-blocking segment. The two faces
   * on both sides of the removed segment are merged.
   *
   * A merge is represented by the (half)edge to be removed.
   *
   * A merge has several effects on tessellation elements:
   * - Two vertices are removed.
   * - Two pairs of edges are merged.
   * - Two faces are replaced by a single one.
   * - A non-blocking segment is removed.
   * - The status (blocking or not) of segments incident to 
   *   the removed segment may change.
   */
  class Merge : public Modification { 
  public:
    Merge();
    Merge(Halfedge_handle);
    virtual                ModList modified_elements();
    /** \brief Return a handle to the suppressed edge*/
    inline Halfedge_handle get_e(){return e;}
    bool                   is_valid();
  private:
    Halfedge_handle e;
    inline Halfedge_handle e1(){return e->twin()->next();}
    inline Halfedge_handle e2(){return e->next();}
    inline Point2           p1(){return e->source()->point();}
    inline Point2           p2(){return e->target()->point();}

  };
  /** \brief A flip that can be applied to a T-tessellation
   *
   * A flip is represented by the halfedge at the end of a blocking
   * segment that is removed.
   */
  class Flip : public Modification {
  public:
    Flip();
    Flip(Halfedge_handle);
    virtual                ModList modified_elements();
    /** \brief Return a handle to the suppressed edge*/
    inline Halfedge_handle get_e1(){return e1;}
    /** \brief Return a handle to the split edge*/
    inline Halfedge_handle get_e2(){return e2;}
    /** \brief Return a handle to the new vertex*/
    inline Point2          get_p2(){return p2;}
  private:
    Halfedge_handle        e1; // Suppressed edge
    Halfedge_handle        e2; // Split edge
    Point2                 p2; // New vertex
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

  TTessel();
  TTessel(LineTes&);
  virtual void         insert_window(Rectangle);
  virtual void         insert_window(Polygon&);
  Halfedge_handle      update(Split);
  Face_handle          update(Merge);
  Halfedge_handle      update(Flip);
  void                 clear();
  bool                 is_valid(bool verbose=false);
  Split                propose_split();
  Merge                propose_merge();
  Flip                 propose_flip();
  Split_list           split_sample(Size);
  Split_list           poisson_splits(double);
  Merge_list           all_merges();
  Flip_list            all_flips();
  Size                 number_of_blocking_segments();
  Size                 number_of_non_blocking_segments();
  Seg_sublist_iterator blocking_segments_begin();
  Seg_sublist_iterator blocking_segments_end();
  Seg_sublist_iterator non_blocking_segments_begin();
  Seg_sublist_iterator non_blocking_segments_end();
  /** \brief Sum of all halfedge lengths
   *
   * The sum of halfedge lengths is equal to the sum of cell perimeters.
   */
  inline double        get_total_internal_length(){return int_length;}
  Polygons             all_faces();
  void                 print_blocking_segments();
  void                 print_non_blocking_segments();
  void                 printRCALI(std::ostream&);
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
  std::vector<double (*)(Polygon,TTessel*)>             faces;///< face featuress
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
 */
class Energy{

public:

  Energy();  
  /** \brief Return the pointer to the associated TTessel object */
  inline TTessel*          get_ttessel(){return ttes;}
  /** \brief Return the current value of the energy for the associated
      T-tessellation */
  inline double            get_value(){return value;}
  /** \brief Return the theta parameter */
  inline CatVector         get_theta(){return theta;}
  /** \brief Set the theta parameter */
  inline void              set_theta(CatVector par){theta = par;}
  /** \brief Increment the current value of the energy */
  inline void              add_value(double delta){value +=delta;}
  /** Associate the model with a T-tessellation */
  void                     set_ttessel(TTessel*);
  double                   variation(TTessel::Modification&);
  CatVector                statistic_variation(TTessel::Modification&);
  void                     add_theta_vertices(double);
  void                     add_theta_edges(double);
  void                     add_theta_faces(double);
  void                     add_theta_segs(double);
  void                     del_theta_vertices();
  void                     del_theta_edges();
  void                     del_theta_faces();
  void                     del_theta_segs();
  void                     add_features_vertices(double (*)(Point2,TTessel*));
  void                     add_features_edges(double (*)(Segment,TTessel*));
  void                     add_features_faces(double (*)(Polygon,TTessel*) );
  void                     add_features_segs(double (*)(std::vector<Point2>,TTessel*) );
  void                     del_features_vertices();
  void                     del_features_edges();
  void                     del_features_faces();
  void                     del_features_segs();

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

/** \brief Markovian T-tessellation
 *
 * Updates are splits, merges and flips. Transitions ruled according a Metropolis-Hastings-Green construction.
 */
class SMFChain {
public:
  SMFChain() {};
  SMFChain(Energy *,double, double);
  /** \brief Set the energy function of the model to be simulated
   */
  inline void    set_energy(Energy *e) {engy = e;};
  /** \brief Get the energy function of the model to be simulated
   */
  inline Energy* get_energy() {return engy;};
  void           set_smf_prob(double,double);
  /** \enum ModType
   * \brief Types of modifications applicable to a T-tessellation*/
  enum ModType {
    SPLIT = 1,///< split (of a face)
    MERGE ,///< merge (of two faces separated by a non-blocking segment)
    FLIP ///< a segment is shortened, another one is lengthened
  };
  ModType        propose_modif_type();
  double         Hasting_ratio(TTessel::Split&,double*);
  double         Hasting_ratio(TTessel::Merge&,double*);
  double         Hasting_ratio(TTessel::Flip&,double*);
  ModifCounts    step(unsigned long int=1);
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
 * The Gibbs model to be considered is stored in a private data member as an
 * Energy object. The pseudolikelihood involves an integral over the space
 * of splits that can be applied to the observed T-tessellation. The discrete
 * approximation consists in replacing the integral by a discrete sum computed
 * on a sample of dummy splits. Dummy splits are drawn independently and
 * uniformly.
 */
class PseudoLikDiscrete {
 public:
  PseudoLikDiscrete() {};
  PseudoLikDiscrete(Energy*);
  /** \brief Get the Energy data member
   */
  inline Energy* GetEnergy() {return engy;};
  /** \brief Set the Energy data member
   */
  inline void SetEnergy(Energy *e) {engy = e;};
  void AddSplits(Size);
  void AddGivenSplit(TTessel::Split); // only for debugging?
  /** \brief Return the split term of the discrete approximation of
   * the log-pseudolikelihood*/
  inline std::vector<CatVector> GetSplitStatistics() {return stat_splits;};
  /** \brief Return the flip term of the discrete approximation of
   * the log-pseudolikelihood*/
  inline std::vector<CatVector> GetFlipStatistics() {return stat_flips;};
  void ClearSplits();
  double GetValue(CatVector);
  CatVector GetGradient(CatVector);
  CatMatrix GetHessian(CatVector);
  friend std::ostream& operator<<(std::ostream &, const PseudoLikDiscrete &);
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

/** \brief Pseudo-likelihood inference of a Gibsian T-tessellation with alternate updates
 *
 * Numerical maximization of the pseudo-likelihood. Each step consists of a 
 * double update. First, dummy splits are added to the discrete approximation
 * of the pseudo-likelihood (see class PseudoLikDiscrete). Second, the current
 * parameter estimate is updated according to Newton's method based on the
 * Hessian of the approximated pseudo-likelihood.
 */
class PLInferenceNOIS : public PseudoLikDiscrete{
 public:
  PLInferenceNOIS(Energy*,bool=false,double=1.0);
  /** \brief Return step size used in Newton's method */
  inline double GetStepSize() {return stepsize;};
  /** \brief Set step size used in Newton's method */
  inline void SetStepSize(double lambda) { stepsize = lambda;};
  void Step(unsigned int=1);
  CatVector GetEstimate();
  unsigned int Run(double=0.05, unsigned int=100);
  /** \brief Return all intermediate parameter estimates
   */
  inline std::vector<CatVector> GetEstimates() {return estimates;};
  /** \brief Return all intermediate log-pseudolikelihood values
   */
  inline std::vector<double>    GetValues() {return values;};
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
Segment                  clip_segment_by_convex_polygon(Segment, Polygon);
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
Polygon                  face2poly(TTessel::Face_handle);
/** \defgroup features Features of T-tessellations
 *
 * Functions that can be used by Energy objects for defining a Gibbs model
 * of T-tessellation. %Features are measured on vertices, edges, faces 
 * or segments.*/
double                   face_number(Polygon,TTessel*);
double                   face_area_2(Polygon,TTessel*);
double                   face_sum_of_angles(Polygon,TTessel*);
double                   face_shape(Polygon,TTessel*);
double                   angle_between_vectors(Vector v1,Vector v2);
double                   edge_length(Segment,TTessel*);
double                   is_point_inside_window(Point2, LineTes*);
double                   is_point_inside_window(Point2, TTessel*);
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
double                   min_angle(Polygon, TTessel*);
double                   sum_of_faces_squared_areas(TTessel*);
double                   sum_of_min_angles(TTessel*);
double                   sum_of_angles_obt(TTessel*);
double                   segment_size_2(std::vector<Point2>,TTessel*);
double                   sum_of_segment_squared_sizes(TTessel* t);

/*  FONCTIONS POUR CHRONOMETRER */
//inline std::clock_t      tic(void){return std::clock();}
//inline std::clock_t      toc(std::clock_t start){return std::clock()-start;}
//inline double            ticks2sec(std::clock_t nticks){return 1.0*nticks/CLOCKS_PER_SEC;}
