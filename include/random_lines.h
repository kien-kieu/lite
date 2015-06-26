/* Line Tessellation (LiTe) library
   |||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
   |||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

#include <CGAL/Aff_transformation_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/generators.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/intersections.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils_classes.h>
#include <CGAL/Random.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

template<class Kernel> class IUR_line_process {
  typedef CGAL::Point_2<Kernel>                           Point;
  typedef CGAL::Vector_2<Kernel>                          Vector;
  typedef CGAL::Line_2<Kernel>                            Line;
  typedef CGAL::Iso_rectangle_2<Kernel>                   Rectangle;
  typedef std::vector<Line>                               Lines;
private:
  Rectangle rect; // Rectangle that the lines must hit
  double radius;  // Radius of the smallest circle containing the rectangle
  Vector centre;  // Centre of the smallest circle containing the rectangle
  CGAL::Random& rnd;

  inline void rect_radius(void) {
    radius = pow(CGAL::to_double((rect.xmax()-rect.xmin()))/2,2.0);
    radius += pow(CGAL::to_double((rect.ymax()-rect.ymin()))/2,2.0);
    radius = sqrt(radius);
  };
  inline void rect_centre(void) {
    centre = Vector((rect.xmax()+rect.xmin())/2,(rect.ymax()+rect.ymin())/2);
  };
public:
  // Constructor
  IUR_line_process(Rectangle r,CGAL::Random& rnd_gen=CGAL::default_random) : rect(r), rnd(rnd_gen) {
    rect_radius();
    rect_centre();
  };


  virtual ~IUR_line_process() {};
  Lines run(void); // Line generator 
  virtual int nb_lines(void) = 0;
};

template<class Kernel> 
std::vector<CGAL::Line_2<Kernel> > IUR_line_process<Kernel>::run(void) {
  typedef CGAL::Point_2<Kernel>                           Point;
  typedef CGAL::Vector_2<Kernel>                          Vector;
  typedef CGAL::Line_2<Kernel>                            Line;
  typedef std::vector<Line>                               Lines;
  typedef CGAL::Aff_transformation_2<Kernel>              Transformation;
  typedef CGAL::Random_points_in_disc_2<Point>            Rnd_pt_disk; 
  typedef typename Kernel::FT                             NT;
  typedef typename Kernel::RT                             RT;
  
  Lines            lines; // Lines hitting the rectangle
  Point            lfoot; // Line foot
  Line             rl;
  Transformation   scale;
  RT               r;
  int              i=0,n;

  Rnd_pt_disk glft(radius,rnd);
  n = nb_lines();
  r = RT(radius);
  while (i<n) {
    lfoot = *glft++;
    double sf = CGAL::to_double(CGAL::square<NT>(lfoot.x())
				+CGAL::square<NT>(lfoot.y()));
    sf = sqrt(sf);
    scale = Transformation(CGAL::SCALING,RT(sf),r);
    lfoot = scale(lfoot);
    // Prendre la droite dont le pied est egal au point genere
    // Calculer la direction orthogonale
    rl = Line(lfoot+centre,Vector(lfoot.y(),-lfoot.x()));
    // Test if the line hit the rectangle
    if (CGAL::do_intersect(rl,rect)) {
      lines.push_back(rl);
      i++;
    }
  }
  return lines;
};

template<class Kernel> class Binomial_line_process: public IUR_line_process<Kernel> {
  typedef CGAL::Iso_rectangle_2<Kernel> Rectangle;
private:
  int number_of_lines;
public:
  // Constructor
  Binomial_line_process(Rectangle r,int nl,CGAL::Random& rnd_gen) : 
    IUR_line_process<Kernel>(r,rnd_gen), number_of_lines(nl) {};
  int nb_lines(void) {
    return number_of_lines;
  };
};
    
template<class Kernel> class Poisson_line_process: public IUR_line_process<Kernel> {
  typedef CGAL::Iso_rectangle_2<Kernel> Rectangle;
private:
  double mean_nb_lines;
  gsl_rng *random_gen;
public:
  Poisson_line_process(Rectangle r, double inty, gsl_rng *rng) : 
    IUR_line_process<Kernel>(r), random_gen(rng) {
    double perimeter = 2*CGAL::to_double(r.xmax()-r.xmin()+r.ymax()-r.ymin());
    mean_nb_lines = inty*perimeter/CGAL_PI;
  }; // Constructor
  int nb_lines(void) {
    return gsl_ran_poisson(random_gen,mean_nb_lines);
  };
};

