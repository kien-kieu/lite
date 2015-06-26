/* Line Tessellation (LiTe) library
   |||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
   |||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

/* Arrangement de droites aléatoires à l'intérieur d'un carré */

// Define shorter names to please linker (g++)
#include "short_names.h"

#include <ctime>

#include <CGAL/Cartesian.h>
/* J'ai essayé plusieurs types de nombres. J'aurais aimé pouvoir
   utiliser les doubles pour des raisons de rapidité de calcul.  Mais
   le calcul d'arrangements avec des doubles semble vraiment donner
   n'importe quoi. Ce que j'ai trouvé de mieux, ce sont les flottants
   CORE. */
#include <CGAL/CORE/Expr.h>
#include <CGAL/double.h>
/*#include <CGAL/Arr_2_default_dcel.h>*/
#include <CGAL/Arr_default_dcel.h>
/* J'ai aussi envisagé d'utiliser des arrangements de polylignes. En
   effet, le bord de la fenêtre peut être considéré comme une
   polyligne. Mais c'est la seule et j'ai l'impression qu'en se
   limitant aux segments, on gagne un peu de temps.*/
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include "random_lines.h"
#include <CGAL/intersections.h>
#include <CGAL/number_utils_classes.h>

#define NB_LINES_MAX 30

typedef CORE::Expr                                      NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;

typedef Traits::Point_2                                 Point2;
typedef std::vector<Point2>                             Points;
typedef Traits::Segment_2                               Segment;
typedef Traits::Curve_2                                 Curve;
typedef Traits::X_monotone_curve_2                      X_monotone_curve;

typedef CGAL::Arr_default_dcel<Traits>                  Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr_2;

typedef CGAL::Line_2<Kernel>                            Line;
typedef std::vector<Line>                               Lines;
typedef CGAL::Iso_rectangle_2<Kernel>                   Rectangle;
typedef CGAL::Random_points_in_disc_2<Point2>           Rnd_pt_disk; 

typedef Binomial_line_process<Kernel>                   Blp;

struct noeud {
  int contrainte;
  int croise;
  int no_sur_croise; /* Comme les 2 derniers champs sont de même type,
			un seul suffit peut-être ? */
};

struct ext_contr{
  int p2;
  int d2;
  int p3;
  int d3;
};

typedef std::vector<noeud>                              droite;
typedef std::vector<droite>                             config_droites;
typedef std::vector<std::vector<int> >                  comp_seg;


// Variables globales :

int last_stl[NB_LINES_MAX][3];

class CroftonTessellation: public Arr_2 {
private:
  Rectangle window;
public:
  void insert_window(Rectangle w);
  int number_of_lines(void);
  void insert_random_lines(int nb_lines,CGAL::Random& rnd=CGAL::default_random);
  void remove_internal_lines(void);
  config_droites init_partial_stl();
}; 

void CroftonTessellation::insert_window(Rectangle w) {
  window = w;
  for (int i=0;i<4;i++) 
    insert(Curve(window[i], window[i+1]));
}

int CroftonTessellation::number_of_lines(void) {
  CroftonTessellation::Curve_iterator c1,c2;
  int nb=0;
  c1 = curve_node_begin();
  c1++;c1++;c1++;c1++;
  c2 = curve_node_end();
  while(c1++!=c2)
    nb++;
  return(nb);
}
  
void CroftonTessellation::insert_random_lines(int nb_lines,CGAL::Random& rnd) {
  Line line;
  Segment liw; // line in window (segment)

  Blp line_gen(window,nb_lines,rnd);
  Lines lines = line_gen.run();
  Lines::iterator first = lines.begin();
  Lines::iterator last = lines.end();
  //int count_lines = 4;
  while(first!=last) {
    //count_lines++;
    line = *first++;
    /* La simplification/approximation ci-dessous permet de calculer
       les arrangements plus rapidement (environ un facteur 10 pour
       une trentaine de droites). */
    line = Line(Kernel::RT(CGAL::to_double(line.a())),
		Kernel::RT(CGAL::to_double(line.b())),
		Kernel::RT(CGAL::to_double(line.c())));
    // Calculer pour chaque droite l'intersection avec le bord de la fenetre
    CGAL::Object inter = CGAL::intersection(line,window);
    if(CGAL::assign(liw,inter))
      insert(Curve(liw[0],liw[1]));
  }
}

void CroftonTessellation::remove_internal_lines() {
  /* Remove all curves except the boundary curves from the
     tessellation. BE CAUTIOUS: removing internal lines does not
     remove vertices along the boundary. */
  CroftonTessellation::Curve_iterator c1,c2;
  c1 = curve_node_begin();
  c1++;c1++;c1++;c1++;
  c2 = curve_node_end();
  while(c1!=c2)
    remove_curve(c1++);
}

config_droites CroftonTessellation::init_partial_stl() {
  droite d;
  noeud n;
  config_droites stl;
  CroftonTessellation::Curve_iterator c_it,c_it0,buf;
  CroftonTessellation::Edge_iterator e_it;
  CroftonTessellation::Halfedge_around_vertex_circulator havc,havc0;
  int i,j;

  i = 5;
  c_it0=curve_node_begin();
  c_it0++;c_it0++;c_it0++;c_it0++;
  n.contrainte = 0;
  // Parcours des droites
  for(c_it=c_it0;c_it!=curve_node_end();++c_it) {
    // Parcours des arêtes de la droite courante
    for(e_it=c_it->edges_begin();e_it!=c_it->edges_end();e_it++) {
      /* Identification de la droite qui croise la droite courante à
	 l'extrémité (début) de l'arête courante*/
      // havc0 : première demi-arête incidente
      havc0 = e_it->halfedge()->source()->incident_halfedges();
      /* On parcourt toutes les autres jusqu'à ce qu'on tombe sur une
	 demi-arête incidente contenue dans une droite distincte de la
	 droite courante */
      havc = havc0;
      while(++havc!=havc0)
	if(havc->edge_node()->curve_node()!=c_it)
	  break;
      /* Pour obtenir le numéro de la droite contenant la
	 demi-arête havc, on parcourt à l'aide de l'itérateur buf
	 l'ensemble des droites jusqu'à ce qu'on retombe sur la
	 droite considérée. */
      j = 1;
      for(buf=curve_node_begin();buf!=curve_node_end();buf++) {
	if(buf==havc->edge_node()->curve_node())
	  break;
	else
	  j++;
      }
      n.croise = j;
      d.push_back(n);
    }
    e_it--;
    havc0 = e_it->halfedge()->target()->incident_halfedges();
    havc = havc0;
    while(++havc!=havc0)
      if(havc->edge_node()->curve_node()!=c_it)
	break;
    j = 1;
    for(buf=curve_node_begin();buf!=curve_node_end();buf++) {
      if(buf==havc->edge_node()->curve_node())
	break;
      else
	j++;
    }
    n.croise = j;
    d.push_back(n);
    stl.push_back(d);
    d.clear();
  }
  // Relier les mêmes noeuds sur 2 droites distinctes
  config_droites::iterator d1 = stl.begin();
  config_droites::iterator d2 = stl.end();
  i = 5;
  while(d1!=d2) {
    droite::iterator n1 = d1->begin()+1;
    droite::iterator n2 = d1->end()-1;
    while(n1!=n2) {
      config_droites::iterator droite_croisee;
      droite_croisee = stl.begin()+n1->croise-5;
      droite::iterator na = droite_croisee->begin()+1;
      j = 1;
      while(na->croise!=i) {
	na++; j++;
      }
      n1->no_sur_croise =  j;
      n1++;
    }
    d1++; i++;
  }
  return(stl);
}


inline std::clock_t tic(void) {
  return std::clock();
}

inline std::clock_t toc(std::clock_t start) {
  return std::clock()-start;
}

inline double ticks2sec(std::clock_t nticks) {
  return 1.0*nticks/CLOCKS_PER_SEC;
}


void print_crossing_stl(config_droites &stl) {
    config_droites::iterator d1 = stl.begin();
    config_droites::iterator d2 = stl.end();
    int i = 5;
    while(d1!=d2) {
      droite::iterator n1 = d1->begin();
      droite::iterator n2 = d1->end();
      std::cout << "droite ";
      std::cout << i;
      std::cout << " : ";
      while(n1!=n2) {
	std::cout << n1->croise;
	if(n1!=d1->begin() & n1!=d1->end()-1)
	  std::cout << "(" << n1->no_sur_croise << ")";
	std::cout << " ";
	n1++;
      }
      std::cout << std::endl;
      d1++; i++;
    }
}

void afficher_contraintes(config_droites stl)
{
  int i=5;
  config_droites::iterator d1,d2;
  droite::iterator n1, n2;
  d1=stl.begin();
  d2=stl.end();
  
  while (d1!=d2) {
    n1=d1->begin();
    n2=d1->end();
    std::cout<< "droite " << i << " : ";
    while (n1!=n2){
      std::cout << n1->contrainte<< " "; 
      n1++;
    }
    std::cout<<std::endl;
    d1++; i++;
  }
}

void verifier_stl(config_droites &stl)
{
  int i;
  config_droites::iterator d1,d2;
  config_droites::iterator droite_croisee;
  droite::iterator no_curr, no_end;
  
  std::cout << "Vérification des intersections..." << std::endl;
  std::cout << "Nombre de droites : " << stl.size() <<  std::endl;
  i = 5;
  d1 = stl.begin();
  d2 = stl.end();

  while(d1!=d2) {
    std::cout << "droite " << i << " : ";
    droite::iterator n1 = d1->begin()+1;
    droite::iterator n2 = d1->end()-1;

    while(n1!=n2) {
      droite_croisee = stl.begin()+n1->croise-5;
      droite::iterator na = droite_croisee->begin()+n1->no_sur_croise;
      std::cout << na->croise<<" ";
      n1++;
    }
    std::cout << std::endl;
    d1++; i++;
    }
}

inline void insert_contrainte(config_droites& stl, droite& d, int ind_droite, int ind_no, int contrainte)
{
  int no_dc = d[ind_no].croise;
  if (no_dc > ind_droite)
    {
      int no_sur_croise = d[ind_no].no_sur_croise;
      stl[no_dc-ind_droite-1][no_sur_croise].contrainte = contrainte;
    }
}

void insert_segment(droite& d,config_droites& stl, int start, int end, int ind_droite) 
{
 
/* Algorithme 3 : 
   ajout d'un segment à une sous-tessellation partielle */ 

  int ind_no=1, last_no;
  last_no = d.size()-1;

  // On remplit le contenu des noeuds sur la droite avant "start" 
  // et on propage les 4
  while (ind_no < start)
    {
      insert_contrainte(stl,d,ind_droite,ind_no,4);
      ind_no++;

    }
 

  // On remplit "start" et on propage 3
  
  if (ind_no==start)
    {
      insert_contrainte(stl,d,ind_droite,ind_no,3);
      ind_no++;
    }


  // On remplit l'intérieur du segment et on propage 2

  while (ind_no<end)
    {
      insert_contrainte(stl,d,ind_droite,ind_no,2);
      ind_no++;
    }

  if (ind_no<last_no)
    {
      // On remplit "end" et on propage 3
      insert_contrainte(stl,d,ind_droite,ind_no,3);
      ind_no++; 
    
      // On remplit les 0 après "end" et on propage 4
      while (ind_no<last_no)
	{
	  insert_contrainte(stl,d,ind_droite,ind_no,4);
	  ind_no++;
	}
    }
}

void find_contr(droite& d,ext_contr& e)
{
  int n=d.size();

  for (int ind_no=0;ind_no<n;ind_no++)
    {
      if (d[ind_no].contrainte==2) 
	{
	  if (ind_no<e.p2) e.p2 = ind_no;
          if (ind_no>e.d2) e.d2 = ind_no;
	}

      if (d[ind_no].contrainte==3) 
	{
	  if (ind_no<e.p3) e.p3 = ind_no;
          if (ind_no>e.d3) e.d3 = ind_no;
	}
    }
      
}


comp_seg compatible_segments(droite& d)
{
  ext_contr e; 
  int n;
  droite::iterator start, end, end_droite;
  int ind_start, ind_end;
  std::vector<int> start_end;
  std::vector<std::vector<int> > segments;

  /* Algorithme 2 : 
     calcul de tous  les segments d'une droite compatibles 
     avec une sous-tessellation partielle */    
    
  start = d.begin(); ind_start = 0;
  end_droite = d.end()-1;
  
  n=d.size();
  e.p2 = n;
  e.d2 = -1;
  e.p3 = n;
  e.d3 = -1;

  find_contr(d,e); 


  /* si la droite contient à la fois des 2 et des 3
     la recherche des segments compatibles se restreint :
     -> à un segment delimité par le premier (le dernier) 2 et le bord,
        lorsque tous les trois se trouvent d'un côté de tous les 2,
     -> à un segment delimité par le premier et le dernier 2 
        lorsque tous les 3 sont situés entre deux 2 consécutifs
  */

  if (e.p3<n) // au moins une contrainte 3
    {
     if (e.p2<n)  // au moins une contrainte 2
	{ 
         // il y a des 3 separées par 2 :
	  if ((e.p3<e.p2 & e.p2<e.d3) || (e.p3<e.d2 & e.d2<e.d3)) { } 
          else 
            // tous les 3 à droite du dernier 2 :
            if (e.d2<e.p3) {start = start + e.d2; ind_start = ind_start+e.d2;}
            else
	      {
                // tous les 3 à gauche du premier 2 :   
		if (e.d3<e.p2)  {end_droite = start + e.p2;}
                // tous les 3 encadrées par deux 2 consécutifs :
                else {start = start+e.p2; ind_start = ind_start+e.p2;
		       end_droite = start + e.d2;}
	      }
          	  
	}
  
    }
  
  while (start < end_droite)
  {
    // -> pas de 3 avant "start"
    // -> pas de 4 en "start"
    if ((ind_start<e.p3) & start->contrainte != 4 )
      {
	end = start; ind_end = ind_start; 
	do
	  {
	    end++;ind_end++;
	    // -> pas de 3 après "end"
	    // -> pas de 4 en "end"  
            if ((ind_end > e.d3) & end->contrainte != 4)   
      	      {
		start_end.push_back(ind_start);
		start_end.push_back(ind_end);
		segments.push_back(start_end);
		start_end.clear();
	      }
	  }
          
        // pas de 2 entre "start" et "end"
	while (end->contrainte!=2 & end <= end_droite-1);
      }
    start++; ind_start++;
  }
return(segments);
}

unsigned long long int nb_subtessellations(config_droites &stl,int i) {
  /* Implémente Algorithme 1 : calcul du nombre de sous-tessellations
     compatibles avec une sous-tessellation parrtielle donnée. */


  /* Déterminer tous les segments de la droite numéro i compatibles
     avec la sous-tessellation partielle. */


  comp_seg all_seg_i = compatible_segments(stl[0]);
  if((int) stl.size()==1) // On est sur la dernière droite
    {    
      if (all_seg_i.size() > 0) 
	//{
	  comp_seg::iterator seg=all_seg_i.end()-1;
	  //last_stl[i][0] = (*seg)[0];
	  //last_stl[i][1] = (*seg)[1];
	  //last_stl[i][2] = stl[0].size(); 
	  //std::ofstream tracer("tracer.txt",std::ios::out); 
          //for (int j=0;j<=i;j++)
	  //  tracer << j+5 << " : (" << last_stl[j][0] << "," <<  last_stl[j][1] <<")" << " taille : " << last_stl[j][2] << std::endl;
	    //tracer.close();
	    //}
      
      return((unsigned long long int) all_seg_i.size()); 
    }
   
  unsigned long long int nb = 0; // Nombre de sous-tessellations à calculer.
  for(comp_seg::iterator seg=all_seg_i.begin();seg<all_seg_i.end();seg++) {
    config_droites new_stl(stl.begin()+1,stl.end());
    insert_segment(stl[0],new_stl,(*seg)[0],(*seg)[1],i+5); // Insertion segment    
    // Incrémentation du nbre de sous-tessellations. */
    //last_stl[i][0] = (*seg)[0];
    //last_stl[i][1] = (*seg)[1];
    //last_stl[i][2] = stl[0].size(); 
    nb += nb_subtessellations(new_stl,i+1);
  }
return(nb);
}


void usage(void) {
  std::cout << "Calcul du nombre de sous-tessellations en T associées à des" << std::endl;
  std::cout << "configurations aléatoires de droites" << std::endl;
  std::cout << std::endl;
  std::cout << "Usage : count_stl [-s s] [-d] [-r r] -n n" << std::endl;
  std::cout << std::endl;
  std::cout << "Arguments :" << std::endl;
  std::cout << "  -s s, s est la semence utilisée par le générateur aléatoire." << std::endl;
  std::cout << "   Si la semence n'est pas fournie, la date est utilisée pour initialiser ";
  std::cout << "le générateur aléatoire." << std::endl;
  std::cout << "  -r r, génère r configurations de droites. Par défaut, r=1." << std::endl;
  std::cout << "  -n n, n est le nombre de droites par configuration." << std::endl;
}

int parse_arguments(int argc, char** argv, int* nb_lines, int* nb_rep, unsigned int *seed) {
  /* Error codes:
   - 0, OK.
   - 1, number of lines is missing
   - 2, unexpected argument */
  int opt,err=0;
  extern char* optarg;
  *nb_rep = 1; // default value
  *nb_lines = -1;
  std::time_t s;
  CGAL_CLIB_STD::time( &s);
  *seed = s;
  while((opt=getopt(argc,argv,"s:dr:n:"))!=EOF) {
    switch(opt) {
    case 's':
      *seed = atoi(optarg);
      break;
    case 'r':
      *nb_rep = atoi(optarg);
      break;
    case 'n':
      *nb_lines = atoi(optarg);
      break;
    default:
      err = 2;
    }
  }
  if(*nb_lines<0)
    err = 1;
  return err;
}

int main(int argc, char** argv)
{
  CroftonTessellation  ctel;
  Rectangle            win;
  config_droites       stl;

  int                  nb_lines,nb_rep;
  unsigned int         semence;

  if(parse_arguments(argc,argv,&nb_lines,&nb_rep,&semence)!=0) {
    usage();
    return 1;
  }
  CGAL::Random our_random(semence);

  // Define window
  win = Rectangle(Point2(0,0),Point2(100,100));

  for(int j=0;j<nb_rep;j++) {
    ctel.clear();
    ctel.insert_window(win);
    stl.clear();


  // Generate a random tessellation
    ctel.insert_random_lines(nb_lines,our_random);

    std::cout << ctel.number_of_lines() << " ";
#ifdef BAVARD
    std::cout << "droites" << std::endl;
#endif
    std::cout << ctel.number_of_vertices()-2*ctel.number_of_lines()-4 << " ";
#ifdef BAVARD
    std::cout << "internal vertices" << std::endl;
#endif
    std::cout << ctel.number_of_halfedges()/2-4-2*nb_lines << " ";
#ifdef BAVARD
    std::cout << "internal edges" << std::endl;
#endif
    std::cout << ctel.number_of_faces()-1 << " ";
#ifdef BAVARD
    std::cout << "faces" << std::endl;
#endif
  
    /* Get a void representation of a partial subtessellation. */
    
    stl = ctel.init_partial_stl();
     
    
#ifdef BAVARD
    print_crossing_stl(stl);

#endif
    int nb_stl = nb_subtessellations(stl,0);
    std::cout << nb_stl;
#ifdef BAVARD
    std::cout << " sous-tessellations"  << std::endl;
#endif
    if(nb_stl==0) {
      std::cout << std::endl << "AUCUNE SOUS-TESSELLATION ! " << std::endl;
      print_crossing_stl(stl);
      break;
    }
    std::cout << std::endl;
  }


  return 0;
}
