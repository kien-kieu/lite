/* Line Tessellation (LiTe) library
|||Development version
   Authors: Katarzyna Adamczyk and Kiên Kiêu.
|||Copyright INRA 2006-yyyy.
   Interdeposit Certification: IDDN.FR.001.030007.000.R.P.2015.000.31235
   License: GPL v3. */

 # include "ttessel.h"
 # include <algorithm>
 
typedef std::vector<Segment>                              Segments;  
typedef std::vector<Polygon>                              Polygons;

struct TesselList {
  Points  vertices;
  Segments   edges;
  Polygons   faces;
  std::vector<Points>   segs;
}; 

/*--------------------------*
 *    TYPE CONVERSION		*
 *--------------------------*/

TesselList TTessel2TesselList(TTessel &tesl){
 
/* Extracts vertices, edges, faces and segments of tesl into a structure "TesselList" */	

TesselList elements;
     
    Points pp;  
	for (TTessel::Vertex_iterator i=tesl.vertices_begin();i!=tesl.vertices_end();i++){
		 	elements.vertices.push_back(i->point());
	} //vertices
	
	
	for (TTessel::Edge_iterator e=tesl.edges_begin();e!=tesl.edges_end();e++){
		elements.edges.push_back(Segment(e->source()->point(),e->target()->point()));
        elements.edges.push_back(Segment(e->target()->point(),e->source()->point()));
	} // edges 
		
	for (TTessel::Face_iterator f=tesl.faces_begin();f!=tesl.faces_end();f++){
		if (!f->is_unbounded()){ 
			  Polygon poly;
			  TTessel::Ccb_halfedge_circulator e;
			  e =f->outer_ccb();
              TTessel::Halfedge_handle e_first=e;
              do {
                  poly.push_back(e->target()->point());
              e++;
              }
            while (e!=e_first);
			elements.faces.push_back(poly);}
	} // faces
	
	for (TTessel::Seg_list_iterator si=tesl.segments_begin();si!=tesl.segments_end();si++) {
   		elements.segs.push_back((*si)->list_of_points());	
	} // segments
return elements;
} 



TesselList ModList2TesselList(ModList el,bool added){
	/* extracts the added (resp. deleted) elements from ModList to TesselList */
TesselList ex_elements;

	if (added){
		
		ex_elements.vertices = el.add_vertices;	
		ex_elements.edges = el.add_edges;
		ex_elements.faces = el.add_faces;
		ex_elements.segs = el.add_segs;
	}
	else{
		ex_elements.vertices = el.del_vertices;
		ex_elements.edges = el.del_edges;
		ex_elements.faces = el.del_faces;
		ex_elements.segs = el.del_segs;	
	}
return ex_elements;	
}

void clear_list(TesselList &elements){
	/* delets TesselList elements */
elements.vertices.clear();
elements.edges.clear();
elements.faces.clear();
elements.segs.clear();
}

/*--------------------------*
 *      PRINTING COMMANDS   *
 *--------------------------*/


inline void print_hf(TTessel::Halfedge_handle e){
/* prints halfedge */

std::cout << "[(" << e->source()->point()<< ")] -[(" << e->target()->point() << ")]"<< std::endl;
}

void print_points(std::vector<Points> vopp){
/* prints a vector of points */

	for (unsigned int i=0;i!=vopp.size();i++){
	std::cout << i << ": " ;
	std::cout << "[";	
	for (unsigned int j=0;j!=(vopp[i].size()-1);j++){
	std::cout << "(" << vopp[i][j] << ")";
	std::cout <<",";
	}
	std::cout << "(" << vopp[i][vopp[i].size()-1] << ")";
	std::cout << "]"<< std::endl;
	}
}


void print_elements(TesselList elements){
/* prints  elements of TesselList*/

std::vector<Points> vopp;
	
	if (elements.vertices.size()>0){
	std::cout << "Vertices:" << std::endl;
	std::cout << std::endl;
	vopp.push_back(elements.vertices);
	print_points(vopp);
	std::cout << std::endl;
	vopp.clear();
	}
	if (elements.edges.size()>0){
	std::cout << "Edges:" << std::endl;
	std::cout << std::endl;
	Points pp;
    for (unsigned int i=0;i!=elements.edges.size();i++){
			for (unsigned int j=0;j!=2;j++)
     		pp.push_back(elements.edges[i][j]);
     		vopp.push_back(pp);
	 		pp.clear();
		}
	print_points(vopp);
	std::cout << std::endl;
	vopp.clear();
	}

	if (elements.faces.size()>0){
	std::cout << "Faces:" << std::endl;
	std::cout << std::endl;
	Points pp;
		for (unsigned int i=0;i!=elements.faces.size();i++){
			for (int j=0;j!=elements.faces[i].size();j++)
     		pp.push_back(elements.faces[i][j]);
 	 		vopp.push_back(pp);
	 		pp.clear();
		}
	print_points(vopp);
    std::cout << std::endl;
    vopp.clear();
	}
	if (elements.segs.size()>0){
    std::cout << "Segments:" << std::endl;
	std::cout << std::endl;
	print_points(elements.segs);
	}
}


/*-------------------------------------*
 *      CHECKUP OF TTESSEL ATTRIBUTES   *
 *-------------------------------------*/

/*  NEXT_HF , PREV_HF */


bool check_next_and_prev(TTessel &tesl){
/* returns "true" if attributes "next_hf" and "prev_hf" are correctly updated for all the halfedges */	

for (TTessel::Halfedge_iterator e=tesl.halfedges_begin();e!=tesl.halfedges_end();e++){
 		 	if (e->next_hf!=compute_next_hf(e)){
 		 	std::cout << "Next_hf not correctly updated for a halfedge:"<< std::endl;
 		 	print_hf(e);
 		 	return false;
 		 	}
 		 	if (e->prev_hf!=compute_prev_hf(e)){
 		 	std::cout << "Prev_hf not correctly updated for a halfedge:"<< std::endl;
 		 	print_hf(e);
 		 	return false;
 		 	}
 	}

 return true;	
}


/* DIR */

bool check_dir(TTessel &tesl){
/* returns "true" if the direction of every halfedge 
 * is equal to the direction of its prev_hf AND next_hf (if at least one of them exists) */
 
 	for (TTessel::Halfedge_iterator e=tesl.halfedges_begin();e!=tesl.halfedges_end();e++){
    		if (compute_next_hf(e)!=NULL_HALFEDGE_HANDLE && (e->dir != e->next_hf->dir)){
    		std::cout << "Error : different direction for next_hf of the halfedge:" <<std::endl;
    		print_hf(e);
    		return false;
 			}
 	 
    		if (compute_prev_hf(e)!=NULL_HALFEDGE_HANDLE && (e->dir != e->prev_hf->dir)){
    		std::cout << "Error : different direction for prev_hf of the halfedge:" <<std::endl;
    		print_hf(e);
    		return false;
 			}
    	}
  return true;  
}


/* LENGTH */

bool check_length(TTessel &tesl,double tolerance=1e-10){
/* returns "true" if the absolute difference between the calculated length *
 * and the attribute "length" is less than "tolerance" for every halfedge */
   
   Vector v;
   double diff=0;
    for (TTessel::Halfedge_iterator e=tesl.halfedges_begin();e!=tesl.halfedges_end();e++){
    	v = Vector(e->source()->point(),e->target()->point());
    	diff = abs(sqrt(CGAL::to_double(v.x()*v.x()+v.y()*v.y())) - e->length);
    	   	if (diff > tolerance){
    		std::cout << "Difference between 'length' and reference length greater than "
    		<< tolerance << " for the halfedge: " <<std::endl;
    		print_hf(e);
    		return false;
    	   	}
    	}
 return true;
}

/* TOTAL INTERNAL LENGTH */

double window_length(Rectangle r) {
/* returns the length of tessellation window */	
	
double win_length=0;
  		for (int i=0;i<4;i++) {
    	Vector v(r[i],r[i+1]);
    	win_length += sqrt(CGAL::to_double(v.x()*v.x()+v.y()*v.y()));
  	}
return win_length;
}

bool check_total_internal_length(TTessel &tesl,double tolerance=1e-10){
/* returns "true" if the difference between total length of internal edges 
 * and the attribute "total_internal_length" is less than "tolerance"*/

 double l=0, diff_l; 
 	for (TTessel::Halfedge_iterator e=tesl.halfedges_begin();e!=tesl.halfedges_end();e++){
 	    	if (tesl.is_on_boundary(e)==0){
    			Vector v(e->source()->point(),e->target()->point());
    			l += sqrt(CGAL::to_double(v.x()*v.x()+v.y()*v.y()));
 	    	}
    }   
 
    diff_l = abs(l- (tesl.get_total_internal_length()-window_length(tesl.get_window()))); 
    if (diff_l > tolerance) return false;
return true;	
}


bool check_blockings_1(TTessel &tesl) {
  for (TTessel::Seg_sublist_iterator si=tesl.non_blocking_segments_begin();
  si!=tesl.non_blocking_segments_end();si++) {
    TTessel::Halfedge_handle e = (*si)->get_halfedge_handle();
    if (e->prev_hf!=NULL_HALFEDGE_HANDLE || e->next_hf!=NULL_HALFEDGE_HANDLE) {
      std::clog << "Segment ";
      (*si)->print();
      std::clog << " considered as non_blocking" << std::endl;
      return false;
    }
  }
  
  for (TTessel::Seg_sublist_iterator si=tesl.blocking_segments_begin();
       si!=tesl.blocking_segments_end();si++) {
    TTessel::Halfedge_handle e = (*si)->get_halfedge_handle();
    if (e->prev_hf==NULL_HALFEDGE_HANDLE && e->next_hf==NULL_HALFEDGE_HANDLE) {
      std::clog << "Segment ";
      (*si)->print();
      std::clog << " considered as blocking" << std::endl;
      return false;
    }
  }
  
  return true;
}




/* CHECKUP OF ALL THE ATTRIBUTES*/

bool check_ttessel_attributes(TTessel &tesl, double one_edge_length_tol,
                               double all_edges_length_tol){
                               	
	/* returns "true" if the following attributes of ttessel are valid :
	 * next_hf, prev_hf, dir, length, total_internal_length, 
	 * blocking_segments, non_blocking_segments */
	 
	if  (!check_next_and_prev(tesl)) return false;	
	if  (!check_dir(tesl)) return false;	
	if  (!check_length(tesl,one_edge_length_tol)) return false;
	if  (!check_total_internal_length(tesl)) return false;
	if  (!check_blockings_1(tesl)) return false;
	
return true;
}



/*-----------------------------------------*
 *      CHECKUP OF TTESSEL MODIFICATIONS   *
 *-----------------------------------------*/

template < class T >

bool v1_in_v2(std::vector< T > v1,std::vector< T > v2){
/*  returns "true" if each element of vector v1
 *  belongs to vector v2 */
 
 for (unsigned int i=0;i!=v1.size();i++){
     bool v1i_in_v2 = false;	
    	for (unsigned int j=0;j!=v2.size();j++){
    		if (v1[i]==v2[j]) v1i_in_v2  = true;
    	}
	 if (!v1i_in_v2) return false; 
}
 return true;
}


template < class T >

bool v1_out_of_v2(std::vector< T > v1,std::vector< T > v2){
/* returns "true" if each element of vector v1 
 * is out of vector v2*/

	for (unsigned int i=0;i!=v1.size();i++){
    	 for (unsigned int j=0;j!=v2.size();j++){
    		if (v1[i]==v2[j]) return false;
    	 }
	}
	return true;
}

/*  Checkup of a feature type
 *  (vertices,edges,faces or segments) */

template < class T >

bool check_feature(std::vector< T > set_in,std::vector< T > add_el,std::vector< T > del_el,std::vector< T >  set_out,char element_name[30]){

/* returns "true" if set_out\set_in = add_el and set_in\set_out = del_el */

/* set_in : features  in the input tesselation */
/* add_el : features added to the input tessellation */
/* del_el : feature deleted from the input tessellation */
/* set_in : features in the output tesselation */
   
	if (!v1_in_v2(add_el,set_out)){
		std::cout << "Error in " << element_name << " update:" << std::endl; 
		std::cout << "Added " << element_name << " not in the output set!" << std::endl;  
		return false;
	}
	if (!v1_out_of_v2(add_el,set_in)){
		std::cout << "Error in " << element_name << " update:" << std::endl;
		std::cout << "Added " << element_name << " in the input set!" << std::endl; 	
		return false;
	}
	if (!v1_in_v2(del_el,set_in)){
		std::cout << "Error in " << element_name << " update:" << std::endl; 
		std::cout << "Deleted  " << element_name << " not in the input set!" << std::endl;  
		return false;
	}
	if (!v1_out_of_v2(del_el,set_out)){
		std::cout << "Error in " << element_name << " update:" << std::endl; 
		std::cout << "Deleted  " << element_name << " in the output set!" << std::endl; 
		return false;
	}
return true;
}


bool check_modifs(TesselList tesl_in,TesselList add_el,TesselList del_el,TesselList tesl_out){

/* returns "true" if the modifications of the features of "tesl_in" are valid */

if (!check_feature(tesl_in.vertices,add_el.vertices,del_el.vertices,tesl_out.vertices,"vertices")) return false;
if (!check_feature((tesl_in.edges),(add_el.edges),(del_el.edges),(tesl_out.edges),"edges")) return false;
if (!check_feature(tesl_in.faces,add_el.faces,del_el.faces,tesl_out.faces,"faces")) return false;
if (!check_feature(tesl_in.segs,add_el.segs,del_el.segs,tesl_out.segs,"segments")) return false;

return true;
}


/*---------*
 *  MAIN   *
 *---------*/

int main(int argc, char** argv) {
  //Point_location           pl;
TTessel                  tesl;
Rectangle                w=Rectangle(Point2(0,0),Point2(100,100));
int                      n_i=1000;
int                      seed=5;
double                   edge_length_tol=1e-10,edges_length_tol=1e-10;
TTessel::Split 			 spl;
TTessel::Merge 			 mrg;
TTessel::Flip 			 flp;
TesselList 				 tesl_in,add_el,del_el,tesl_out;

rnd = new CGAL::Random(seed); 
tesl.insert_window(w);


  for (int i=0;i<n_i;i++) {
    //std::cout <<"Iteration " <<  i << std::endl;
  	/* ------------------------ */
  	/* check of tesl attributes */
  	/* ------------------------ */
  	
  	if (!check_ttessel_attributes(tesl,edge_length_tol,edges_length_tol)){
	  std::cout << "Iteration " << i << std::endl;
  	break;
  	}
  	
  	/* -------------------------- */
  	/* check of tesl modification */
  	/* -------------------------- */
  	
	tesl_in = TTessel2TesselList(tesl);
    int type_modif = rnd->get_int(1,4);
    
   /* SPLIT */
    if (type_modif==1){
        spl = tesl.propose_split();
    	add_el = ModList2TesselList(spl.modified_elements(),true);
        del_el = ModList2TesselList(spl.modified_elements(),false);
        tesl.update(spl);
    	tesl_out = TTessel2TesselList(tesl);
  		if (!check_modifs(tesl_in,add_el,del_el,tesl_out)){
  		 	std :: cout << "Split - iteration " << i << std::endl;
  		 	std :: cout << std::endl;
    	    return 0 ;
  		 	}   
    }
    
    /* MERGE */
    
    if  ((type_modif==2) && (tesl.number_of_non_blocking_segments()>0)) {
      	 mrg = tesl.propose_merge();	
      	 add_el = ModList2TesselList(mrg.modified_elements(),true);
         del_el = ModList2TesselList(mrg.modified_elements(),false);
		 tesl.update(mrg);   
		 tesl_out = TTessel2TesselList(tesl);	
		 
		
  		 if (!check_modifs(tesl_in,add_el,del_el,tesl_out)){
  		 	std :: cout << "Merge - iteration " << i << std::endl;
  		 	std :: cout << std::endl;
  		 	return 0 ;
  		 	}     
    }
    
    
    /* FLIP */
   
    if ((type_modif==3) && (tesl.number_of_blocking_segments()>0)){ 
      		flp=tesl.propose_flip();
      		add_el = ModList2TesselList(flp.modified_elements(),true);
            del_el = ModList2TesselList(flp.modified_elements(),false);
  		 	tesl.update(flp);
  		 	tesl_out = TTessel2TesselList(tesl);
  		 	if (!check_modifs(tesl_in,add_el,del_el,tesl_out)){
  		 	std :: cout << "Flip - iteration " << i << std::endl;
  		 	std :: cout << std::endl;
  		 	return 0 ;
  		 	} 	
    }
    
    clear_list(tesl_in);
    clear_list(tesl_out);
    clear_list(add_el);
    clear_list(del_el);
     
  }
 
 tesl.display(); 
 delete rnd;

 return 0;
}

