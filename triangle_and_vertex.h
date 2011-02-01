 #ifndef __triangle_and_vertex__
 #define __triangle_and_vertex__
 
 #include <gmp.h>
 
 /*****************************************************************************/
/* this is a (3d probably) rational vector                                   */
/*****************************************************************************/
typedef struct {
  mpq_t* coord;
  int dim;
  int malloced;
} rvector;

/*****************************************************************************/
/* this is a 3d triangle                                                     */
/*****************************************************************************/
typedef struct {
  int* verts;    //this is the vertex indices in the vertex list
                 //the number of vertices should be implicit
  //rvector normal[3];
  double area;
  int is_scl_linear; //1 means yes, 0 means unknown, -1 means no
} triangle;

/*****************************************************************************/
/* this is a list  of triangles                                              */
/*****************************************************************************/
typedef struct {
  triangle* tris;
  int num_tris;
  int first_nonlinear_triangle;   //this records the index of the first triangle
                                  //which isn't linear (which we can work on)
} tri_list;


/*****************************************************************************/
/* this is a list  of vertices                                               */
/*****************************************************************************/
typedef struct {
  rvector* verts;
  int num_verts;
} vert_list;
 
typedef struct {
  double** verts;
  int num_verts;
} vert_list_d;
 
 double triangle_area(vert_list* V, triangle* T);
 void tri_list_add_copy(tri_list* T, triangle* t);
 void tri_list_print(tri_list* T, vert_list* V);
 void tri_list_print_d(tri_list* T, vert_list_d* V);
 void tri_list_delete_index(tri_list* T, int index);
 void tri_list_delete_indices(tri_list* T, int ind1, int ind2);
 void rvector_init(rvector* v, int len) ;
 void rvector_free(rvector* v);
 void rvector_sub(rvector* dest, rvector* a, rvector* b);
 void rvector_dot(mpq_t ans, rvector* v1, rvector* v2);
 void rvector_print(rvector* v);
 void rvector_cross(rvector* c, rvector* v1, rvector* v2);
 void vert_list_add_copy(vert_list* V, rvector* v);
 void vert_list_d_add_copy(vert_list_d* V, double* v);
 void vert_list_d_delete_index(vert_list_d* V, int ind);
 
 #endif
