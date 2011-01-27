/*
 *  triangle_and_vertex.c
 *  
 
 for clearness, these functions which deal with triangle and vertex
 lists, and triangles, etc
 
 *
 */


#include "scabble3d.h"   


double triangle_area(vert_list* V, triangle* T) {
  double base[3];
  double height[3];
  double base_len = 0;
  double height_len = 0;
  int i;
  
  for (i=0; i<3; i++) {
    base[i] = mpq_get_d(V->verts[ T->verts[1] ]->coord[i]) -  
              mpq_get_d(V->verts[ T->verts[0] ]->coord[i]);
    height[i] = mpq_get_d(V->verts[ T->verts[2] ]->coord[i]) - (0.5*base[i]);
    
    base_len += base[i]*base[i];
    height_len += height[i]*height[i];
  }
  
  return 0.5 * sqrt(base_len) * sqrt(height_len);
}



void tri_list_add_copy(tri_list* T, triangle* t) {
  int i;
  T->tris = (triangle*)realloc((void*)(T->tris), (T->num_tris+1)*sizeof(triangle));
  T->tris[T->num_tris].verts = (int*)malloc(3*sizeof(int));
  for (i=0; i<3; i++) {
    T->tris[T->num_tris].verts[i] = t->verts[i];
  }
  T->tris[T->num_tris].area = t->area;
  T->tris[T->num_tris].is_scl_linear = t->is_scl_linear;
}


  
  
  
  
                          









