/*
 *  triangle_and_vertex.c
 *  
 
 for clearness, these functions which deal with triangle and vertex
 lists, and triangles, etc
 
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>

#include "triangle_and_vertex.h"

double triangle_area(vert_list* V, triangle* T) {
  double base[3];
  double height[3];
  double base_len = 0;
  double height_len = 0;
  int i;
  
  for (i=0; i<3; i++) {
    base[i] = mpq_get_d(V->verts[ T->verts[1] ].coord[i]) -  
              mpq_get_d(V->verts[ T->verts[0] ].coord[i]);
    height[i] = mpq_get_d(V->verts[ T->verts[2] ].coord[i]) - (0.5*base[i]);
    
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


void tri_list_delete_index(tri_list* T, int index) {
  int i;
  free(T->tris[index].verts);
  for (i = index; i<T->num_tris-1; i++) {
    T->tris[i] = T->tris[i+1];
  }
  T->num_tris--;
  T->tris = (triangle*)realloc((void*)(T->tris),
                               (T->num_tris)*sizeof(triangle));
}

void tri_list_delete_indices(tri_list* T, int ind1, int ind2) {
  int min_ind = (ind1 > ind2 ? ind2 : ind2);
  int max_ind = (ind1 > ind2 ? ind1 : ind2);
  tri_list_delete_index(T, max_ind);
  tri_list_delete_index(T, min_ind);
}

void rvector_init(rvector* v, int len) {
  v->coord = (mpq_t*)malloc(len*sizeof(mpq_t));
  int i;
  for (i=0; i<len; i++) {
    mpq_init(v->coord[i]);
  }
  v->dim = len;
  v->malloced = 1;
}

void rvector_free(rvector* v) {
  int i;
  for (i=0; i<v->dim; i++) {
    mpq_clear(v->coord[i]);
  }
  free(v->coord);
  v->malloced = 0;
}

void rvector_dot(mpq_t ans, rvector* v1, rvector* v2) {
  int i;
  mpq_t temp;
  mpq_init(temp);
  mpq_set_si(ans, 0, 1);
  for (i=0; i<v1->dim; i++) {
    mpq_mul(temp, v1->coord[i], v2->coord[i]);
    mpq_add(ans, ans, temp);
  }
  mpq_clear(temp);
}

void rvector_cross(rvector* c, rvector* v1, rvector* v2) {
  mpq_t temp;
  mpq_init(temp);
  
  mpq_mul(c->coord[0], v1->coord[1], v2->coord[2]);
  mpq_mul(temp, v1->coord[2], v2->coord[1]);
  mpq_sub(c->coord[0], c->coord[0], temp);
  
  mpq_mul(c->coord[1], v1->coord[2], v2->coord[0]);
  mpq_mul(temp, v1->coord[0], v2->coord[2]);
  mpq_sub(c->coord[1], c->coord[1], temp);
  
  mpq_mul(c->coord[2], v1->coord[0], v2->coord[1]);
  mpq_mul(temp, v1->coord[1], v2->coord[0]);
  mpq_sub(c->coord[2], c->coord[2], temp);
  
  mpq_clear(temp);
}
  

void rvector_sub(rvector* dest, rvector* a, rvector* b) {
  int i;
  for (i=0; i<dest->dim; i++) {
    mpq_sub(dest->coord[i], a->coord[i], b->coord[i]);
  }
}
 
void rvector_print(rvector* v) {
  int i;
  printf("[");
  for (i=0; i<v->dim-1; i++) {
    mpq_out_str(NULL, 10, v->coord[i]);
    printf(", ");
  }
  mpq_out_str(NULL, 10, v->coord[v->dim-1]);
  printf("]");
}

void vert_list_add_copy(vert_list* V, rvector* v) {
  int i;
  printf("I'm adding a vertex\n");
  V->num_verts ++;
  V->verts = (rvector*)realloc((void*)(V->verts), (V->num_verts)*sizeof(rvector));
  rvector_init(&(V->verts[V->num_verts-1]), v->dim);
  for (i=0; i<v->dim; i++) {
    mpq_set(V->verts[V->num_verts-1].coord[i], v->coord[i]);
  }
  printf("Created vertex:\n");
  //rvector_print(&(V->verts[V->num_verts-1]));
  //printf("\n -- current list at %x\n", (int)V);
  //for (i=0; i<V->num_verts; i++) {
  //  rvector_print(&(V->verts[i])); printf("\n");
  // }
   
}

  

  
  
                          









