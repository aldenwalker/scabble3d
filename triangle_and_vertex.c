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

/*****************************************************************************/
/* this is a handy cross product   (only for 3d)                             */
/*****************************************************************************/
void dvector_cross(double* dest, double* a, double* b) {
  dest[0] = a[1]*b[2] - a[2]*b[1];
  dest[1] = a[2]*b[0] - a[0]*b[2];
  dest[2] = a[0]*b[1] - a[1]*b[0];
}

double dvector_dot(double* a, double* b) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}







/*****************************************************************************/
/* triangle and triangle list functions                                      */
/*****************************************************************************/
double triangle_area(vert_list* V, triangle* T) {
  double side1[3];
  double side2[3];
  double side1_len = 0;
  double side2_len = 0;
  double dp;
  double angle = 0;
  int i;
  
  for (i=0; i<3; i++) {
    side1[i] = mpq_get_d(V->verts[ T->verts[1] ].coord[i]) -  
              mpq_get_d(V->verts[ T->verts[0] ].coord[i]);
    side2[i] = mpq_get_d(V->verts[ T->verts[2] ].coord[i]) - 
               mpq_get_d(V->verts[ T->verts[0] ].coord[i]);
    side1_len += side1[i]*side1[i];
    side2_len += side2[i]*side2[i];
  }
  side1_len = sqrt(side1_len);
  side2_len = sqrt(side2_len);
  dp = dvector_dot(side1, side2);
  angle = acos( dp / (side1_len*side2_len) );  
  
  return 0.5 * side1_len * side2_len * sin(angle);
}


void tri_list_add_copy(tri_list* T, triangle* t) {
  int i;
  T->tris = (triangle*)realloc((void*)(T->tris), (T->num_tris+1)*sizeof(triangle));
  T->tris[T->num_tris].verts = (int*)malloc(3*sizeof(int));
  for (i=0; i<3; i++) {
    T->tris[T->num_tris].verts[i] = t->verts[i];
    //rvector_copy(&T->tris[T->num_tris].normal[i], &t->normal[i]);
  }
  T->tris[T->num_tris].area = t->area;
  T->tris[T->num_tris].is_scl_linear = t->is_scl_linear;
  T->num_tris++;
}

void tri_list_delete_index(tri_list* T, int index) {
  int i;
  printf("I'm deleting the triangle at position %d\n", index);
  free(T->tris[index].verts);
  for (i = index; i<T->num_tris-1; i++) {
    T->tris[i] = T->tris[i+1];
  }
  T->num_tris--;
  T->tris = (triangle*)realloc((void*)(T->tris),
                               (T->num_tris)*sizeof(triangle));
}

void tri_list_delete_indices(tri_list* T, int ind1, int ind2) {
  int min_ind = (ind1 > ind2 ? ind2 : ind1);
  int max_ind = (ind1 > ind2 ? ind1 : ind2);
  tri_list_delete_index(T, max_ind);
  tri_list_delete_index(T, min_ind);
}

void tri_list_free(tri_list* T) {
  int i;
  for (i=0; i<T->num_tris; i++) {
    free(T->tris[i].verts);
  }
  free(T->tris);
}

void tri_list_print(tri_list* T, vert_list* V) {
  int i;
  for (i=0; i<T->num_tris; i++) {
    if (T->tris[i].is_scl_linear == 1) {
      printf("*");
    }
    printf("(area %f) [%d,%d,%d] = ", T->tris[i].area, 
                            T->tris[i].verts[0],
                            T->tris[i].verts[1],
                            T->tris[i].verts[2]);
    rvector_print(&V->verts[T->tris[i].verts[0]]);
    printf(", ");
    rvector_print(&V->verts[T->tris[i].verts[1]]);
    printf(", ");
    rvector_print(&V->verts[T->tris[i].verts[2]]);
    printf("\n");
  }
} 
  
void tri_list_print_d(tri_list* T, vert_list_d* V) {
  int i;
  for (i=0; i<T->num_tris; i++) {
    if (T->tris[i].is_scl_linear == 1) {
      printf("*");
    }
    printf("(area %f) [%d,%d,%d] = ", T->tris[i].area, 
                            T->tris[i].verts[0],
                            T->tris[i].verts[1],
                            T->tris[i].verts[2]);
    printf("[%f,%f,%f], ", V->verts[T->tris[i].verts[0]][0], 
                          V->verts[T->tris[i].verts[0]][1], 
                          V->verts[T->tris[i].verts[0]][2]);
    printf("[%f,%f,%f], ", V->verts[T->tris[i].verts[1]][0], 
                          V->verts[T->tris[i].verts[1]][1], 
                          V->verts[T->tris[i].verts[1]][2]);
    printf("[%f,%f,%f]\n", V->verts[T->tris[i].verts[2]][0], 
                          V->verts[T->tris[i].verts[2]][1], 
                          V->verts[T->tris[i].verts[2]][2]);

  }
}                         

/*****************************************************************************/
/* rational vector functions                                                 */
/*****************************************************************************/
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



/*****************************************************************************/
/* vertex list functions                                                     */
/*****************************************************************************/
void vert_list_add_copy(vert_list* V, rvector* v) {
  int i;
  printf("I'm adding a vertex\n");
  V->num_verts ++;
  V->verts = (rvector*)realloc((void*)(V->verts), (V->num_verts)*sizeof(rvector));
  rvector_init(&(V->verts[V->num_verts-1]), v->dim);
  for (i=0; i<v->dim; i++) {
    mpq_set(V->verts[V->num_verts-1].coord[i], v->coord[i]);
  }
  //printf("Created vertex:\n");
  //rvector_print(&(V->verts[V->num_verts-1]));
  //printf("\n -- current list at %x\n", (int)V);
  //for (i=0; i<V->num_verts; i++) {
  //  rvector_print(&(V->verts[i])); printf("\n");
  // }
   
}

void vert_list_d_add_copy(vert_list_d* V, double* v) {
  int i;
  V->num_verts ++;
  V->verts = (double**)realloc((void*)(V->verts), (V->num_verts)*sizeof(double*));
  V->verts[V->num_verts-1] = (double*)malloc(3*sizeof(double));
  for (i=0; i<3; i++) {
    V->verts[V->num_verts-1][i] = v[i];
  }
}
 
void vert_list_d_delete_index(vert_list_d* V, int ind) {
  int i;
  free(V->verts[ind]);
  for (i = ind; i<V->num_verts-1; i++) {
    V->verts[i] = V->verts[i+1];
  }
  V->num_verts--;
  V->verts = (double**)realloc((void*)(V->verts),
                               (V->num_verts)*sizeof(double*));
}
  
  
void vert_list_free(vert_list* V) {
  int i;
  for (i=0; i<V->num_verts; i++) {
    rvector_free(&V->verts[i]);
  }
  free(V->verts);
}

void vert_list_d_free(vert_list_d* V) {
  int i;
  for (i=0; i<V->num_verts; i++) {
    free(V->verts[i]);
  }
  free(V->verts);
}  
  
  
  
  










