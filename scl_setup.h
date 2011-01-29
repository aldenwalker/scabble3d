/*
 *  scl_setup.h
 *  
 *
 *
 */

#ifndef __scl_setup__
#define __scl_setup__

#include "scabble3d.h"
#include "matrix.h"

void generate_arcs(arc** arc_list, 
                   int* num_arcs,
                   char** word_list,
                   int word_list_len);
void generate_polygons(char** word_list,
                       int word_list_len,
                       arc* arc_list,
                       int num_arcs,
                       polygon** poly_list,
                       int* num_polys,
                       int maxjun);
void create_constraint_matrix(char*** chains, 
                              int num_chains,
                              int* chain_lens,
                              arc* arc_list,
                              int num_arcs,
                              polygon* polygon_list,
                              int num_polys,
                              int* weights,
                              RatMat* constraints,
                              int** equalityType);
void scl_problem_init(scl_problem* scl_prob, 
                      char*** chains,
                      int num_chains,
                      int* chain_lens,
                      char** word_list,
                      int num_words,
                      int* weights,
                      int maxjun);
#endif
