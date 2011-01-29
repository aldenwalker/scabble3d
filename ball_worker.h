/*
 *  ball_worker.h
 *  
 *
 */

#ifndef __ball_worker__
#define __ball_worker__

#include "scabble3d.h"

void point_scl(scl_problem* scl_prob,
               rvector* point,
               mpq_t scl,
               scallop_lp_solver solver);

int min_scl_over_triangle(scl_problem* scl_prob, 
                          vert_list* V,
                          triangle* t,
                          rvector* new_vertex,
                          scallop_lp_solver solver);

int find_undone_triangle(tri_list* T, 
                         double tolerance);

void split_triangles(vert_list* V, tri_list* T, int split_ind, int new_vert);

int one_orthant_step(orthant_problem* orth, double tolerance);

void one_computation_step(ball_problem* ball);

void run_execution(execution* E);

void init_orthant_problem(orthant_problem* orth, 
                          int index, 
                          char*** chains, 
                          int* chain_lens,
                          int* weights,
                          int num_words,
                          scallop_lp_solver solver);
                          
void init_computation(execution* E,
                      chains*** chains,
                      int* chain_lens,
                      int* weights,
                      int num_words,
                      double tolerance,
                      scallop_lp_solver solver);


#endif




















