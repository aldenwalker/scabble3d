/*
 *  ball_worker.c
 *  
 
 these are the backend functions which compute the scl ball;
 they are intended to be run in a separate worker thread
  
 *
 */

#include "ball_worker.h"





/*****************************************************************************/
/* compute the scl of a particular point in the span of the three (or        */
/* however many) chains                                                      */
/*****************************************************************************/
void point_scl(scl_problem* scl_prob,
               rvector* point,
               mpq_t scl,
               scallop_lp_solver solver) {
  int i,j,k;
  int first_word_index;
  mpq_t* coef = (mpq_t*)malloc(3*sizeof(mpq_t));
  for (i=0; i<3; i++) {
    mpq_init(coef[i]);
  }
  mpq_t temp_mpq;
  mpq_init(temp_mpq);
  
  RatMat_change_num_rows(scl_prob->constraints, 
                         scl_prob->constraints->nR+3);
  scl_prob->equality_type = (int*)realloc((void*)(scl_prob->equality_type),
                                          scl_prob->constraints->nR);
  
  for (i=0; i<scl_prob->num_polys; i++) {
    for (j=0; j<3; j++) {
      mpq_set_si(coef[j], 0, 1);
    }
    for (j=0; j<scl_prob->poly_list[i].num_arcs; j++) {
      first_word_index = 0;
      for (k=0; k<3; k++) {
        if (scl_prob->arc_list[scl_prob->poly_list[i].arc[j]].first_word == first_word_index
            && scl_prob->arc_list[scl_prob->poly_list[i].arc[j]].first == 0) {
          mpq_set_si(temp_mpq, 1, weights[first_word_index]);
          mpq_add(coef[k], coef[k], temp_mpq);
        }
        first_word_index += scl_prob->chain_lens[k];
      }
    }
    for (j=0; j<3; j++) {
      RatMat_set(scl_prob->constraints,
                 scl_prob->constraints->nR-(1+j),
                 i,
                 coef[j]);
    }
  }
  
  //set the RHS for these rows
  for (i=0; i<3; i++) {
    RatMat_set(scl_prob->constraints, 
               scl_prob->constraints->nR-(1+j),
               scl_prob->constraints->nC-1,
               point->coord[i]);
    scl_prob->equality_type[scl_prob->constraints->nR-(1+j)] = 0;
  }
  
  //do the linear program
  rvecctor solution_vector;
  rvector_init(&solution_vector, scl_prob->constraints->nC-1);
  
  linear_program_from_ratmat(scl_prob->poly_list, 
                             &solution_vector,
                             scl,
                             scl_prob->constraints,
                             scl_prob->equality_type,
                             solver);
  
  //get rid of the extra stuff we added
  RatMat_change_num_rows(scl_prob->constraints, scl_prob->constraints->nR-3);
  scl_prob->equality_type = (int*)realloc((void*)(scl_prob->equality_type),
                                          scl_prob->constraints->nR);
  mpq_clear(temp_mpq);
  for (i=0; i<3; i++) {
    mpq_clear(coef[i]);
  }
  free(coef);
  rvector_free(&solution_vector);
}
               
          
               




/*****************************************************************************/
/* compute the minimum scl over a triangle.  if scl is linear over the       */
/* triangle (so its min is 1), this returns 1.  Note this assume that the    */
/* vertices of the triangle have scl = 1                                     */
/*****************************************************************************/
int min_scl_over_triangle(scl_problem* scl_prob, 
                          vert_list* V,
                          triangle* t,
                          rvector* new_vertex,
                          scallop_lp_solver solver) {
  int i,j,k,l;
  int first_word_index;
  mpq_t temp_mpq;
  mpq_init(temp_mpq);
  mpq_t coef;
  mpq_init(coef);
  
  //we need to add the constraints for both the hyperplane and for the 
  //triangle (3 for the triangle)
  
  /************* hyperplane **************************************************/
  RatMat_change_num_rows(scl_prob->constraints, scl_prob->constraints->nR+1);
  
  rvector spanning_vector1;
  rvector spanning_vector2;
  rvector normal_vector;
  mpq_t normal_value;
  rvector_init(&spanning_vector1, 3);
  rvector_init(&spanning_vector2, 3);
  rvector_init(&normal_vector, 3);
  mpq_init(normal_value);
  
  rvector_sub(&spanning_vector1, &(V->verts[t->verts[1]]), &(V->verts[t->verts[0]]));
  rvector_sub(&spanning_vector2, &(V->verts[t->verts[2]]), &(V->verts[t->verts[0]]));
  rvector_cross(&normal_vector, &spanning_vector1, &spanning_vector2);
  rvector_dot(normal_value, &normal_vector, &(V->verts[t->verts[0]]));
  
  //now go through the polygons and say that the image in the 3d space, dotted
  //with normal_vector, gives normal_value
  for (i=0; i<scl_prob->num_polys; i++) {
    mpq_set_si(coef, 0, 1);
    for (j=0; j<scl_prob->poly_list[i].num_arcs; j++) {
      first_word_index = 0;
      for (k=0; k<scl_prob->num_chains; k++) {
        if (scl_prob->arc_list[scl_prob->poly_list[i].arc[j]].first_word == first_word_index
            && scl_prob->arc_list[scl_prob->poly_list[i].arc[j]].first == 0) {
          mpq_set_si(temp_mpq, scl_prob->weights[first_word_index]);
          mpq_div(temp_mpq, normal_vector.coord[k], temp_mpq);
          mpq_add(coef, coef, temp_mpq);
        }
        first_word_index += scl_prob->chain_lens[k];
      }
    }
    RatMat_set(scl_prob->constraints, scl_prob->constraints->nR-1, i, coef);
  }
  RatMat_set(scl_prob->constraints, 
             scl_prob->constraints->nR-1, 
             scl_prob->constraints->nC-1, 
             normal_value);
  scl_prob->equality_type = (int*)realloc((void*)(scl_prob->equality_type),
                                          scl_prob->constraints->nR*sizeof(int));
  scl_prob->equality_type[scl_prob->constraints->nR-1] = 0;
  
  
  /*********************** inequality constraints ****************************/
  //we need to cut with hyperplanes which cut out the triangle
  
  rvector parallel_vector;
  mpq_t parallel_value1;
  mpq_t parallel_value2;
  revector_init(parallel_vector);
  mpq_init(parallel_value1);
  mpq_init(parallel_value2);
  
  for (i=0; i<3; i++) {
    
    RatMat_change_num_rows(scl_prob->constraints, scl_prob->constraints->nR+2);
    
    rvector_sub(&spanning_vector1, 
                &(V->verts[t->verts[(i+1)%3]]), 
                &(V->verts[t->verts[i]]));
    rvector_cross(parallel_vector, normal_vector, spanning_vector1);
    rvector_dot(parallel_value1, parallel_vector, &(V->verts[t->verts[i]]));
    rvector_dot(parallel_value2, parallel_vector, &(V->verts[t->verts[(i+2)%3]]));
    
    if (mpq_cmp(parallel_value1, parallel_value2) > 0) {
      mpq_swap(parallel_value1, parallel_value2);
    }
    
    for (j=0; j<scl_prob->num_polys; j++) {
      mpq_set_si(coef, 0, 1);
      for (k=0; k<scl_prob->polys[j].num_arcs; k++) {
        first_word_index = 0;
        for (l=0; l<scl_prob->num_chains; l++) {
          if (scl_prob->arc_list[scl_prob->poly_list[j].arc[k]].first_word == first_word_index
              && scl_prob->arc_list[scl_prob->poly_list[j].arc[k]].first == 0) {
            mpq_set_si(temp_mpq, scl_prob->weights[first_word_index]);
            mpq_div(temp_mpq, parallel_vector.coord[l], temp_mpq);
            mpq_add(coef, coef, temp_mpq);
          }
          first_word_index += scl_prob->chain_lens[l];
        }
      }
      RatMat_set(scl_prob->constraints, scl_prob->constraints->nR-2, j, coef);
      mpq_neg(coef, coef);
      RatMat_set(scl_prob->constraints, scl_prob->constraints->nR-1, j, coef);
      scl_prob->equality_type = (int*)realloc((void*)(scl_prob->equality_type),
                                              scl_prob->constraints->nR*sizeof(int));
      scl_prob->equality_type[scl_prob->constraints->nR-2] = -1;
      scl_prob->equality_type[scl_prob->constraints->nR-1] = -1;
    }
    //set the right hand side for these two rows
    RatMat_set(scl_prob->constraints, 
               scl_prob->constraints->nR-2, 
               scl_prob->constraints->nC-1,
               parallel_value2);
    mpq_neg(parallel_value1, parallel_value1);
    RatMat_set(scl_prob->constraints, 
               scl_prob->constraints->nR-2, 
               scl_prob->constraints->nC-1,
               parallel_value1);
  }
  
  
  /**************** linear programming ***************************************/
  rvector solution_vector;
  rvector_init(&solution_vector, scl_prob->constraints->nC-1);
  mpq_t scl;
  mpq_init(scl);
  
  linear_program_from_ratmat(scl_prob->poly_list, 
                             &solution_vector,
                             scl,
                             scl_prob->constraints,
                             scl_prob->equality_type,
                             solver);
  
  
  /***************** read solution vector ************************************/
  for (i=0; i<3; i++) {
    mpq_set_si(new_vertex->coord[i], 0, 1);
  }
  
  if (0==mpq_cmp_si(scl, 1, 1)) {
    goto FREE_EVERYTHING;
  }
  
  for (i=0; i<scl_prob->num_polys; i++) {
    for (j=0; j<scl_prob->poly_list[i].num_arcs; j++) {
      first_word_index = 0;
      for (k=0; k<3; k++) {
        if (scl_prob->arc_list[scl_prob->poly_list[j].arc[k]].first_word == first_word_index
            && scl_prob->arc_list[scl_prob->poly_list[j].arc[k]].first == 0) {
          mpq_set_si(temp_mpq, scl_prob->weights[first_word_index]);
          mpq_div(temp_mpq, solution_vector->coord[i], temp_mpq);
          mpq_add(new_vertex->coord[k], new_vertex->coord[k], temp_mpq);
        }
        first_word_index += scl_prob->chain_lens[k];
      }
    }
  }
  
  //now scale so that the point has scl =1
  for (i=0; i<3; i++) {
    mpq_set_si(temp_mpq, 1, 1);
    mpq_div(temp_mpq, temp_mpq, scl);
    mpq_mul(new_vertex->coord[i], new_vertex->coord[i], temp_mpq);
  }
  
  //remove those rows and stuff
  RatMat_change_num_rows(scl_prob->constraints, 
                         scl_prob->constraints->nR-(2*3+1));
  scl_prob->equality_type = (int*)realloc((void*)(scl_prob->equality_type),
                                          scl_prob->constraints->nR*sizeof(int));
  
  
FREE_EVERYTHING:
  mpq_clear(temp_mpq);
  mpq_clear(coef);
  mpq_clear(normal_value);
  mpq_clear(parallel_value1);
  mpq_clear(parallel_value2);
  rvector_free(spanning_vector1);
  rvector_free(spanning_vector2);
  rvector_free(normal_vector);
  rvector_free(solution_vector);
  
  if (mpq_cmp_si(scl, 1, 1)==0) {
    mpq_clear(scl);
    return 1;
  } else {
    mpq_clear(scl);
    return 0;
  }

      
}


/*****************************************************************************/
/* this function finds the "worst" (highest area) undone triangle which is   */
/* larger than the tolerance                                                 */
/* it returns the index of this triangle                                     */
/* or -1 if no such triangles exists                                         */
/*****************************************************************************/
int find_undone_triangle(tri_list* T, 
                         double tolerance) {
  int i;
  int max_found_index = -1;
  double max_found_area = 0;
  for (i=T->first_nonlinear_triangle; i<T->num_tris; i++) {
    if (   T->tris[i].is_scl_linear < 1
        && T->tris[i].area > tolerance  
        && T->tris[i].area > max_area_found) {
      max_found_area = T->tris[i].area;
      max_found_index = i;
    }
  }
  return max_found_index;
}

/*****************************************************************************/
/* split triangle at index split_ind with new vertex new_vert -- note that   */
/* we might be splitting an edge, in which case we need to look to find the  */
/* matching edge and split that too                                          */
/*****************************************************************************/
void split_triangles(vert_list* V, tri_list* T, int split_ind, int new_vert) {
  int i,j,k;
  int vert_in_edge;
  int edge_vert1, edge_vert2;
  triangle temp_tri;
  temp_tri.verts = (int*)malloc(3*sizeof(int));
  temp_tri.is_scl_linear = 0;
  
  //find out where the new vertex is in the triangle
  vert_in_edge = find_edge_for_vertex(V, T->tris[split_int], new_vert);
  
  if (vert_in_edge == 3) {  //it's in the interior
    //if the triangle is 0 1 2, and the new vertex is x, then
    //the new triangles are 0 1 x, 1 2 x, 2 0 x
    for (j=0; j<3; j++) {
      temp_tri.verts[0] = T->tris[split_ind].verts[j];
      temp_tri.verts[1] = T->tris[split_ind].verts[(j+1)%3];
      temp_tri.verts[2] = new_vert;
      temp_tri.area = triangle_area(V, &temp_tri);
      tri_list_add_copy(T, &temp_tri);
    }
    //get rid of the triangle that we split
    tri_list_delete_index(T, split_ind);
    
  } else {  //it's in an edge
    //if the edge to split is i and the new vertex is x, then the new triangles
    //are i+1, i+2, x  and  i+2, i, x
    temp_tri.verts[0] = T->tris[split_ind].verts[(vert_in_edge+1)%3];
    temp_tri.verts[1] = T->tris[split_ind].verts[(vert_in_edge+2)%3];
    temp_tri.verts[2] = new_vert;
    temp_tri.area = triangle_area(V, &temp_tri);
    tri_list_add_copy(T, &temp_tri);
    temp_tri.verts[0] = T->tris[split_ind].verts[(vert_in_edge+2)%3];
    temp_tri.verts[1] = T->tris[split_ind].verts[vert_in_edge];
    temp_tri.verts[2] = new_vert;
    temp_tri.area = triangle_area(V, &temp_tri);
    tri_list_add_copy(T, &temp_tri);
    
    //now there's another triangle somewhere which has to be split as well
    edge_vert1 = T->tris[split_ind].verts[vert_in_edge];
    edge_vert2 = T->tris[split_ind].verts[(vert_in_edge+1)%3];
    for (j=T->first_nonlinear_triangle; j<T->num_tris; j++){
      for (k=0; k<3; k++) {
        if (   T->tris[j].verts[k] == edge_vert2 
            && T->tris[j].verts[(k+1)%3] == edge_vert1) {
          temp_tri.verts[0] = T->tris[j].verts[edge_vert1];
          temp_tri.verts[1] = T->tris[j].verts[(edge_vert1+1)%3];
          temp_tri.verts[2] = new_vert;
          temp_tri.area = triangle_area(V, &temp_tri);
          tri_list_add_copy(T, &temp_tri);
          temp_tri.verts[0] = T->tris[split_ind].verts[(edge_vert1+1)%3];
          temp_tri.verts[1] = T->tris[split_ind].verts[edge_vert2];
          temp_tri.verts[2] = new_vert;
          temp_tri.area = triangle_area(V, &temp_tri);
          tri_list_add_copy(T, &temp_tri);
          break;
        }
      }
      //note we can have only one triangle which matches the edge, so if
      //we broke, we found the triangle, and we can get rid of it
      if (k<3) {
        break;
      }
    }
    
    //delete both our triangle and the other one
    tri_list_delete_indices(T, split_ind, j);
    
    //free
    free(temp_tri.verts);
  }
}





/****************************************************************************/
/* do one orthant step                                                      */
/* returns 1 if this orthant is complete (up to the current tolerance)      */
/****************************************************************************/
int one_orthant_step(orthant_problem* orth, double tolerance) {
  int i,j;
  int its_linear;
  triangle temp_tri;
  mpq_t min_scl;
  rvector temp_vert;
  
  //find the next undone triangle
  i = find_undone_triangle(orth->triangles, tolerance);
  if (i == -1) {
    return 1;
  }
  
  //so this triangle can be worked
  rvector_init(&temp_vert, 3);
  mpq_init(min_scl);
  
  its_linear = min_scl_over_triangle(orth->scl_prob, 
                                     orth->vertices,
                                     &(orth->triangles->tris[i])
                                     min_scl,
                                     &temp_vert);
  
  if (its_linear == 1) {
    //move the triangle into the correct position near the top
    orth->triangles->tries[i].is_scl_linear = 1;
    temp_tri = orth->triangles->tris[orth->triangles->first_nonlinear_triangle];
    orth->triangles->tris[orth->triangles->first_nonlinear_triangle] = orth->triangles->tris[i];
    orth->triangles->tris[i] = temp_tri;
    orth->triangles->first_nonlinear_triangle ++;
    return 0;
  }
  
  //otherwise, we need to split the triangle, add new triangles, and split
  //any other triangles which share the split edges
  
  //first, scale the new vertex so that it has scl = 1
  for (j=0; j<3; j++) {
    mpq_div(temp_vert[j], temp_vert[j], min_scl);
  }
  
  //add it to the vertex list
  vert_list_add_copy(orth->vertices, &temp_vert);
  
  //now, split the triangles
  split_triangles(orth->vertices,
                  orth->triangles,
                  i,
                  orth->vertices->num_verts - 1);
  
  //we are now done with this step
  rvector_free(&temp_vert);
  mpq_clear(min_scl);
  return 0;
}




/*****************************************************************************/
/* this function takes a ball problem and just does the next triangle in the */
/* current_working_orthant which is larger than the tolerance -- if it can't */
/* find any, it goes through the orthants; if it finds something to do, it   */
/* just sets current_working_orthant to be correct, and exits                */
/* if there really is nothing to do, it sets is_complete to 1 and exits      */
/*****************************************************************************/
void one_computation_step(ball_problem* ball) {
  int cur_orth = ball->current_working_orthant;
  int i;
  
  if (1 == one_orthant_step(ball->orthants[cur_orth], ball->tolerance)) {  //if 1, we are done
    //try to find something to do
    for (i=1; i<4; i++) {
      if (0 <= find_undone_triangle(ball->orthants[(cur_orth+i)%4]->triangles, 
                                    ball->tolerance)) { //something to do here
        ball->current_working_orthant = (cur_orth+i)%4;
        return;
      }
    }
    //we found nothing to do
    ball->is_complete = 1;
    return;
    
  } else {  //this means we did some work, but there's more to do
    return;
  }
}



/*****************************************************************************/
/* take an execution and run with it -- every time a new triangle is computed*/
/* it g_idle_adds()'s a function to update the display                       */
/* it also checks for input from the GUI thread                              */
/* we expect that current_working_orthant is set, so that we know what to do */
/*****************************************************************************/
void run_execution(execution* E) {
  
  //set the status to running
  sem_wait(&(E->status));
  E->status = 1;
  sem_post(&(E->status));
  
  sem_wait(&(E->running_sem));
  
  //enter the main loop:
  while (1) {
    //compute the next triangle
    one_computation_step(&(E->ball_problem));
    
    //tell the GUI about the new triangles
    sem_wait(&(E->read_data_sem));
    g_idle_add( /*stuff*/ );
    sem_wait(&(E->read_data_sem));  //wait until the gui is done
    sem_post(&(E->read_data_sem));
    
    //check for a new tolerance
    sem_wait(&(E->new_tolerance_sem));
    if (E->new_tolerance_check == 1) {
      E->ball->tolerance = E.new_tolerance;
      E->new_tolerance_check = 0;
    }
    sem_post(&(E->new_tolerance_check));
    
    //check to see if we should skip the current orthant
    sem_wait(&(E->skip_orthant_sem));
    if (E->skip_orthant == 1) {
      E->ball->current_working_orthant = (E->ball->current_working_orthant + 1)%4;
      E->skip_orthant = 0;
    }
    
    //check to see if we should stop
    sem_wait(&(E->status_message_sem));
    if (E->status_message == 1) {
      sem_wait(&(E->status_sem));
      E->status = 0;
      sem_post(&(E->status_sem));
      sem_post(&(E->running_sem));
      return;
    }
    
    //check to see if actually, we're just done (up to tolerance)
    if (E->ball->is_complete == 1) {
      sem_wait(&(E->status_sem));
      E->status = 0;
      sem_post(&(E->status_sem));
      sem_post(&(E->running_sem));
      g_idle_add(/* we're done function */);
      return;
    }
  }
  
}





/*****************************************************************************/
/* initializes an orthant problem                                            */
/*****************************************************************************/
void init_orthant_problem(orthant_problem* orth, 
                          int index, 
                          char*** chains, 
                          int* chain_lens,
                          int* weights,
                          int num_words,
                          scallop_lp_solver solver) {
  int i,j,k;
  char*** these_chains = (char***)malloc(3*sizeof(char**));
  char** all_words = NULL;
  for (i=0; i<3; i++) {
    these_chains[i] = (char**)malloc(chain_lens[i]*sizeof(char*));
    for (j=0; j<chain_lens[i]; j++) {
      these_chains[i][j] = (char*)malloc((strlen(chains[i][j])+1)*sizeof(char));
      strcpy(these_chains[i][j], chains[i][j]);
      if (((index>>i)&1)==1) {
        invert(these_chains[i][j]);
      }
    }
  }
  int index_in_words = 0;
  all_words = (char**)malloc(num_words*sizeof(char*));
  for (i=0; i<num_chains; i++) {
    for (j=0; j<chain_lens[i]; j++) {
      all_words[idex_in_words] = (char*)malloc((strlen(chains[i][j])+1)*sizeof(char));
    }
  }
  
  orth->scl_prob = (scl_problem*)malloc(sizeof(scl_problem));
  
  scl_problem_init(orth->scl_prob, 
                   these_chains, 
                   3, 
                   chain_lens, 
                   all_words, 
                   num_words, 
                   weights, 
                   maxjun);
 
  orth->orthant_num = index;
  orth->is_complete = 0;
  orth->triangles = (tri_list*)malloc(sizeof(tri_list));
  orth->triangles->tris = (triangle*)malloc(sizeof(triangle));
  orth->triangles->first_nonlinear_triangle = 0;
  orth->triangles->num_tris = 1;
  orth->triangles->tris[0].verts = (int*)malloc(3*sizeof(int));
  orth->triangles->tris[0].verts[0] = 0;
  orth->triangles->tris[0].verts[1] = 1;
  orth->triangles->tris[0].verts[2] = 2;
  orth->vertices = (vert_list*)malloc(sizeof(vert_list));
  orth->vertices->verts = NULL;
  orth->vertices->num_verts = 0;
  
  mpq_t scl;
  mpq_init(scl);
  rvector point;
  rvector_init(&point, 3);
  for (i=0; i<3; i++) {
    mpq_set_si(point.coord[0], (i==0?1:0), 1);
    mpq_set_si(point.coord[1], (i==1?1:0), 1);
    mpq_set_si(point.coord[2], (i==2?1:0), 1);
    point_scl(orth->scl_prob,
              point,
              scl,
              solver);
    for (j=0; j<3; j++) {
      mpq_div(point.coord[j], point.coord[j], scl);
    }
    vert_list_add_copy(point);
  }
  
  orth->triangles->tris[0].area = triangle_area(orth->vertices, orth->triangles->tris);
  orth->triangles->tris[0].is_scl_linear = 0;
  
  //now we have one triangle set up
  free(all_words);
  for (i=0; i<3; i++) {
    for (j=0; j<chain_lens[i]; j++) {
      free(these_chains[i][j]);
    }
  }
  mpq_clear(scl);
  rvector_free(&point);
  
}
  
 


/*****************************************************************************/
/* this functions sets up the whole computation                              */
/* this probably won't run in a separate thread                              */
/*****************************************************************************/
void init_computation(execution* E,
                      chains*** chains,
                      int* chain_lens,
                      int* weights,
                      int num_words,
                      double tolerance,
                      scallop_lp_solver solver) {
  int i;
  //set up the execution structure
  E->status = 0;
  E->status_message = 0;
  E->new_tolerance_check = 0;
  E->skip_orthant = 0;
  E->ball = (ball_problem*)malloc(sizeof(ball_problem));
  E->ball->tolerance = tolerance;
  E->ball->num_chains = 3;
  E->ball->chain_lens = (int*)malloc(3*sizeof(int));
  for (i=0; i<3; i++) {
    E->ball->chain_lens[i] = chain_lens[i];
  }
  E->ball->is_complete = 0;
  E->orthants = (orthant_problem**)malloc(4*sizeof(orthant_problem*));
  for (i=0; i<4; i++) {
    E->orthants[i] = (orthant_problem*)malloc(sizeof(orthant_problem));
    init_orthant_problem(E->orthants[i], i, chains, chain_lens, weights, num_words, solver);
  }
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
         
         
         
