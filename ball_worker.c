/*
 *  ball_worker.c
 *  
 
 these are the backend functions which compute the scl ball;
 they are intended to be run in a separate worker thread
  
 *
 */

#include "ball_worker.h"




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
  rvector_init(&temp_vert);
  mpq_init(min_scl);
  
  its_linear = min_scl_over_triangle(orth->scl_prob, 
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
  vert_list_add(orth->vertices, &temp_vert);
  
  //now, split the triangles
  split_triangles(orth->vertices,
                  orth->triangles,
                  i,
                  orth->vertices->num_verts - 1);
  
  //we are now done with this step
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







  
         
         
         
