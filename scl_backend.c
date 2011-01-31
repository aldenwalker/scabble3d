/*
 *  ball_worker.c
 *  
 
 these are the backend functions which compute the scl ball;
 they are intended to be run in a separate worker thread
  
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <semaphore.h>

#include <gtk/gtk.h>

#include "scabble3d.h"
#include "triangle_and_vertex.h"
#include "word.h"
#include "scl_backend.h"

#include "exlp-package/lpstruct.h"
#include "exlp-package/solve_lp.h"
#include "exlp-package/mylib.h"

/*****************************************************************************/
/* linear programming                                                        */
/*****************************************************************************/
void init_lp() {
  mylib_init();
}

void linear_program_from_ratmat(polygon* poly_list,
                                rvector* solution_vector,
                                mpq_t scl,
                                RatMat* constraints,
                                int* equality_type,
                                enum scallop_lp_solver solver) {
  
  int i,j;
  //so I'm just going to implement exlp for now
  
  if (solver == EXLP) {
    
    LP* lp;
    int  result;
    char buf[100];
    int varNum;
    int nCols = constraints->nC-1;
    int nRows = constraints->nR;
    
    //if (VERBOSE==1) 
    //  cout << "About to create a new lp\n";    
    lp = new_lp(NULL);
    
    //if (VERBOSE==1) 
    //  cout << "Done\n";
    
    //if (VERBOSE)
    //  cout << "Init hash\n";
    
    lp_hash_str_init(lp, lp->hash_entries);
    //my_hash_mpq_init(lp->hash_entries);
    
    //if (VERBOSE)
    //  cout << "Done\n";
    
    
    
    sprintf(buf, "scabble");
    lp_set_name(lp, buf);
    
    //the lp by default (set in lpstruct.c)
    //has all the right stuff, I think
    
    //now we input the rows -- it likes to name them
    for (i=0; i<constraints->nR; i++) {
      sprintf(buf, "r%d", i);
      lp_add_row(lp, buf);
      lp_set_row_equality(lp, lp_get_row_num(lp, buf), (equality_type[i] == 0 ? 'E' : 'L'));
    }
    
    //if (VERBOSE==1) {
    //  cout << "Entered the row names and equalities\n";
   // }
    
    //add the objective function row
    sprintf(buf, "obj");
    lp_set_obj_name(lp, buf);
    
    int* columnIndices = (int*)malloc(nCols*sizeof(int));//new int[nCols]; //this is probably useless
    int rowNum;
    int objectiveIndex = lp_get_row_num(lp, (char*)"obj");
    mpq_t entry;
    mpq_init(entry);
    
    for (i=0; i<nCols; i++) {
      //since we only enter entries from a column once, we only need to do this once,
      //as opposed to read_columns in mps.c
      sprintf(buf, "col%d", i);
      varNum = lp_add_var_without_A(lp, buf);
      columnIndices[i] = varNum;
      
      //for (j=0; j<nRows; j++) {
      //we're going to enter something in the varNum column, in the right row
      //with coefficient from constraints, but we will know what to do 
      //}
    }
    
    //if (VERBOSE==1) {
    //  cout << "Added the columns\n";
   // }
    
    //this is from mps.c
    matrix_resize(lp->A, lp->rows, lp->vars);
    vector_resize(lp->b, lp->rows);
    vector_resize(lp->xb, lp->rows);
    vector_resize(lp->cb, lp->rows);
    
    
    //now we actually enter the data
    for (i=0; i<nCols; i++) {
      varNum = columnIndices[i];
      for (j=0; j<nRows; j++) {
        sprintf(buf, "r%d", j);
        rowNum = lp_get_row_num(lp, buf);
        //put the entry in
        RatMat_get(constraints, j,i, entry);
        lp_set_coefficient(lp, entry, rowNum, varNum);
      }
      //set the objective function
      mpq_set_si(entry, poly_list[i].num_arcs-2, 1);
      lp_set_coefficient(lp, entry, objectiveIndex, varNum);
    }
    
    //lp->maximize is false, so we reverse the sign
    vector_rev_sgn(lp->c);
    
    //read the right hand sides
    for (i=0; i<nRows; i++) {
      sprintf(buf, "r%d", i);
      rowNum = lp_get_row_num(lp, buf);
      RatMat_get(constraints, i, nCols, entry);
      lp_set_rhs(lp, rowNum, entry);
    }
    
    //read the bounds on the columns
    mpq_set_si(entry, 0,1);
    for (i=0; i<nCols; i++) {
      varNum = columnIndices[i];
      if (lp->lower.is_valid[varNum] == FALSE) //this will always be the case I think
        mpq_init(lp->lower.bound[varNum]);
      mpq_set(lp->lower.bound[varNum], entry);
      lp->lower.is_valid[varNum] = TRUE;
    }
    
    //if (VERBOSE==1) {
    //  cout << "Rows: " << lp->rows;
    //  cout << "Vars: " << lp->vars;
    //}
    
    
    result = solve_lp(lp);
    
    if (result != LP_RESULT_OPTIMAL) {
      printf("got error code %d\n", result);
    }
    
    lp_get_object_value(lp, &entry);
    //if (mpz_get_si(mpq_denref(entry)) > 1000000) {
    //  cout << "I got a huge denominator?\n";
    //}
    
    mpq_t temp;
    mpq_init(temp);
    mpq_set_si(temp, 1, 4);
    mpq_mul(scl, entry, temp);
    mpq_clear(temp);
    
    //scl = rational(entry)/rational(4,1);
    
    for (i=0; i<nCols; i++) {
      mpq_set(entry, *vector_get_element_ptr(lp->x, columnIndices[i]));
      mpq_set(solution_vector->coord[i], entry);
    }
    
    lp_free(lp);
    
    free(columnIndices);
    mpq_clear(entry);
  
  } else {
    printf("I don't know that solver\n");
  }

  
}









/*****************************************************************************/
/* these compute the arcs and polygons                                       */
/*****************************************************************************/
void generate_arcs(arc** arc_list, 
                   int* num_arcs,
                   char** word_list,
                   int word_list_len) {
  int word1, word2, i, j;
  arc temp_arc;
  int a, b;
  int* word_lens = (int*)malloc(word_list_len*sizeof(int));
  (*arc_list) = NULL;
  (*num_arcs) = 0;
  for (i=0; i<word_list_len; i++) {
    word_lens[i] = strlen(word_list[i]);
  }
  
  for (word1=0; word1<word_list_len; word1++) {
    for (i=0; i<word_lens[word1]; i++) {
      
      a = (int)word_list[word1][i];
      for (word2=word1; word2<word_list_len; word2++) {
        for (j=(word1==word2 ? i : 0);
             j<word_lens[word2]; 
             j++) {
          b = (int)word_list[word2][j];
          if ( (32+a-b)%64 == 0 ) {     //the letters are inverses
            (*arc_list) = (arc*)realloc((void*)(*arc_list),
                                        (*num_arcs+2)*sizeof(arc));
            temp_arc.first = i;
            temp_arc.last = j;
            temp_arc.first_word = word1;
            temp_arc.last_word = word2;
            (*arc_list)[*num_arcs] = temp_arc;
            temp_arc.first = j;
            temp_arc.last = i;
            temp_arc.first_word = word2;
            temp_arc.last_word = word1;
            (*arc_list)[*num_arcs+1] = temp_arc;
            (*num_arcs) += 2;
          }
        }
      }
      
    }
  }
  
  free(word_lens);
  
}
            
  



void generate_polygons(char** word_list,
                       int word_list_len,
                       arc* arc_list,
                       int num_arcs,
                       polygon** poly_list,
                       int* num_polys,
                       int maxjun) {
  int i,j,k,add,close,appeared,l;
  int size;
  polygon testpoly;
  int* word_lens = (int*)malloc(word_list_len*sizeof(int));
  (*poly_list) = NULL;
  (*num_polys) = 0;
  for (i=0; i<word_list_len; i++) {
    word_lens[i] = strlen(word_list[i]);
  }
  
  testpoly.arc = (int*)malloc(maxjun*sizeof(int));
  testpoly.num_arcs = 0;
  
  
  for(i=0;i<num_arcs-1;i++){  // i is initial index of hypothetical polygon
    testpoly.arc[0]=i;
    j=(i+1);
    size=1;
    while(size>0){
      
      add=1;    // haven't yet added
      close=1;  // haven't yet closed
      
      // does arc j glue up to end of last arc of testpoly?
      if(arc_list[j].first_word == arc_list[testpoly.arc[size-1]].last_word && 
         (arc_list[j].first - arc_list[testpoly.arc[size-1]].last-1)%word_lens[arc_list[j].first_word]==0){    
        appeared=1;
        for(k=0;k<size;k++){
          if(j==testpoly.arc[k]){    //has arc j already appeared in testpoly?
            appeared=0;
          }
        }
        
        if(appeared==1){          // if j has not appeared
          testpoly.arc[size]=j;    // then add it to testpoly
          add=0;                  // note that we have added an arc
          
          // does it close up?
          if(arc_list[j].last_word == arc_list[testpoly.arc[0]].first_word && 
             (arc_list[testpoly.arc[0]].first-arc_list[j].last-1)%word_lens[arc_list[j].last_word]==0){    
            testpoly.num_arcs = size+1;
            (*poly_list) = (polygon*)realloc((void*)(*poly_list),
                                             (*num_polys+1)*sizeof(polygon));
            (*poly_list)[*num_polys].num_arcs = size+1;
            (*poly_list)[*num_polys].arc = (int*)malloc((size+1)*sizeof(int));
            for (l=0; l<size+1; l++) {
              (*poly_list)[*num_polys].arc[l] = testpoly.arc[l];
            }
            (*num_polys)++;
            close=0;        // note that we have closed a polygon
          }
        }
      }
      
      
      if(add==0 && close==1){
        size++;
        j=i+1;
        if(size>=maxjun){
          j=num_arcs;
        }
      } else {
        j++;    // increment j
      }
      while(j>=num_arcs){
        size--;
        j=testpoly.arc[size]+1;
      }
    }
  }
  
}

/*****************************************************************************/
/* creates the basic constraint matrix (the restrictions from the arcs, plus */
/* the constraints saying that the boundary is ag + bh + ck (or just two)    */
/* this function works for 2d *and* 3d (actually, for any d)                 */
/*****************************************************************************/
void create_constraint_matrix(char*** chains, 
                              int num_chains,
                              int* chain_lens,
                              arc* arc_list,
                              int num_arcs,
                              polygon* polygon_list,
                              int num_polys,
                              int* weights,
                              RatMat* constraints,
                              int** equalityType) {
  int totalNumWords = 0;
  int offset;  // as we build the matrix, we use this for sanity
  int i,j,k,l;
  int n;
  int myWordNumber, firstWordNumber;
  //vector<int> equalityType; //-1 is <=, 0 is ==
  int coefFirstWord;
  int coefMyWord;
  
  for (i=0; i<num_chains; i++) {
    totalNumWords += chain_lens[i];
  }
  
  
  /*
   build the matrix of constraints for the unit ball:
   each column is a polygon, saying how much of it we want
   rows come from:
   - every arc, saying that that arc must appear the same number of times
   as positive and negative
   - the chains, saying that the boundary decomposes as ag + bh
   The last column is the b vector, as in Ax <= b is the polyhedron
   */
  RatMat_init(constraints,  
              num_arcs/2 + totalNumWords - num_chains, 
              num_polys + 1);
  (*equalityType) = (int*)malloc((num_arcs/2 + totalNumWords - num_chains)*sizeof(int));
  
  RatMat_reset_int(constraints, 0);
  
  /*  this section removed since I believe glpk will do it for us in lp.cc
   //each polygon must appear >= 0 times:
   for (i=0; i<polygon_list_length; i++) {
   RatMat_set_int(constraints, i,i, -1,1);
   RatMat_set_int(constraints, i, polygon_list_length, 0,1);
   equalityType[i] = -1;
   }
   */
  
  offset = 0;
  
  //each arc must appear the same number of times coming and going
  for (i=0; i<num_arcs/2; i++) {  //for each arc
    for (j=0; j<num_polys; j++) {  //for each polygon
      //count the number of times the arc appears in the polygon, with signs
      n=0;
      for (k=0; k<polygon_list[j].num_arcs; k++) {  // for each arc in the polygon
        if(polygon_list[j].arc[k]==2*i){
          n=n+1;
        }
        if(polygon_list[j].arc[k]==2*i+1){  //the same arc, but backwards
          n=n-1;
        }
      }
      // n is now the coefficient in the offset+i row and jth column
      RatMat_set_int(constraints, offset+i, j, n, 1);
    }
    (*equalityType)[offset+i] = 0;
    // note that these rows are == 0, which the last entry in the jth row already is
  }
  
  offset += num_arcs/2;
  
  //each chain must appear as one -- that if we have aw + bv, they must appear
  //as caw + cbv = c(aw+bv) (if the weights are 1, then just the same # of times)
  
  firstWordNumber = 0; //note this is a running total through all words
  for (i=0; i<num_chains; i++) { //for each chain
    for (j=1; j<chain_lens[i]; j++) {  //for each word past the first one
      myWordNumber = firstWordNumber + j;
      for (k=0; k<num_polys; k++) { //for each polygon
        //the number of times our word appears, divided by its weight is the same
        //as the number of times that the first word appears, divided by its weight
        coefFirstWord = 0;
        coefMyWord = 0;
        for (l=0; l<polygon_list[k].num_arcs; l++) {
          if (arc_list[polygon_list[k].arc[l]].first_word == myWordNumber &&
              arc_list[polygon_list[k].arc[l]].first == 0 ) { //for each copy of the word, this is 1
            coefMyWord ++;
          }
          else if (arc_list[polygon_list[k].arc[l]].first_word == firstWordNumber &&
                   arc_list[polygon_list[k].arc[l]].first == 0) {
            coefFirstWord++;
          }
        }
        //the entry to put in is mycoef/myweight - firstcoef/firstweight
        RatMat_set_int(constraints, offset, k, 
                       weights[firstWordNumber]*coefMyWord - weights[myWordNumber]*coefFirstWord,
                       weights[myWordNumber]*weights[firstWordNumber]);
      }
      (*equalityType)[offset] = 0;
      //again, the right hand side is 0, so we don't need to set it
      offset++;
    }
    firstWordNumber += chain_lens[i];
  }  
  
}



/*****************************************************************************/
/* compute the scl of a particular point in the span of the three (or        */
/* however many) chains                                                      */
/*****************************************************************************/
void point_scl(scl_problem* scl_prob,
               rvector* point,
               mpq_t scl,
               enum scallop_lp_solver solver) {
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
                                                  (scl_prob->constraints->nR) * sizeof(int));
  
  for (i=0; i<scl_prob->num_polys; i++) {
    for (j=0; j<3; j++) {
      mpq_set_si(coef[j], 0, 1);
    }
    for (j=0; j<scl_prob->poly_list[i].num_arcs; j++) {
      first_word_index = 0;
      for (k=0; k<3; k++) {
        if (scl_prob->arc_list[scl_prob->poly_list[i].arc[j]].first_word == first_word_index
            && scl_prob->arc_list[scl_prob->poly_list[i].arc[j]].first == 0) {
          mpq_set_si(temp_mpq, 1, scl_prob->weights[first_word_index]);
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
               scl_prob->constraints->nR-(1+i),
               scl_prob->constraints->nC-1,
               point->coord[i]);
    scl_prob->equality_type[scl_prob->constraints->nR-(1+i)] = 0;
  }
  
  //do the linear program
  rvector solution_vector;
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
                          mpq_t scl,
                          rvector* new_vertex,
                          enum scallop_lp_solver solver) {
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
          mpq_set_si(temp_mpq, scl_prob->weights[first_word_index], 1);
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
  rvector_init(&parallel_vector, 3);
  mpq_init(parallel_value1);
  mpq_init(parallel_value2);
  
  for (i=0; i<3; i++) {
    
    RatMat_change_num_rows(scl_prob->constraints, scl_prob->constraints->nR+2);
    
    rvector_sub(&spanning_vector1, 
                &(V->verts[t->verts[(i+1)%3]]), 
                &(V->verts[t->verts[i]]));
    rvector_cross(&parallel_vector, &normal_vector, &spanning_vector1);
    rvector_dot(parallel_value1, &parallel_vector, &(V->verts[t->verts[i]]));
    rvector_dot(parallel_value2, &parallel_vector, &(V->verts[t->verts[(i+2)%3]]));
    
    if (mpq_cmp(parallel_value1, parallel_value2) > 0) {
      mpq_swap(parallel_value1, parallel_value2);
    }
    
    for (j=0; j<scl_prob->num_polys; j++) {
      mpq_set_si(coef, 0, 1);
      for (k=0; k<scl_prob->poly_list[j].num_arcs; k++) {
        first_word_index = 0;
        for (l=0; l<scl_prob->num_chains; l++) {
          if (scl_prob->arc_list[scl_prob->poly_list[j].arc[k]].first_word == first_word_index
              && scl_prob->arc_list[scl_prob->poly_list[j].arc[k]].first == 0) {
            mpq_set_si(temp_mpq, scl_prob->weights[first_word_index],1);
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
          mpq_set_si(temp_mpq, scl_prob->weights[first_word_index],1);
          mpq_div(temp_mpq, solution_vector.coord[i], temp_mpq);
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
  rvector_free(&spanning_vector1);
  rvector_free(&spanning_vector2);
  rvector_free(&normal_vector);
  rvector_free(&solution_vector);
  
  if (mpq_cmp_si(scl, 1, 1)==0) {
    return 1;
  } else {
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
        && T->tris[i].area > max_found_area) {
      max_found_area = T->tris[i].area;
      max_found_index = i;
    }
  }
  return max_found_index;
}

/*****************************************************************************/
/* find the edge (or 3 for interior) in which the point is                   */
/*  (new_vert is an index in the vert list V                                 */
/*****************************************************************************/
int find_edge_for_vertex(vert_list* V, triangle* t, int new_vert) {
  int i;
  rvector cp;
  rvector_init(&cp, 3);
    
  mpq_t dp;
  mpq_init(dp);
  
  for (i=0; i<3; i++) {
    rvector_cross(&cp, &(V->verts[t->verts[i]]), &(V->verts[t->verts[(i+1)%3]]));
    rvector_dot(dp, &(V->verts[new_vert]), &cp);
    if (0 == mpq_cmp_si(dp, 0, 1)) {
      mpq_clear(dp);
      rvector_free(&cp);
      return i;
    }
  }
  mpq_clear(dp);
  rvector_free(&cp);
  return 3;
}
    

/*****************************************************************************/
/* split triangle at index split_ind with new vertex new_vert -- note that   */
/* we might be splitting an edge, in which case we need to look to find the  */
/* matching edge and split that too                                          */
/*****************************************************************************/
void split_triangles(vert_list* V, tri_list* T, int split_ind, int new_vert) {
  int j,k;
  int vert_in_edge;
  int edge_vert1, edge_vert2;
  triangle temp_tri;
  temp_tri.verts = (int*)malloc(3*sizeof(int));
  temp_tri.is_scl_linear = 0;
  
  //find out where the new vertex is in the triangle
  vert_in_edge = find_edge_for_vertex(V, &(T->tris[split_ind]), new_vert);
  
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
int one_orthant_step(orthant_problem* orth, double tolerance, enum scallop_lp_solver solver) {
  int i,j;
  int its_linear;
  triangle temp_tri;
  mpq_t min_scl;
  rvector temp_vert;
  
  printf("*Doing one orthant step\n");
  
  printf("Current triangle list:\n");
  for (i=0; i<orth->triangles->num_tris; i++) {
    printf("[%d,%d,%d] = ", orth->triangles->tris[i].verts[0],
                            orth->triangles->tris[i].verts[1],
                            orth->triangles->tris[i].verts[2]);
    rvector_print(&orth->vertices->verts[orth->triangles->tris[i].verts[0]]);
    printf(", ");
    rvector_print(&orth->vertices->verts[orth->triangles->tris[i].verts[1]]);
    printf(", ");
    rvector_print(&orth->vertices->verts[orth->triangles->tris[i].verts[2]]);
    printf("\n");
  }
    
  
  
  //find the next undone triangle
  i = find_undone_triangle(orth->triangles, tolerance);
  if (i == -1) {
    return 1;
  }
  
  printf("I found that triangle %d was undone\n", i);
  
  //so this triangle can be worked
  rvector_init(&temp_vert, 3);
  mpq_init(min_scl);
  
  its_linear = min_scl_over_triangle(orth->scl_prob, 
                                     orth->vertices,
                                     &(orth->triangles->tris[i]),
                                     min_scl,
                                     &temp_vert,
                                     solver);
  printf("min_scl: ");
  mpq_out_str(NULL, 10, min_scl);
  printf("\n");
  printf("its_linear: %d\n", its_linear);
  
  if (its_linear == 1) {
    printf("It's linear -- moving the triangle up\n");
    //move the triangle into the correct position near the top
    orth->triangles->tris[i].is_scl_linear = 1;
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
    mpq_div(temp_vert.coord[j], temp_vert.coord[j], min_scl);
  }
  
  printf("It's not linear; adding vertex ");
  rvector_print(&temp_vert);
  printf("\n");
  
  //add it to the vertex list
  vert_list_add_copy(orth->vertices, &temp_vert);
  
  //now, split the triangles
  split_triangles(orth->vertices,
                  orth->triangles,
                  i,
                  orth->vertices->num_verts - 1);
  
  i = find_undone_triangle(orth->triangles, tolerance);
  if (i==-1) {
    orth->max_undone_triangle_area = orth->triangles->tris[i].area;
  }
  
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
void one_computation_step(ball_problem* ball, enum scallop_lp_solver solver) {
  int cur_orth = ball->current_working_orthant;
  int i;
  
  printf("Making one orthant step; cur_oth = %d\n", cur_orth);
  printf("(which is orthant at %lx\n", (long int)ball->orthants[cur_orth]);
  printf("with tolerance %f\n", ball->tolerance);
  
  if (1 == one_orthant_step(ball->orthants[cur_orth], 
                            ball->tolerance, 
                            solver)) {  //if 1, we are done
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
void* run_execution(void* E_void) {
  
  execution* E = (execution*)E_void;
  
  printf("*started main execution\n");
  execution_print(E);
  
  printf("*About to wait for semaphore at %lx\n", (long int)&(E->message_sem));
  //set the status to running
  sem_wait(&(E->message_sem));
  E->status = 1;
  sem_post(&(E->message_sem));
  printf("*Done\n");
  
  sem_wait(&(E->running_sem));
  
  //enter the main loop:
  while (1) {
    //compute the next triangle
    one_computation_step(E->ball, E->solver);
    
    printf("*done computation step\n");
    
    //tell the GUI about the new triangles
    sem_wait(&(E->read_data_sem));
    g_idle_add( (GSourceFunc)update_ball_picture_while_running, E );
    printf("*g_idle_added display function\n");
    sem_wait(&(E->read_data_sem));  //wait until the gui is done
    sem_post(&(E->read_data_sem));
    
    //check for a new tolerance
    sem_wait(&(E->message_sem));
    if (E->new_tolerance_check == 1) {
      E->ball->tolerance = E->new_tolerance;
      E->new_tolerance_check = 0;
    }
    sem_post(&(E->message_sem));
    
    //check to see if we should skip the current orthant
    sem_wait(&(E->message_sem));
    if (E->skip_orthant == 1) {
      E->ball->current_working_orthant = (E->ball->current_working_orthant + 1)%4;
      E->skip_orthant = 0;
    }
    sem_post(&(E->message_sem));
    
    //check to see if we should stop
    sem_wait(&(E->message_sem));
    if (E->status_message == 1) {
      E->status = 0;
      sem_post(&(E->running_sem));
      sem_post(&(E->message_sem));
      return NULL;
    }
    
    //check to see if actually, we're just done (up to tolerance)
    if (E->ball->is_complete == 1) {
      sem_wait(&(E->message_sem));
      E->status = 0;
      sem_post(&(E->message_sem));
      sem_post(&(E->running_sem));
      //g_idle_add(/* we're done function */);
      return NULL;
    }
  }
  
}



/****************************************************************************/
/* find the arcs, polygons, and constraints, given a bunch of chains        */
/****************************************************************************/
void scl_problem_init(scl_problem* scl_prob, 
                      char*** chains,
                      int num_chains,
                      int* chain_lens,
                      char** word_list,
                      int num_words,
                      int* weights,
                      int maxjun) {
  int i,j;
  scl_prob->arc_list = NULL;
  scl_prob->poly_list = NULL;
  scl_prob->constraints = (RatMat*)malloc(sizeof(RatMat));
  scl_prob->equality_type = NULL;
  
  //copy everything
  scl_prob->chains = (char***)malloc(num_chains*sizeof(char**));
  scl_prob->num_chains = num_chains;
  scl_prob->chain_lens = (int*)malloc(num_chains*sizeof(int));
  for (i=0; i<num_chains; i++) {
    scl_prob->chain_lens[i] = chain_lens[i];
    scl_prob->chains[i] = (char**)malloc(chain_lens[i]*sizeof(char*));
    for (j=0; j<chain_lens[i]; j++) {
      scl_prob->chains[i][j] = (char*)malloc((strlen(chains[i][j])+1)*sizeof(char));
      strcpy(scl_prob->chains[i][j], chains[i][j]);
    }
  }
  
  printf("Chains:\n");
  for (i=0; i<3; i++) {
    for (j=0; j<scl_prob->chain_lens[i]; j++) {
      printf("%s ", scl_prob->chains[i][j]);
    }
    printf("\n");
  }
    
  scl_prob->word_list = (char**)malloc(num_words*sizeof(char*));
  scl_prob->num_words = num_words;
  scl_prob->weights = (int*)malloc(num_words*sizeof(int));
  for (i=0; i<num_words; i++) {
    scl_prob->word_list[i] = (char*)malloc((strlen(word_list[i])+1)*sizeof(char));
    strcpy(scl_prob->word_list[i], word_list[i]);
    scl_prob->weights[i] = weights[i];
  }
  
  printf("Word list:\n");
  for (i=0; i<scl_prob->num_words; i++) {
    printf("%d%s ", scl_prob->weights[i], scl_prob->word_list[i]);
  }
  printf("\n");
  
  
  generate_arcs(&(scl_prob->arc_list), &(scl_prob->num_arcs), 
                scl_prob->word_list, scl_prob->num_words);
  
  printf("Generated arcs\n");
  
  generate_polygons(scl_prob->word_list, scl_prob->num_words,
                    scl_prob->arc_list, scl_prob->num_arcs,
                    &(scl_prob->poly_list), &(scl_prob->num_polys),
                    maxjun);
                    
  printf("Generated polygons\n");
  
  //printf("scl_prob->equality_type: %x\n", scl_prob->equality_type);
    
  create_constraint_matrix(chains, num_chains, chain_lens,
                           scl_prob->arc_list, scl_prob->num_arcs,
                           scl_prob->poly_list, scl_prob->num_polys, 
                           scl_prob->weights, scl_prob->constraints, &(scl_prob->equality_type));           
               
  //printf("Generated constraints\n");            
  
  //printf("scl_prob->equality_type: %x\n", scl_prob->equality_type);

}


/*****************************************************************************/
/* initializes an orthant problem                                            */
/*****************************************************************************/
void orthant_problem_init(orthant_problem* orth, 
                          int index, 
                          char*** chains, 
                          int* chain_lens,
                          int* weights,
                          int num_words,
                          mpq_t* predone_scls,
                          int maxjun,
                          enum scallop_lp_solver solver) {
  int i,j;
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
  for (i=0; i<3; i++) {
    for (j=0; j<chain_lens[i]; j++) {
      all_words[index_in_words] = (char*)malloc((strlen(these_chains[i][j])+1)*sizeof(char));
      strcpy(all_words[index_in_words], these_chains[i][j]);
      index_in_words++;
    }
  }
  
  printf("Initializing orthant: %d\n", index);
  printf("with chains:\n");
  for (i=0; i<3; i++) {
    for (j=0; j<chain_lens[i]; j++) {
      printf("%s ", these_chains[i][j]);
    }
    printf("\n");
  }
  printf("And word list:\n");
  for (i=0; i<num_words; i++) {
    printf("%d%s ", weights[i], all_words[i]);
  }
  printf("\n");
  
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
  orth->max_undone_triangle_area = -1;
  
  mpq_t scl;
  mpq_init(scl);
  rvector point;
  rvector_init(&point, 3);
  for (i=0; i<3; i++) {
    mpq_set_si(point.coord[0], (i==0?1:0), 1);
    mpq_set_si(point.coord[1], (i==1?1:0), 1);
    mpq_set_si(point.coord[2], (i==2?1:0), 1);
    if (mpq_cmp_si(predone_scls[i], 0, 1) > 0) {
      mpq_set(scl, predone_scls[i]);
    } else {
      point_scl(orth->scl_prob,
                &point,
                scl,
                solver);
      mpq_set(predone_scls[i], scl);
    }
    for (j=0; j<3; j++) {
      mpq_div(point.coord[j], point.coord[j], scl);
    }
    printf("Got vertex: "); rvector_print(&point); printf("\n");
    vert_list_add_copy(orth->vertices, &point);
    
  }
  
  printf("Got initial vertices (%d):\n", orth->vertices->num_verts);
  for (i=0; i<3; i++) {
    rvector_print(&(orth->vertices->verts[i])); printf("\n");
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
void computation_init(execution* E,
                      char*** chains,
                      int* chain_lens,
                      int* weights,
                      int num_words,
                      double tolerance,
                      int maxjun,
                      enum scallop_lp_solver solver) {
  int i,j;
  
  mpq_t predone_scls[3];
  
  printf("Initializing execution\n");
  
  //set up the execution structure
  E->status = 0;
  E->status_message = 0;
  E->new_tolerance_check = 0;
  E->skip_orthant = 0;
  E->solver = solver;
  E->ball = (ball_problem*)malloc(sizeof(ball_problem));
  E->ball->tolerance = tolerance;
  E->ball->num_chains = 3;
  E->ball->chain_lens = (int*)malloc(3*sizeof(int));
  for (i=0; i<3; i++) {
    E->ball->chain_lens[i] = chain_lens[i];
    mpq_init(predone_scls[i]);
    mpq_set_si(predone_scls[i], 0, 1);
  }
  E->ball->chains = (char***)malloc(3*sizeof(char**));
  for (i=0; i<E->ball->num_chains; i++) {
    E->ball->chains[i] = (char**)malloc((E->ball->chain_lens[i])*sizeof(char*));
    for (j=0; j<E->ball->chain_lens[i]; j++) {
      E->ball->chains[i][j] = (char*)malloc((strlen(chains[i][j])+1)*sizeof(char));
      strcpy(E->ball->chains[i][j], chains[i][j]);
    }
  }
  E->ball->is_complete = 0;
  E->ball->orthants = (orthant_problem**)malloc(4*sizeof(orthant_problem*));
  for (i=0; i<4; i++) {
    E->ball->orthants[i] = (orthant_problem*)malloc(sizeof(orthant_problem));
    orthant_problem_init(E->ball->orthants[i], i, chains, chain_lens, weights, num_words, predone_scls, maxjun, solver);
  }
  E->ball->current_working_orthant = 0;
  
  sem_init(&(E->message_sem), 0, 1);
  sem_init(&(E->running_sem), 0, 1);
  sem_init(&(E->read_data_sem), 0, 1);
  
  
  for (i=0; i<3; i++) {
    mpq_clear(predone_scls[i]);
  }
}
  
  
  
  
void scl_problem_print(scl_problem* sp) {

}


void orthant_problem_print(orthant_problem* orth) {
  printf("orthant problem number %d\n", orth->orthant_num);
  printf("scl_problem at %lx\n", (long int)orth->scl_prob);
  printf("tri_list at %lx\n", (long int)orth->triangles);
  printf("vert_list at %lx\n", (long int)orth->vertices);
  printf("max triangle area: %f\n", orth->max_undone_triangle_area);
  printf("is_complete: %d\n", orth->is_complete);
}


void ball_problem_print(ball_problem* ball) {
  int i,j;
  printf("Ball problem at: %lx:\n", (long int)ball);
  printf("Chains:\n");
  for (i=0; i<ball->num_chains; i++) {
    for (j=0; j<ball->chain_lens[i]; j++) {
      printf("%s ", ball->chains[i][j]);
    }
    printf("\n");
  }
  printf("orthants at:\n");
  for (i=0; i<4; i++) {
    printf("%lx:\n", (long int)ball->orthants[i]);
    orthant_problem_print(ball->orthants[i]);
  }
  printf("\n");
  printf("is_complete: %d\n", ball->is_complete);
  printf("current_working_orthant: %d\n", ball->current_working_orthant);
  printf("tolerance: %f\n", ball->tolerance);
}
  

void execution_print(execution* E) {
  printf("Execution at %lx:\n", (long int)E);
  printf("status: %d\n", E->status);
  printf("status_message: %d\n", E->status_message);
  printf("new_tolerance_check: %d\n", E->new_tolerance_check);
  printf("skip_orthant: %d\n", E->skip_orthant);
  printf("solver: %d\n", E->solver);
  printf("ball at %lx\n", (long int)E->ball);
  
  ball_problem_print(E->ball);

}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
         
         
         
