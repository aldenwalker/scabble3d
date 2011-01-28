/*
 *  lp.c
 *  
 *
 *
 */


#include "exlp-package/lpstruct.h"
#include "exlp-package/solve_lph"
#include "exlp-package/mylib.h"

#include "lp.h"

void linear_program_from_ratmat(polygon* poly_list,
                                rvector* solution_vector,
                                mpq_t scl,
                                RatMat* constraints,
                                int* equality_type,
                                scallop_lp_solver solver) {
  
  int i,j;
  int dim = constraints->nC;
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
      lp_set_row_equality(lp, lp_get_row_num(lp, buf), (equalityType[i] == 0 ? 'E' : 'L'));
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
      mpq_set(solutionVector[i], entry);
    }
    
    lp_free(lp);
    
    free(columnIndices);
    mpq_clear(entry);
  
  } else {
    printf("I don't know that solver\n");
  }

  
}


