#include <vector>
#include <iostream>
#include <stdlib.h>

#include <glpk.h>

#include "lp.h"
//#include "sssgmp.h"
#include <gmp.h>
//#include "gmp/gmp-exec/include/gmp.h"



#define NOT_DOING_QSOPT

#ifndef NOT_DOING_QSOPT
extern "C" {
#include "QSopt_ex/EGlib.h"
}
extern "C" {
#include "QSopt_ex/QSopt_ex.h"
}
#endif


extern "C" {
#include "exlp-package/lpstruct.h"
}
extern "C" {
#include "exlp-package/solve_lp.h"
}
extern "C" {
#include "exlp-package/mylib.h"
}



	
	
	
//note that the scl is a rational, which doesn't currently use gmp, but
//if the answer is that weird....
void linear_program_from_ratmat( vector<polygon>& polygon_list,
                                 vector<rational>& solutionVector,
                                 rational& scl,
                                 RatMat* constraints,
                                 vector<int>& equalityType,
                                 scallop_lp_solver solver,
                                 int VERBOSE ) {
  int i,j;
  
  
  if (solver == GLPK_DOUBLE || solver == GLPK_EXACT) {   //GLPK
  
    int* ind = new int[constraints->nC-1 +1];
    double* val = new double[constraints->nC-1 +1];
    
   	glp_prob *lp;
    glp_smcp parm;
    lp = glp_create_prob();
    glp_init_smcp(&parm);
    parm.presolve=GLP_ON;
    parm.msg_lev=GLP_MSG_OFF;
    glp_set_prob_name(lp, "scl");
    glp_set_obj_dir(lp,GLP_MIN);

    glp_add_rows(lp, constraints->nR );
    mpq_t entry;
    mpq_init(entry);
    for(i=0;i<constraints->nR;i++){
      if (equalityType[i] == 0) {
        RatMat_get(constraints, i, constraints->nC-1, entry);
	      glp_set_row_bnds(lp,i+1, GLP_FX, mpq_get_d(entry), mpq_get_d(entry));
	    } else {
	      RatMat_get(constraints, i, constraints->nC-1, entry);
	      glp_set_row_bnds(lp,i+1, GLP_UP, mpq_get_d(entry), mpq_get_d(entry));
	    }
    }

    glp_add_cols(lp, constraints->nC-1);
    for(i=0;i<constraints->nC-1;i++){
	    glp_set_col_bnds(lp,i+1, GLP_LO, 0.0, 0.0);
	    glp_set_obj_coef(lp,i+1, (polygon_list[i].size-2));
    }

    for (i=0; i<constraints->nR; i++) {
      for (j=0; j<constraints->nC-1; j++) {
        RatMat_get(constraints, i, j, entry);
        ind[j+1] = j+1;
        val[j+1] = mpq_get_d(entry);
      }
      glp_set_mat_row(lp, i+1, constraints->nC-1, ind, val);
    }
    glp_simplex(lp,&parm);

    //////////////////////  this line might not compile with older glpk ////////
    if (solver == GLPK_EXACT) {
      //glp_exact(lp, &parm);  //this is disabled 
    }
    ///////////////////////////////////////////////////////////////////////////

    scl = approxRat(glp_get_obj_val(lp)/4.0);	

    for (i=0; i<constraints->nC-1; i++) {
      solutionVector[i] = approxRat(glp_get_col_prim(lp,i+1));
    }	

    glp_delete_prob(lp);

    delete[] ind;
    delete[] val;
  
  
  
  
  
  } else if (solver == QSOPT_EXACT) {     ///QSopt_exact 


  #ifndef NOT_DOING_QSOPT    /// this will only be compiled probably never
     
    int rval = 0;
    int status = 0;
    int nCols = constraints->nC-1;
    int nRows = constraints->nR;
    int* cmatcnt = new int[nCols];  //[3] = { 2, 2, 1 };
    int* cmatbeg = new int[nCols];  //[3] = { 0, 2, 4 };
    int* cmatind = new int[nCols * nRows]; // [5] = { 0, 1, 0, 1, 0 };
    char* sense = new char[nRows];  //[2] = { 'L', 'E' };
    //const char *colnames[3] = { "x", "y", "z" };
    //const char *rownames[2] = { "c1", "c2"};
    mpq_t* cmatval = new mpq_t[nCols * nRows];
    mpq_t* obj = new mpq_t[nCols];
    mpq_t* rhs = new mpq_t[nRows];
    mpq_t* lower = new mpq_t[nCols];
    mpq_t* upper = new mpq_t[nCols]; 
    mpq_t temp;
    mpq_t* x = new mpq_t[nCols];
    
    
    mpq_QSprob p;
    ///char **colnames = (char **) NULL;
    //char **rownames = (char **) NULL;
   
 
    for (i = 0; i<nCols; i++) mpq_init (obj[i]);
    for (i = 0; i<nCols; i++) 
      for (j=0; j<nRows; j++)
        mpq_init (cmatval[i*nRows + j]);
    for (i = 0; i < nRows; i++)  mpq_init(rhs[i]);
    for (i = 0; i < nCols; i++) {
      mpq_init(lower[i]);
      mpq_init(upper[i]);
    }
    mpq_t value;
    mpq_init(value);
    
    mpq_init(temp);

    for (i = 0; i<nCols; i++) {
      cmatcnt[i] = nRows;
      cmatbeg[i] = i*nRows;
      for (j=0; j<nRows; j++) {
        RatMat_get(constraints, j,i, temp);
        mpq_set(cmatval[i*nRows + j], temp);
        cmatind[i*nRows + j] = j;
      }
    }
    
    for (i = 0; i<nCols; i++) {
      mpq_set_si(obj[i], polygon_list[i].size-2, 1);
    }

    for (i = 0; i < nRows; i++) {
      RatMat_get(constraints, i, nCols, temp);
      mpq_set(rhs[i], temp);
      sense[i] = (equalityType[i] == 0 ? 'E' : 'L');
    }
   
    for (i = 0; i < nCols; i++) {
      mpq_set_si(lower[i], 0,1);
      mpq_set(upper[i], mpq_ILL_MAXDOUBLE);
    }

    
  

    /*  CPXcopylpwnames  */
    cout << "starting load\n";
    p = mpq_QSload_prob ("small", nCols, nRows, cmatcnt, cmatbeg, cmatind, cmatval,
                      QS_MIN, obj, rhs, sense, lower, upper, NULL,
                      NULL);
    cout << "Done -- starting linear programming\n";
    rval = QSexact_solver (p, NULL, NULL, NULL, DUAL_SIMPLEX, &status);
    cout << "Done\n";
    
    
    if (rval) {
        cerr <<  "QSexact_solver failed\n";
    }
    if (status != QS_LP_OPTIMAL) {
        cerr << "Did not find an optimal solution.\n";
        cerr << "Status code: " << status <<"\n";
    }
    
    
    /*  CPXgetobjval  */
    rval = mpq_QSget_objval (p, &value);
    if (rval) {
      cerr << "Could not get obj value, error code" << rval << "\n";;
    } else {
      scl = rational( mpz_get_si(mpq_numref(value)), 4*mpz_get_si(mpq_denref(value)));
    }
    
    

    for (i = 0; i < nCols; i++) mpq_init (x[i]);
    rval = mpq_QSget_x_array (p, x);
    if (rval) {
        cerr << "Could not get x-vector, error code " << rval << "\n";;
    } else {
      for (i=0; i<nCols; i++) {
        solutionVector[i] = rational( mpz_get_si(mpq_numref(x[i])), mpz_get_si(mpq_denref(x[i])));
      }
    }
    
    
    mpq_QSfree_prob (p);
    
    
    //for (i = 0; i < nCols; i++) mpq_clear (x[i]);
    //for (i = 0; i<nCols*nRows; i++) mpq_clear(cmatval[i]);
    
    //delete[] cmatcnt;
    //delete[] cmatbeg;
    //delete[] cmatind;
    //delete[] x;

      
   #endif    
      
  } else if (solver == EXLP) {
    
    LP* lp;
    int  result;
    char buf[100];
    int varNum;
    int nCols = constraints->nC-1;
    int nRows = constraints->nR;

    if (VERBOSE==1) 
      cout << "About to create a new lp\n";    
    lp = new_lp(NULL);
    
    if (VERBOSE==1) 
      cout << "Done\n";
    
    if (VERBOSE)
      cout << "Init hash\n";
        
    lp_hash_str_init(lp, lp->hash_entries);
    //my_hash_mpq_init(lp->hash_entries);
    
    if (VERBOSE)
      cout << "Done\n";

    
    
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
    
    if (VERBOSE==1) {
      cout << "Entered the row names and equalities\n";
    }
    
    //add the objective function row
    sprintf(buf, "obj");
    lp_set_obj_name(lp, buf);
    
    int* columnIndices = new int[nCols]; //this is probably useless
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
    
    if (VERBOSE==1) {
      cout << "Added the columns\n";
    }
    
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
      mpq_set_si(entry, polygon_list[i].size-2, 1);
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
    
    if (VERBOSE==1) {
      cout << "Rows: " << lp->rows;
      cout << "Vars: " << lp->vars;
    }
    
    
    result = solve_lp(lp);
      
    if (result != LP_RESULT_OPTIMAL) {
      cout << "got error code " << result << "\n";
    }
    
    lp_get_object_value(lp, &entry);
    if (mpz_get_si(mpq_denref(entry)) > 1000000) {
      cout << "I got a huge denominator?\n";
    }
    
    scl = rational(entry)/rational(4,1);
    
    for (i=0; i<nCols; i++) {
      mpq_set(entry, *vector_get_element_ptr(lp->x, columnIndices[i]));
      solutionVector[i] = rational(entry);
    }
    
    lp_free(lp);
    
  }    

}
	


void init_lp(){

  #ifndef NOT_DOING_QSOPT
	QSexactStart();
  /*  Must call QSexact_set_precision before any other QS function       */
  /*  Cannot use mpq_ILL_MAXDOUBLE or mpq_ILL_MINDOUBLE before this call */
  QSexact_set_precision (128); 
  #endif
  
  //for exlp
  mylib_init();
  
}
	
	
	
	
	
	
	
	
	
	
