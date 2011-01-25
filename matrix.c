#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"


void Mat_init(Mat* M, int nR, int nC) {
  if (nR*nC == 0) {
    M->mat = NULL;
  } else {
    M->mat = (int*)malloc(nR*nC*sizeof(int));
  }
  M->nR = nR;
  M->nC = nC;
  //printf("Inited matrix at %x with %d rows and %d columns and matrix list at %x\n", (int)M, M->nR,M->nC, (int)(M->mat)); fflush(stdout); 
}

void Mat_set(Mat* M, int r, int c, int val) {
  (M->mat)[r*(M->nC) + c] =  val;
}

void Mat_copy_row(Mat* M, int r, Mat* source, int sourceRow) {
  int i;
  for (i=0; i<M->nC; i++) {
    (M->mat)[r*(M->nC) + i] = (source->mat)[sourceRow*(source->nC) + i];
  }
}

void Mat_copy_neg_row(Mat* M, int r, Mat* source, int sourceRow) {
  int i;
  for (i=0; i<M->nC; i++) {
    (M->mat)[r*(M->nC) + i] = -(source->mat)[sourceRow*(source->nC) + i];
  }
}




Mat* Mat_deepcopy(Mat* source) {
  int i;
  Mat* ans = (Mat*)malloc(sizeof(Mat));
  Mat_init(ans, source->nR, source->nC);
  for (i=0; i<(source->nR)*(source->nC); i++) {
    ans->mat[i] = source->mat[i];
  }
  return ans;
}



void Mat_change_num_rows(Mat* M, int nr) {
  M->mat = (int*)realloc(M->mat, (nr)*(M->nC)*sizeof(int)); 
  M->nR = nr;
}


void Mat_copy_list_row(Mat* M, int r, int* L) {
  int i;
  for (i=0; i<M->nC; i++) {
    (M->mat)[r*(M->nC) + i] = L[i];
  }
}

void Mat_append_mat(Mat* M1, Mat* M2) {
  int oldLen = (M1->nR)*(M1->nC);
  int M2Len = (M2->nR)*(M2->nC);
  Mat_change_num_rows(M1, (M1->nR)+(M2->nR));
  int i;
  for (i=0; i<M2Len; i++) {
    M1->mat[i+oldLen] = M2->mat[i];
  }
}

void Mat_reset(Mat* M, int val) {
  int i;
  for (i=0; i<(M->nR)*(M->nC); i++)
    M->mat[i] = val;
}

inline int Mat_get(Mat* M, int r, int c) {
  return M->mat[r*(M->nC) + c];
}

void Mat_increment(Mat* M, int r, int c, int incr) {
  (M->mat[r*(M->nC) + c]) += incr;
}


void Mat_print(Mat* M, int doSmall) {
  int i,j;
  for (i=0; i<M->nR; i++) {
    printf("[");
    for (j=0; j<M->nC; j++) {
      printf("%d", M->mat[i*(M->nC)+j]);
      if (j<M->nC-1)
        printf(",");
      if (doSmall!=1)
        printf(" ");
    }
    printf("]\n");
  }
}

void Mat_print_row(Mat* M, int row, int doSmall) {
  int j;
  printf("[");
  for (j=0; j<M->nC; j++) {
    printf("%d", M->mat[row*(M->nC)+j]);
    if (j<M->nC-1)
      printf(",");
    if (doSmall!=1)
      printf(" ");
  }
  printf("]\n");
}


Mat* Mat_matrix_from_rows(Mat* source, int* rows, int len) {
  int i;
  Mat* ans = (Mat*)malloc(sizeof(Mat));
  Mat_init(ans, len, source->nC);
  for (i=0; i<len; i++) {
    Mat_copy_row(ans, i, source, rows[i]);
  }
  return ans;
}


void Mat_swap_rows(Mat* M, int r1, int r2) {
  int temp, i;
  for (i=0; i<M->nC; i++) {
    temp = M->mat[r1*(M->nC)+i];
    M->mat[r1*(M->nC)+i] = M->mat[r2*(M->nC)+i];
    M->mat[r2*(M->nC)+i] = temp;
  }
}



typedef struct {
  int* rowPointer;
  int numCols;
} rowData;

int cmpRows(const void* r1, const void* r2) {
  rowData* t1 = (rowData*)r1;
  rowData* t2 = (rowData*)r2;
  int i;
  for (i=0; i<t1->numCols; i++) {
    if (t1->rowPointer[i] - t2->rowPointer[i] != 0) {
      return t1->rowPointer[i] - t2->rowPointer[i];
    }
  }
  return 0;
}
    

void Mat_sort_rows(Mat* M) {
  int i,j;
  Mat* newMat = (Mat*)malloc(sizeof(Mat));;
  rowData* rows = (rowData*)malloc((M->nR)*sizeof(rowData));
  for (i=0; i<M->nR; i++) {
    rows[i].rowPointer = &(M->mat[i*(M->nC)]);
    rows[i].numCols = M->nC;
  }
  qsort((void*)rows, M->nR, sizeof(rowData), cmpRows);
  Mat_init(newMat, M->nR, M->nC);
  for (i=0; i<M->nR; i++) {
    for (j=0; j<M->nC; j++) {
      newMat->mat[i*(newMat->nC) + j] = rows[i].rowPointer[j];
    }
  }
  free(M->mat);
  M->mat = newMat->mat;
  free(newMat);
  free(rows);
}

void Mat_remove_duplicate_rows(Mat* M) {
  int i,j;
  int numAdded=0;
  Mat* newMat = (Mat*)malloc(sizeof(Mat));
  Mat_init(newMat, 0, M->nC);
  //sort the matrix
  Mat_sort_rows(M);
  //go through and compare adjacent rows
  for (i=0; i<M->nR-1; i++) {
    for (j=0; j<M->nC; j++) {
      if (M->mat[i*(M->nC)+j] != M->mat[(i+1)*(M->nC)+j]) {
        //copy the row over
        if (numAdded >= newMat->nR) {
          Mat_change_num_rows(newMat, newMat->nR + 50);
        }
        Mat_copy_row(newMat, numAdded, M, i);
        numAdded++;
        break;
      }
    }
  }
  //copy the last row--it isn't a duplicate
  if (numAdded >= newMat->nR) {
    Mat_change_num_rows(newMat, newMat->nR +1);
  }
  Mat_copy_row(newMat, numAdded, M, (M->nR)-1);
  numAdded++;
  
  Mat_change_num_rows(newMat, numAdded);
  free(M->mat);
  M->nR = newMat->nR;
  M->nC = newMat->nC;
  M->mat = newMat->mat;
  free(newMat);
}



int Mat_rows_are_equal(Mat* M, int i, int j) {
  int k;
  for (k=0; k<M->nC; k++) {
    if (M->mat[i*(M->nC)+k] != M->mat[j*(M->nC)+k]) {
      return 0;
    }
  }
  return 1;
}



int Mat_dot_prod_with_row(Mat* M, int row, int* L) {
  int sum = 0;
  int i;
  for (i=0; i<M->nC; i++) {
    sum += M->mat[row*(M->nC)+i] * L[i];
  }
  return sum;
}
    
    
    
   
    
    
//this returns a matrix whose columns are a basis for the nullspace of M
//it converts the matrix to rational and then clears denominators, so the
//result is an integer matrix
Mat* Mat_nullspace(Mat* M) {
  int i,j;
  int* pivotCols;
  int numPivots;
  int numColsAdded;
  int itsAPivot;
  
  mpq_t temp;
  mpq_init(temp);
  mpq_t temp2;
  mpq_init(temp2);
  
  mpz_t tempZ;
  mpz_init(tempZ);
  
  RatMat* copyM = Mat_get_ratmat(M);
  RatMat* NS = (RatMat*)malloc(sizeof(RatMat));
  Mat* intNS = (Mat*)malloc(sizeof(Mat));
  
  RatMat_rref(copyM, &pivotCols, &numPivots);
  
  RatMat_init(NS, M->nC, M->nC - numPivots);
  RatMat_reset_int(NS,0);
  
  Mat_init(intNS, M->nC, M->nC - numPivots);
  
  
  //go through the columns -- if it's a pivot, ignore it
  //otherwise, put the coefficients in positions which match pivot columns
  numColsAdded=0;
  for (i=0; i<copyM->nC; i++) {
    //check if it's a pivot -- this is messy and could be improved
    itsAPivot = 0;
    for (j=0; j<numPivots; j++) { 
      if (pivotCols[j] == i) {
        itsAPivot = 1;
        break;
      }
    }
    if (itsAPivot) {
      continue;
    }
    
    for (j=0; j<numPivots; j++) {
      RatMat_get(copyM, j, i, temp);
      mpq_neg(temp, temp);
      RatMat_set(NS, pivotCols[j], numColsAdded, temp);
    }
    mpq_set_si(temp, 1, 1);
    RatMat_set(NS, i, numColsAdded, temp);
    numColsAdded++;
  }
  
  //ok now clear all the columns of M
  for (i=0; i<NS->nC; i++) {
    //get lcm
    RatMat_get(NS, 0, i, temp);
    mpz_set(tempZ, mpq_denref(temp));
    for (j=1; j<NS->nR; j++) {
      RatMat_get(NS, j, i, temp);
      mpz_lcm(tempZ, tempZ, mpq_denref(temp));
    }
    //now clear and set the entries of the integer matrix
    mpq_set_z(temp2, tempZ);
    for (j=0; j<NS->nR; j++) {
      RatMat_get(NS, j, i, temp);
      mpq_mul(temp, temp, temp2);  //the denom should be 1
      if ( 0 != mpz_cmp_si(mpq_denref(temp), 1) ) {
        printf("I should have cleared but it's not?\n");
        exit(1);
      }
      Mat_set(intNS, j, i, mpz_get_si(mpq_numref(temp)));
    }
  }
      

  RatMat_free(copyM);
  free(copyM);
  RatMat_free(NS);
  free(NS);
  
  mpq_clear(temp); mpq_clear(temp2); mpz_clear(tempZ);
  free(pivotCols);
  return intNS; 
}  
    



int Mat_rank(Mat* M) {
  int numPivots;
  int* pivotCols = NULL;
  RatMat* RM = Mat_get_ratmat(M);
  RatMat_rref(RM, &pivotCols, &numPivots);
  RatMat_free(RM);
  free(RM);
  free(pivotCols);
  return numPivots;
} 

RatMat* Mat_get_soln_to_aug(Mat* M) {
  int numPivots;
  int i;
  int* pivotCols = NULL;
  RatMat* RM = Mat_get_ratmat(M);
  RatMat_rref(RM, &pivotCols, &numPivots);
  RatMat* soln = (RatMat*)malloc(sizeof(RatMat));
  if (numPivots < M->nC-1) {
    RatMat_init(soln, M->nC-1, 0);
    return soln;
  } else {
    for (i=M->nC; i<M->nR; i++) {
      if (mpq_cmp_si(RM->mat[i*(M->nC)+(M->nC-1)], 0, 1) != 0) {
        RatMat_init(soln, M->nC-1, 0);  
        RatMat_free(RM);
        free(RM);
        free(pivotCols);
        return soln;
      }
    }
    RatMat_init(soln, M->nC-1, 1);
    for (i=0; i< M->nC-1; i++) {
      if (pivotCols[i] != i) {
        printf("I'm confused in getting solution\n");
        exit(1);
      }
      mpq_set(soln->mat[i], RM->mat[i*(M->nC)+(M->nC-1)]); // the last entry
    }  
    RatMat_free(RM);
    free(RM);
    free(pivotCols);
    return soln;
  }
  
  RatMat_free(RM);
  free(RM);
  free(pivotCols);
  return soln;
}


Mat* Mat_multiply(Mat* M1, Mat* M2) {
  int i,j,k,sum;
  if (M1->nC != M2->nR) {
    printf("Matrices have incompatible dimensions\n");
    exit(1);
  }
  Mat* ans = (Mat*)malloc(sizeof(Mat));
  Mat_init(ans, M1->nR, M2->nC);
  for (i=0; i<M1->nR; i++) {
    for (j=0; j<M2->nC; j++) {
      sum =0;
      for (k=0; k<M1->nC; k++) {
        sum += Mat_get(M1, i,k)*Mat_get(M2, k,j);
      }
      ans->mat[i*(ans->nC) + j] = sum;
    }
  }
  return ans;
}
  
  
Mat* Mat_transpose(Mat* M) {
  int i,j;
  Mat* ans = (Mat*)malloc(sizeof(Mat));
  Mat_init(ans, M->nC, M->nR);
  for (i=0; i<M->nR; i++) {
    for (j=0; j<M->nC; j++) {
      ans->mat[j*(ans->nC)+i] = M->mat[i*(M->nC)+j];
    }
  }
  return ans;
}




double Mat_dot_prod_with_float_row(Mat* M, int row, double* L) {
  double sum = 0;
  int i;
  for (i=0; i<M->nC; i++) {
    sum += (double)(M->mat[row*(M->nC)+i]) * L[i];
  }
  return sum;
}
    



void Mat_free(Mat* M) {
  free(M->mat);
}







/*****************************************/
/*****************************************/

void RatMat_init(RatMat* M, int nR, int nC) {
  int i,j;
  M->mat = (mpq_t*)malloc(nR*nC*sizeof(mpq_t));
  M->nR = nR;
  M->nC = nC;
  for (i=0; i<nR; i++) {
    for(j=0; j<nC; j++) {
      mpq_init(M->mat[i*nC+j]);
    }
  }
  //printf("Inited matrix at %x with %d rows and %d columns and matrix list at %x\n", (int)M, M->nR,M->nC, (int)(M->mat)); fflush(stdout); 
}

void RatMat_set(RatMat* M, int r, int c, mpq_t val) {
  mpq_set(M->mat[r*(M->nC) + c], val);
}


void RatMat_get(RatMat* M, int r, int c, mpq_t a) {
  mpq_set(a, M->mat[r*(M->nC) + c]);
}  
  
  
void RatMat_reset_int(RatMat* M, int val) {
  int i;
  for (i=0; i<(M->nR)*(M->nC); i++)
    mpq_set_si( M->mat[i], val, 1);
}

void RatMat_set_int(RatMat* M, int r, int c, int n, int d) {
  mpq_set_si(M->mat[r*(M->nC)+c], n, d);
}

void RatMat_free(RatMat* M) {
  int i,j;
  for (i=0; i<M->nR; i++) {
    for (j=0; j<M->nC; j++) {
      mpq_clear(M->mat[i*(M->nC) + j]);
    }
  }
  free(M->mat);
}

RatMat* Mat_get_ratmat(Mat* M) {
  int i;
  RatMat* ans = (RatMat*)malloc(sizeof(RatMat));
  RatMat_init(ans, M->nR, M->nC);
  for (i=0; i<(M->nR)*(M->nC); i++) {
    mpq_set_si(ans->mat[i], M->mat[i], 1);
  }
  return ans;
}

void RatMat_swap_rows(RatMat* M, int i, int j) {
  int k;
  int dim = M->nC;
  for (k=0; k<dim; k++) {
    mpq_swap(M->mat[i*dim+k], M->mat[j*dim+k]);
  }
}

int RatMat_rank(RatMat* M) {
  int numPivots;
  int* pivotCols = NULL;
  RatMat* RM = RatMat_deepcopy(M);
  RatMat_rref(RM, &pivotCols, &numPivots);
  RatMat_free(RM);
  free(RM);
  free(pivotCols);
  return numPivots;
} 

RatMat* RatMat_get_soln_to_aug(RatMat* M) {
  int numPivots;
  int i;
  int* pivotCols = NULL;
  RatMat* RM = RatMat_deepcopy(M);
  RatMat_rref(RM, &pivotCols, &numPivots);
  //printf("The rref is:\n");
  //RatMat_print(RM,1);
  RatMat* soln = (RatMat*)malloc(sizeof(RatMat));
  if (numPivots < M->nC-1) {
    RatMat_init(soln, M->nC-1, 0);
    return soln;
  } else {
    for (i=M->nC; i<M->nR; i++) {
      if (mpq_cmp_si(RM->mat[i*(M->nC)+(M->nC-1)], 0, 1) != 0) {
        RatMat_init(soln, M->nC-1, 0);  
        RatMat_free(RM);
        free(RM);
        free(pivotCols);
        return soln;
      }
    }
    RatMat_init(soln, M->nC-1, 1);
    for (i=0; i< M->nC-1; i++) {
      if (pivotCols[i] != i) {
        printf("I'm confused in getting solution\n");
        exit(1);
      }
      mpq_set(soln->mat[i], RM->mat[i*(M->nC)+(M->nC-1)]); // the last entry
    }  
    RatMat_free(RM);
    free(RM);
    free(pivotCols);
    return soln;
  }
  
  RatMat_free(RM);
  free(RM);
  free(pivotCols);
  return soln;
}






void RatMat_rref(RatMat* M, int** pivotCols, int* numPivots) {
  int i,j;
  int currentRowLowerLimit = 0;
  int currentColumnLowerLimit = 0;
  int nonZeroEntryInRow;
  int nonZeroEntryInCol;
  mpq_t temp;
  int dim = M->nC;
  mpq_init(temp);
  mpq_t temp2;
  mpq_init(temp2);
  (*pivotCols) = (int*)malloc(M->nC * sizeof(int));
  (*numPivots) = 0;
  
  while (currentRowLowerLimit < M->nR  && currentColumnLowerLimit < dim) {  
    //find a column in which there is a one
    nonZeroEntryInRow = -1;
    for (i=currentColumnLowerLimit; i<dim; i++) {
      for (j=currentRowLowerLimit; j<M->nR; j++) {
        if ( mpq_cmp_si(M->mat[j*dim + i], 0, 1) != 0 ) { //we found a nonzero entry
          nonZeroEntryInRow = j;
          nonZeroEntryInCol = i;
          (*pivotCols)[*numPivots] = nonZeroEntryInCol;
          (*numPivots)++;
          goto breakAll;
        }
      }
    }
    breakAll:
    if (nonZeroEntryInRow == -1) {  //we are done row reducing
      break;
    }
    //printf("current pivot position (%d, %d):\n", nonZeroEntryInRow, nonZeroEntryInCol);
    //RatMat_print(M, 0);
    
    //divide that row
    mpq_set(temp, M->mat[nonZeroEntryInRow*dim + nonZeroEntryInCol]);
    //printf("The multiplier is: ");
    //mpq_out_str(NULL, 10, temp);
    //if (mpz_cmp_si(mpq_denref(temp), 16) > 0) {
    //  exit(1);
    //}
    for (i=nonZeroEntryInCol; i<dim; i++) {   
      //printf("dividing ");
      //mpq_out_str(NULL, 10, M->mat[nonZeroEntryInRow*dim + i]);
      //printf(" by ");
      //mpq_out_str(NULL, 10, temp);
      //printf(" and got :");
      mpq_div(M->mat[nonZeroEntryInRow*dim + i], 
              M->mat[nonZeroEntryInRow*dim + i],
              temp);
      //mpq_out_str(NULL, 10, M->mat[nonZeroEntryInRow*dim + i]);
      //printf("\n");
    }
      

      
    //swap that row
    RatMat_swap_rows(M, currentRowLowerLimit, nonZeroEntryInRow);
    
    //printf("Cleared the row and swapped\n");
    //RatMat_print(M, 0);  
        
    //clear everything in that column
    for (i=0; i<M->nR; i++) {
      if (i==currentRowLowerLimit) {
        continue;
      }
      if (mpq_cmp_si(M->mat[i*dim+nonZeroEntryInCol], 0, 1)==0 ) {  //it's 0
        continue;
      }
      mpq_set(temp, M->mat[i*dim+nonZeroEntryInCol]);
      mpq_neg(temp, temp);
      //printf("clearing row %d (the factor is ", i);
      //mpq_out_str(NULL, 10, temp);
      //printf(")\n");
      for (j=nonZeroEntryInCol; j<dim; j++) {
        //printf("computing ");
        //mpq_out_str(NULL, 10, M->mat[currentRowLowerLimit*dim+j]);
        //printf(" * ");
        //mpq_out_str(NULL, 10, temp);
        //printf(" + ");
        //mpq_out_str(NULL, 10,  M->mat[i*dim+j]);
        //printf("\n");
        mpq_set(temp2, M->mat[currentRowLowerLimit*dim+j]);
        mpq_mul(temp2, temp, temp2);
        mpq_add(M->mat[i*dim+j], M->mat[i*dim+j], temp2);
        //printf("got ");
        //mpq_out_str(NULL, 10, M->mat[i*dim+j]); printf("\n");
      }
    }
    //printf("Cleared the column\n");
    //RatMat_print(M, 0);
    
    //set the new lower limits
    currentRowLowerLimit++;
    currentColumnLowerLimit = nonZeroEntryInCol + 1;
  }
  
  (*pivotCols) = (int*)realloc((void*)(*pivotCols), (*numPivots)*sizeof(int));
}  



RatMat* RatMat_inverse(RatMat* RM) {
  int i,j;
  if (RM->nR != RM->nC) {
    printf("I can't invert a nonsquare matrix\n");
    exit(1);
  }
  RatMat* RMCopy = RatMat_deepcopy(RM);  
  RatMat* Id = RatMat_identity(RM->nC);
  RatMat_append_cols_ratmat(RMCopy, Id);
  int* pivotCols;
  int numPivots;
  RatMat_rref(RMCopy, &pivotCols, &numPivots);
  if (numPivots < RM->nC) {
    printf("Matrix is not invertible\n");
    exit(1);
  }
  RatMat* ans = (RatMat*)malloc(sizeof(RatMat));
  RatMat_init(ans, RM->nR, RM->nC);
  for (i=0; i<RM->nR; i++) {
    for (j=0; j<RM->nC; j++) {
      mpq_set( ans->mat[i*(ans->nC)+j], RMCopy->mat[i*(RMCopy->nC)+(j+RM->nC)]);
    }
  }
  RatMat_free(RMCopy); free(RMCopy);
  RatMat_free(Id); free(Id);
  return ans;
}


void RatMat_append_ratmat(RatMat* M1, RatMat* M2) {
  int oldLen = (M1->nR)*(M1->nC);
  int M2Len = (M2->nR)*(M2->nC);
  RatMat_change_num_rows(M1, (M1->nR)+(M2->nR));
  int i;
  for (i=0; i<M2Len; i++) {
    mpq_set(M1->mat[i+oldLen], M2->mat[i]);
  }
}

void RatMat_append_cols_ratmat(RatMat* M1, RatMat* M2) {
  int i,j;
  int oldLen = (M1->nR)*(M1->nC);
  int oldNumCols = M1->nC;
  int M2Len = (M2->nR)*(M2->nC);
  mpq_t* temp = (mpq_t*)malloc(oldLen*sizeof(mpq_t));
  for (i=0; i<oldLen; i++) {
    mpq_init(temp[i]);
    mpq_set( temp[i], M1->mat[i] );
  }
  M1->mat = (mpq_t*)realloc((void*)(M1->mat), oldLen*M2Len*sizeof(mpq_t));
  M1->nC = oldNumCols + M2->nC;
  for (i=0; i<M1->nR; i++) {
    for (j=0; j<oldNumCols; j++) {
      mpq_set( M1->mat[i*(M1->nC)+j], temp[i*oldNumCols + j] );
    }
    for (j=0; j<M2->nC; j++) {
      mpq_set( M1->mat[i*(M1->nC) + oldNumCols + j], M2->mat[i*(M2->nC) + j] );
    }
  }
  for (i=0; i<oldLen; i++) {
    mpq_clear(temp[i]);
  }
  free(temp);      
}


RatMat* RatMat_identity(int dim) {
  int i,j;
  RatMat* M = (RatMat*)malloc(sizeof(RatMat));
  RatMat_init(M, dim, dim);
  for (i=0; i<dim; i++) {
    for (j=0; j<dim; j++) {
      mpq_set_si( M->mat[i*(M->nC)+j], (i==j ? 1 : 0), 1 );
    }
  }
  return M;
}
      


RatMat* RatMat_transpose(RatMat* RM) {
  int i,j;
  RatMat* ans = (RatMat*)malloc(sizeof(RatMat));
  RatMat_init(ans, RM->nC, RM->nR);
  for (i=0; i<RM->nR; i++) {
    for (j=0; j<RM->nC; j++) {
      mpq_set(ans->mat[j*(ans->nC)+i], RM->mat[i*(RM->nC)+j]);
    }
  }
  return ans;
}

RatMat* RatMat_multiply(RatMat* M1, RatMat* M2) {
  if (M1->nC != M2->nR) {
    printf("Matrices don't have compatible sizes\n");
    exit(1);
  }  
  int i,j,k;
  mpq_t sum;
  mpq_t temp;
  RatMat* ans = (RatMat*)malloc(sizeof(RatMat));
  RatMat_init(ans, M1->nR, M2->nC);
  mpq_init(sum);
  mpq_init(temp);
  for (i=0; i<M1->nR; i++) {
    for (j=0; j<M2->nC; j++) {
      mpq_set_si(sum, 0, 1);
      for (k=0; k<M1->nC; k++) {
        mpq_mul(temp, M1->mat[i*(M1->nC)+k], M2->mat[k*(M2->nC)+j]);
        mpq_add(sum, sum, temp);
      }
      mpq_set( ans->mat[i*(ans->nC)+j], sum);
    }
  }
  mpq_clear(sum);
  mpq_clear(temp);
  return ans;
}




void RatMat_change_num_rows(RatMat* M, int nr) {
  int i,j;
  if (nr == M->nR) {
    return ;
  }
  if (nr < M->nR) {
    for (i=nr; i<M->nR; i++) {
      for(j=0; j<M->nC; j++) {
        mpq_clear(M->mat[i*(M->nC)+j]);
      }
    }
  }
  M->mat = (mpq_t*)realloc((void*)(M->mat), (nr)*(M->nC)*sizeof(mpq_t)); 
  if (nr > M->nR) {
    for (i=M->nR; i<nr; i++) {
      for(j=0; j<M->nC; j++) {
        mpq_init(M->mat[i*(M->nC)+j]);
      }
    }
  }
  M->nR = nr;  
}
 
void RatMat_copy_row(RatMat* M, int r, RatMat* source, int sourceRow) {
  int i;
  for (i=0; i<M->nC; i++) {
    mpq_set( M->mat[r*(M->nC) + i], source->mat[sourceRow*(source->nC) + i] );
  }
}



RatMat* RatMat_deepcopy(RatMat* M) {
  int i;
  RatMat* answer = (RatMat*)malloc(sizeof(RatMat));
  RatMat_init(answer, M->nR, M->nC);
  for (i=0; i<M->nR; i++) {
    RatMat_copy_row(answer, i, M, i);
  }
  return answer;
}

void RatMat_print(RatMat* M, int doSmall) {
  int i,j;
  for (i=0; i<M->nR; i++) {
    printf("[");
    for (j=0; j<M->nC; j++) {
      mpq_out_str(NULL, 10, M->mat[i*(M->nC)+j]);
      if (j<M->nC-1)
        printf(",");
      if (doSmall!=1)
        printf(" ");
    }
    printf("]\n");
  }
}

void RatMat_mult_float_column(double* ans, RatMat* RM, double* func) {
  int i, j;
  double sum;
  for (i=0; i<RM->nR; i++) {
    sum = 0;
    for (j=0; j<RM->nC; j++) {
      sum += func[j] * mpq_get_d(RM->mat[i*(RM->nC)+j]);
    }
    ans[i] = sum;
  }
}
   
   
   
   
   
    
int RatMat_floatRank(RatMat* M) {
  int i,j;
  int currentRowLowerLimit = 0;
  int currentColumnLowerLimit = 0;
  int nonZeroEntryInRow;
  int nonZeroEntryInCol;
  double temp;
  int dim = M->nC;
  double temp2;
  int numPivots;
  
  double* matCopy = (double*)malloc((M->nR * M->nC)*sizeof(double));
  for (i=0; i<(M->nR)*(M->nC); i++) {
    matCopy[i] = mpq_get_d(M->mat[i]);
  }
  
  numPivots = 0;
  
  while (currentRowLowerLimit < M->nR  && currentColumnLowerLimit < dim) {  
    //find a column in which there is a one
    nonZeroEntryInRow = -1;
    for (i=currentColumnLowerLimit; i<dim; i++) {
      for (j=currentRowLowerLimit; j<M->nR; j++) {
        if ( matCopy[j*dim + i] > 0.00001 || matCopy[j*dim + i] < -0.00001 ) { //we found a nonzero entry
          nonZeroEntryInRow = j;
          nonZeroEntryInCol = i;
          numPivots++;
          goto breakAll;
        }
      }
    }
    breakAll:
    if (nonZeroEntryInRow == -1) {  //we are done row reducing
      break;
    }
    //printf("current pivot position (%d, %d):\n", nonZeroEntryInRow, nonZeroEntryInCol);
    //RatMat_print(M, 0);
    
    //divide that row
    temp = matCopy[nonZeroEntryInRow*dim + nonZeroEntryInCol];
    for (i=nonZeroEntryInCol; i<dim; i++) {   
      matCopy[nonZeroEntryInRow*dim + i] = matCopy[nonZeroEntryInRow*dim + i] / temp;
    }
      
    //swap that row
    for (i=0; i<M->nC; i++) {
      temp = matCopy[currentRowLowerLimit*(M->nC) + i];
      matCopy[currentRowLowerLimit*(M->nC) + i] = matCopy[nonZeroEntryInRow*(M->nC) + i];
      matCopy[nonZeroEntryInRow*(M->nC) + i] = temp;
    }
    
    //printf("Cleared the row and swapped\n");
    //RatMat_print(M, 0);  
        
    //clear everything in that column
    for (i=0; i<M->nR; i++) {
      if (i==currentRowLowerLimit) {
        continue;
      }
      if (matCopy[i*dim+nonZeroEntryInCol] < 0.000001 &&  matCopy[i*dim+nonZeroEntryInCol] > -0.000001) {  //it's 0
        continue;
      }
      temp = -matCopy[i*dim+nonZeroEntryInCol];

      //printf("clearing row %d (the factor is ", i);
      for (j=nonZeroEntryInCol; j<dim; j++) {
        temp2 = matCopy[currentRowLowerLimit*dim+j];
        temp2 *= temp;
        matCopy[i*dim+j] = matCopy[i*dim+j] + temp2;
      }
    }
    //printf("Cleared the column\n");
    //RatMat_print(M, 0);
    
    //set the new lower limits
    currentRowLowerLimit++;
    currentColumnLowerLimit = nonZeroEntryInCol + 1;
  }
  
  //for (i=0; i<M->nR; i++) {
  //  for (j=0; j<M->nC; j++) {
  //    printf("%f ", matCopy[i*M->nC + j]);
  //  }
  //  printf("\n");
  //}
  
  
  
  free(matCopy);
  return numPivots;
}



double* RatMat_get_soln_to_aug_float(RatMat* M) {
  int i,j;
  int currentRowLowerLimit = 0;
  int currentColumnLowerLimit = 0;
  int nonZeroEntryInRow;
  int nonZeroEntryInCol;
  double temp;
  int dim = M->nC;
  double temp2;
  int numPivots;
  
  double* matCopy = (double*)malloc((M->nR * M->nC)*sizeof(double));
  for (i=0; i<(M->nR)*(M->nC); i++) {
    matCopy[i] = mpq_get_d(M->mat[i]);
  }
  
  numPivots = 0;
  
  while (currentRowLowerLimit < M->nR  && currentColumnLowerLimit < dim) {  
    //find a column in which there is a one
    nonZeroEntryInRow = -1;
    for (i=currentColumnLowerLimit; i<dim; i++) {
      for (j=currentRowLowerLimit; j<M->nR; j++) {
        if ( matCopy[j*dim + i] > 0.00001 || matCopy[j*dim + i] < -0.00001 ) { //we found a nonzero entry
          nonZeroEntryInRow = j;
          nonZeroEntryInCol = i;
          numPivots++;
          goto breakAll;
        }
      }
    }
    breakAll:
    if (nonZeroEntryInRow == -1) {  //we are done row reducing
      break;
    }
    //printf("current pivot position (%d, %d):\n", nonZeroEntryInRow, nonZeroEntryInCol);
    //RatMat_print(M, 0);
    
    //divide that row
    temp = matCopy[nonZeroEntryInRow*dim + nonZeroEntryInCol];
    for (i=nonZeroEntryInCol; i<dim; i++) {   
      matCopy[nonZeroEntryInRow*dim + i] = matCopy[nonZeroEntryInRow*dim + i] / temp;
    }
      
    //swap that row
    for (i=0; i<M->nC; i++) {
      temp = matCopy[currentRowLowerLimit*(M->nC) + i];
      matCopy[currentRowLowerLimit*(M->nC) + i] = matCopy[nonZeroEntryInRow*(M->nC) + i];
      matCopy[nonZeroEntryInRow*(M->nC) + i] = temp;
    }
    
    //printf("Cleared the row and swapped\n");
    //RatMat_print(M, 0);  
        
    //clear everything in that column
    for (i=0; i<M->nR; i++) {
      if (i==currentRowLowerLimit) {
        continue;
      }
      if (matCopy[i*dim+nonZeroEntryInCol] < 0.000001 &&  matCopy[i*dim+nonZeroEntryInCol] > -0.000001) {  //it's 0
        continue;
      }
      temp = -matCopy[i*dim+nonZeroEntryInCol];

      //printf("clearing row %d (the factor is ", i);
      for (j=nonZeroEntryInCol; j<dim; j++) {
        temp2 = matCopy[currentRowLowerLimit*dim+j];
        temp2 *= temp;
        matCopy[i*dim+j] = matCopy[i*dim+j] + temp2;
      }
    }
    //printf("Cleared the column\n");
    //RatMat_print(M, 0);
    
    //set the new lower limits
    currentRowLowerLimit++;
    currentColumnLowerLimit = nonZeroEntryInCol + 1;
  }
  
  //for (i=0; i<M->nR; i++) {
  //  for (j=0; j<M->nC; j++) {
  //    printf("%f ", matCopy[i*M->nC + j]);
  //  }
  //  printf("\n");
  //}
  
  if (numPivots < M->nC-1) {
    free(matCopy);
    return NULL;
  } else {
    for (i=M->nC; i<M->nR; i++) {
      if (matCopy[i*(M->nC)+(M->nC-1)] > 0.00001 || matCopy[i*(M->nC)+(M->nC-1)] < -0.00001) {
        free(matCopy);
        return NULL;
      }
    }
    double* soln = (double*)malloc((M->nC-1)*sizeof(double));
    for (i=0; i< M->nC-1; i++) {
      //if (pivotCols[i] != i) {
      //  printf("I'm confused in getting solution\n");
      //  exit(1);
      //}
      soln[i] = matCopy[i*(M->nC)+(M->nC-1)]; // the last entry
    }
    free(matCopy);
    return soln;
  }
  

}



