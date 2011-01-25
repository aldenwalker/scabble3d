/*****************************************************************************/
/*  A matrix class                                                           */
/*****************************************************************************/
#ifndef _mat_defined
#define _mat_defined

//#include "gmp/gmp-exec/include/gmp.h"
#include <gmp.h>
//#include "sssgmp.h"

typedef struct {
  int* mat;
  int nR;
  int nC;
} Mat;

typedef struct {
  mpq_t* mat;
  int nR;
  int nC;
} RatMat;

void Mat_init(Mat* M, int nR, int nC);
void Mat_set(Mat* M, int r, int c, int val);
Mat* Mat_deepcopy(Mat* source);
void Mat_reset(Mat* M, int val);
inline int Mat_get(Mat* M, int r, int c);
void Mat_increment(Mat* M, int r, int c, int incr);
void Mat_print(Mat* M, int doSmall);
void Mat_free(Mat* M);
void Mat_print_row(Mat* M, int row, int doSmall);
void Mat_copy_row(Mat* M, int r, Mat* source, int sourceRow);
void Mat_copy_neg_row(Mat* M, int r, Mat* source, int sourceRow);
void Mat_copy_list_row(Mat* M, int r, int* L);
void Mat_change_num_rows(Mat* M, int nr);
Mat* Mat_matrix_from_rows(Mat* source, int* rows, int len);
void Mat_remove_duplicate_rows(Mat* M);
void Mat_sort_rows(Mat* M);
int Mat_rows_are_equal(Mat* M, int i, int j);
int Mat_dot_prod_with_row(Mat* M, int row, int* L);
double Mat_dot_prod_with_float_row(Mat* M, int row, double* L);
Mat* Mat_nullspace(Mat* M); 
Mat* Mat_multiply(Mat* M1, Mat* M2);
Mat* Mat_transpose(Mat* M);
int Mat_rank(Mat* M);
void Mat_append_mat(Mat* M1, Mat* M2);
RatMat* Mat_get_soln_to_aug(Mat* M);



RatMat* Mat_get_ratmat(Mat* M);

void RatMat_init(RatMat* M, int nR, int nC);
void RatMat_set(RatMat* M, int r, int c, mpq_t val);
void RatMat_set_int(RatMat* M, int r, int c, int n, int d);
void RatMat_reset_int(RatMat* M, int val);
void RatMat_get(RatMat* M, int r, int c, mpq_t a);
void RatMat_free(RatMat* M);
void RatMat_rref(RatMat* M, int** pivotCols, int* numPivots);
void RatMat_swap_rows(RatMat* M, int i, int j);  
RatMat* RatMat_inverse(RatMat* RM);
RatMat* RatMat_transpose(RatMat* RM);
RatMat* RatMat_multiply(RatMat* M1, RatMat* M2);
void RatMat_append_ratmat(RatMat* M1, RatMat* M2);
void RatMat_append_cols_ratmat(RatMat* M1, RatMat* M2);
RatMat* RatMat_identity(int dim);
void RatMat_change_num_rows(RatMat* M, int nr);
void RatMat_mult_float_column(double* ans, RatMat* RM, double* func);  
void RatMat_copy_row(RatMat* M, int r, RatMat* source, int sourceRow);
RatMat* RatMat_deepcopy(RatMat* M);
void RatMat_print(RatMat* M, int doSmall);
int RatMat_rank(RatMat* M);
RatMat* RatMat_get_soln_to_aug(RatMat* M);
int RatMat_floatRank(RatMat* M);
double* RatMat_get_soln_to_aug_float(RatMat* M);


#endif

