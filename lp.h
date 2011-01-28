
/*
 *  lp.h
 *  
 *
 *
 */

#ifndef __lp__
#define __lp__

#include <gmp.h>
#include "scabble3d.h"
#include "matrix.h"



enum scallop_lp_solver {GLPK_DOUBLE, GLPK_EXACT, QSOPT_EXACT, EXLP};



void linear_program_from_ratmat(polygon* poly_list,
                                rvector* solution_vector,
                                mpq_t scl,
                                RatMat* constraints,
                                int* equality_type,
                                enum scallop_lp_solver solver);

#endif