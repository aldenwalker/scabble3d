/*
 *  lp.h
 *  
 *
 *
 */

#include <gmp.h>
#include "scabble3d.h"



enum scallop_lp_solver {GLPK_DOUBLE, GLPK_EXACT, QSOPT_EXACT, EXLP};



void linear_program_from_ratmat(polygon* poly_list,
                                rvector* solution_vector,
                                mpq_t scl,
                                RatMat* constraints,
                                int* equality_type,
                                scallop_lp_solver solver)