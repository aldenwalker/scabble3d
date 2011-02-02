/*
 *  ball_worker.h
 *  
 *
 */

#ifndef __ball_worker__
#define __ball_worker__

#include <gmp.h>
#include <gtk/gtk.h>
#include <semaphore.h>

#include "matrix.h"
#include "triangle_and_vertex.h"


enum scallop_lp_solver {GLPK_DOUBLE, GLPK_EXACT, QSOPT_EXACT, EXLP};


/*****************************************************************************/
/* these are the arc and polygon types                                       */
/*****************************************************************************/
typedef struct {
  int first_word;
  int last_word;
  int first;
  int last;
} arc;

typedef struct {
  int* arc;
  int num_arcs;
} polygon;


/*****************************************************************************/
/* this describes an scl problem (the chains, words, weights, etc)           */
/*****************************************************************************/
typedef struct {
  char*** chains;
  int num_chains;
  int* chain_lens;
  char** word_list;
  int num_words;
  int* weights;    //this is a list through all the words
  RatMat* constraints; //these are the basic constraints
  int* equality_type;  //
  arc* arc_list;
  int num_arcs;
  polygon* poly_list;
  int num_polys;
} scl_problem;

/*****************************************************************************/
/* this is a particular orthant calculation                                  */
/*****************************************************************************/
typedef struct {
  int orthant_num;
  scl_problem* scl_prob;
  tri_list* triangles;
  vert_list* vertices;
  double max_undone_triangle_area;
  int is_complete;
} orthant_problem;

/*****************************************************************************/
/* this type holds the information about a particular unit ball calculation  */
/*****************************************************************************/
typedef struct {
  char*** chains; //duplicated inside the orthant_problems, but that's ok
                  //these are the chains for orthant 0
  int num_chains;
  int* chain_lens;
  orthant_problem** orthants;  //there are always 4 of these
  int is_complete;
  int current_working_orthant;
  double tolerance;
} ball_problem;


/*****************************************************************************/
/* this is the information about a currently running (or paused) execution   */
/*****************************************************************************/
typedef struct {
  ball_problem* ball;       //this is the problem we are working on
  enum scallop_lp_solver solver;
  
  sem_t running_sem;         //this is so the worker can block if necessary;
  sem_t message_sem;         //this blocks for message passing
  
  //message stuff (should be blocked by message_sem)
  int one_step;              //note really a message -- must be set at the beginning -- only do one step
  int status;                //1 means currently running, 0 means not
  int status_message;        //1 means stop, 0 means continue
  int new_tolerance_check;  //1 means there is a new tolerance 0 means nope
  double new_tolerance;
  int skip_orthant;        //1 means yes, skip, 0 means no keep going
  
  sem_t read_data_sem;     //this tells the gui thread to read the data
                           //note that the worker thread is waiting for this
  char* initial_arguments[3];                         
   
  GtkWidget* target_drawing_area; //this is the widget to draw to
  
} execution;


void init_lp();
void linear_program_from_ratmat(polygon* poly_list,
                                rvector* solution_vector,
                                mpq_t scl,
                                RatMat* constraints,
                                int* equality_type,
                                enum scallop_lp_solver solver);


void generate_arcs(arc** arc_list, 
                   int* num_arcs,
                   char** word_list,
                   int word_list_len);
void generate_polygons(char** word_list,
                       int word_list_len,
                       arc* arc_list,
                       int num_arcs,
                       polygon** poly_list,
                       int* num_polys,
                       int maxjun);
void create_constraint_matrix(char*** chains, 
                              int num_chains,
                              int* chain_lens,
                              arc* arc_list,
                              int num_arcs,
                              polygon* polygon_list,
                              int num_polys,
                              int* weights,
                              RatMat* constraints,
                              int** equalityType);
void scl_problem_init(scl_problem* scl_prob, 
                      char*** chains,
                      int num_chains,
                      int* chain_lens,
                      char** word_list,
                      int num_words,
                      int* weights,
                      int maxjun);



void point_scl(scl_problem* scl_prob,
               rvector* point,
               mpq_t scl,
               enum scallop_lp_solver solver);

int min_scl_over_triangle(scl_problem* scl_prob, 
                          vert_list* V,
                          triangle* t,
                          mpq_t scl,
                          rvector* new_vertex,
                          enum scallop_lp_solver solver);


int find_undone_triangle(tri_list* T, 
                         double tolerance);

void split_triangles(vert_list* V, tri_list* T, int split_ind, int new_vert);

int find_edge_for_vertex(vert_list* V, triangle* t, int new_vert);

int one_orthant_step(orthant_problem* orth, double tolerance, enum scallop_lp_solver solver);

void one_computation_step(ball_problem* ball, enum scallop_lp_solver solver);

void* run_execution(void* E_void);

void orthant_problem_init(orthant_problem* orth, 
                          int index, 
                          char*** chains, 
                          int* chain_lens,
                          int* weights,
                          int num_words,
                          mpq_t* predone_scls,
                          int maxjun,
                          enum scallop_lp_solver solver);
                          
void computation_init(execution* E,
                      char*** chains,
                      int* chain_lens,
                      int* weights,
                      int num_words,
                      double tolerance,
                      int maxjun,
                      enum scallop_lp_solver solver,
                      GtkWidget* target_drawing_area);


void scl_problem_print(scl_problem* sp);
void orthant_problem_print(orthant_problem* orth);
void ball_problem_print(ball_problem* ball);
void execution_print(execution* E);




#endif




















