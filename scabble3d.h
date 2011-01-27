

/*****************************************************************************/
/* this is a (3d probably) rational vector                                   */
/*****************************************************************************/
typedef struct {
  mpq_t* coord;
  int dim;
  int malloced;
} rvector;

/*****************************************************************************/
/* this is a 3d triangle                                                     */
/*****************************************************************************/
typedef struct {
  int* verts;    //this is the vertex indices in the vertex list
                 //the number of vertices should be implicit
  double area;
  int is_scl_linear; //1 means yes, 0 means unknown, -1 means no
} triangle;
  
/*****************************************************************************/
/* this is a list  of triangles                                              */
/*****************************************************************************/
typedef struct {
  triangle* tris;
  int num_tris;
  int first_nonlinear_triangle;   //this records the index of the first triangle
                                  //which isn't linear (which we can work on)
} tri_list;

/*****************************************************************************/
/* this is a list  of vertices                                               */
/*****************************************************************************/
typedef struct {
  rvector* verts;
  int num_verts;
} vert_list;


/*****************************************************************************/
/* this describes an scl problem (the chains, words, weights, etc)           */
/*****************************************************************************/
typedef struct {
  char*** chains;
  int num_chains;
  int* chain_lens;
  char** all_words;
  int num_words;
  int* weights;    //this is a list through all the words
  RatMat* constraints; //these are the basic constraints
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
  
  sem_t status_sem;
  int status;                //1 means currently running, 0 means not
  
  sem_t running_sem;         //this is so the worker can block if necessary;
  
  sem_t status_message_sem;
  int status_message;        //1 means stop, 0 means continue
  
  sem_t new_tolerance_sem;
  int new_tolerance_check;  //1 means there is a new tolerance 0 means nope
  double new_tolerance;
  
  sem_t skip_orthant_sem;
  int skip_orthant;        //1 means yes, skip, 0 means no keep going
  
  sem_t read_data_sem;     //this tells the gui thread to read the data
                           //note that the worker thread is waiting for this
} execution;
  
  
  









