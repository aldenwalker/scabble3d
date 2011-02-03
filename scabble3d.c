#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <semaphore.h>
#include <pthread.h>

#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <GL/gl.h>
#include <gmp.h>

#include "scabble3d.h"
#include "triangle_and_vertex.h"
#include "matrix.h"
#include "scl_backend.h"



/*****************************************************************************/
/* the data structure for holding all the widget info                        */
/*****************************************************************************/
typedef struct {
  GtkEntry* entry1;
  GtkEntry* entry2;
  GtkEntry* entry3;
  GtkEntry* tol_entry;
  GtkEntry* file;
  GtkWidget* m5_check;
  GtkWidget* drawing;
} fieldList;

/*****************************************************************************/
/* this gives our current status in terms of rotation, etc                   */
/*****************************************************************************/
struct {
  double button1_down_x;
  double button1_down_y;
  double button3_down_x;
  double button3_down_y;
  double current_rotation_x;
  double current_rotation_y;
  double constant_rotation_x;
  double constant_rotation_y;
} dstatus;

//these are the opengl globals that I think I need
GdkPixmap* pixmap = NULL;
GdkGLPixmap* GLPixmap = NULL;
GdkGLContext* GLContext = NULL;
ball_mesh* GLmesh = NULL;

//this is the global execution going on 
execution* EGlobal = NULL;

//I've got so many globals, I might as well throw in verbose too, for
//simplicity
int VERBOSE = 0;


/*****************************************************************************/
/* this function draws the opengl scene                                      */
/* note that the GLmesh holds some triangles and normals which must be       */
/* flipped in orientation                                                    */
/*****************************************************************************/
void draw_mesh() {
  int i,j,k,orth;
  GdkGLDrawable* gldrawable = gdk_pixmap_get_gl_drawable(pixmap);
  double normal[3];
  double triangle[3][3];
  double temp_d;
  
  
  if (!gdk_gl_drawable_gl_begin(gldrawable, GLContext)) {
    printf("Couldn't start opengl drawing\n");
  }
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  
  GLfloat blueish[] = {0.2, 0.2, 1, 1};
  GLfloat redish[] = {1, 0.2, 0.2, 1};
  GLfloat red[] = {1,0,0,1};
  GLfloat green[] = {0,1,0,1};
  GLfloat blue[] = {0,0,1,1};
  GLfloat spec[] = {0.1, 0.1, 0.1, 1};
  
  glRotatef(dstatus.current_rotation_x + dstatus.constant_rotation_x, 0, 1, 0);
  glRotatef(dstatus.current_rotation_y + dstatus.constant_rotation_y, -1, 0, 0);
  
  //draw the axis
  glMaterialfv(GL_FRONT, GL_DIFFUSE, red);
  glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(4,0,0);
  glEnd();
  glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
  glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(0,4,0);
  glEnd();
  glMaterialfv(GL_FRONT, GL_DIFFUSE, blue);
  glBegin(GL_LINES);
    glVertex3f(0,0,0);
    glVertex3f(0,0,4);
  glEnd();
  
  
  //go through the triangles and print them
  //I want to make this more efficient
  glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
  for (orth=0; orth<4; orth++) {
    //printf("For orthant %d, I'm displaying the following triangles:\n", orth);
    //tri_list_print_d(GLmesh->triangles[orth], GLmesh->vertices[orth]);
    for (i=0; i<GLmesh->triangles[orth]->num_tris; i++) {
      if (GLmesh->triangles[orth]->tris[i].is_scl_linear == 1) {
        glMaterialfv(GL_FRONT, GL_DIFFUSE, blueish);
      } else {
        glMaterialfv(GL_FRONT, GL_DIFFUSE, redish);
      }
      
      //load the normal and triangle
      for (j=0; j<3; j++) {
        normal[j] = GLmesh->normals[orth]->verts[i][j];
        for (k=0; k<3; k++) {
          triangle[j][k] =  GLmesh->vertices[orth]->verts[ GLmesh->triangles[orth]->tris[i].verts[j] ][k];
        }
      }
      //if we're doing orthants 1 or 3, then we need to swap the orientation of
      //the triangle and reverse the normal
      if (orth==1 || orth==2) {
        for (j=0; j<3; j++) {
          normal[j] = -normal[j];
          temp_d = triangle[1][j];
          triangle[1][j] = triangle[2][j];
          triangle[2][j] = temp_d;
        }
      }
      
      glBegin(GL_TRIANGLES);
      //FRONT
      //printf("Normal: [%f,%f,%f]\n", normal[0], 
      //                               normal[1], 
      //                               normal[2]);
      glNormal3f(normal[0], 
                 normal[1], 
                 normal[2]);
      for (j=0; j<3; j++) {
        //printf("\t\tVertex: [%f,%f,%f]\n", triangle[j][0],
        //                                   triangle[j][1],
        //                                   triangle[j][2] );
        glVertex3f( triangle[j][0],
                    triangle[j][1],
                    triangle[j][2] );
      }
      //BACK
      //printf("Normal: [%f,%f,%f]\n", -GLmesh->normals[orth]->verts[i][0], 
      //                               -GLmesh->normals[orth]->verts[i][1], 
      //                               -GLmesh->normals[orth]->verts[i][2]);      
      glNormal3f(-normal[0], 
                 -normal[1], 
                 -normal[2]);
      for (j=0; j<3; j++) {
        //printf("\t\tVertex: [%f,%f,%f]\n", -GLmesh->vertices[orth]->verts[ GLmesh->triangles[orth]->tris[i].verts[j] ][0],
        //                                -GLmesh->vertices[orth]->verts[ GLmesh->triangles[orth]->tris[i].verts[j] ][1],
        //                                -GLmesh->vertices[orth]->verts[ GLmesh->triangles[orth]->tris[i].verts[j] ][2] );
        glVertex3f( -triangle[j][0],
                    -triangle[j][1],
                    -triangle[j][2] );
      }
      glEnd();
          
    }
  }
 
  glPopMatrix();
  
  if(gdk_gl_drawable_is_double_buffered(gldrawable)) {
    gdk_gl_drawable_swap_buffers(gldrawable);
  } else {
	  glFlush();
	}

	gdk_gl_drawable_gl_end(gldrawable);
  
  gdk_window_invalidate_rect(EGlobal->target_drawing_area->window, 
                             &EGlobal->target_drawing_area->allocation, 
                             FALSE);
	gdk_window_process_updates(EGlobal->target_drawing_area->window, FALSE);
}


/*****************************************************************************/
/* this function draws the polygon ball                                      */
/*****************************************************************************/
void update_ball_picture_while_running(execution* E) {
  int i,j,k,l;
  double new_vert[3];
  double point[3][3];
  double diff1[3];
  double diff2[3];
  double normal[3];
  
  if (VERBOSE) 
    printf("Hey I'm drawing the ball now\n");
  
   
  tri_list* T = NULL;
  vert_list* V = NULL;
  for (i=0; i<4; i++) {
    T = E->ball->orthants[i]->triangles;
    V = E->ball->orthants[i]->vertices;
    
    //add the new vertices
    for (j=GLmesh->vertices[i]->num_verts; j<V->num_verts; j++) {
      for (k=0; k<3; k++) {
        new_vert[k] = mpq_get_d(V->verts[j].coord[k]);
      }
      new_vert[0] = ((i&1)==1 ? -new_vert[0] : new_vert[0]);
      new_vert[1] = (((i>>1)&1)==1 ? -new_vert[1] : new_vert[1]);
      vert_list_d_add_copy(GLmesh->vertices[i], new_vert);
    }
    
    //delete triangles in the mesh which have disappeared
    j = 0;
    while (j<GLmesh->triangles[i]->num_tris) {
      for (k=0; k<T->num_tris; k++) {
        if (T->tris[k].verts[0] == GLmesh->triangles[i]->tris[j].verts[0]
           && T->tris[k].verts[1] == GLmesh->triangles[i]->tris[j].verts[1]
           && T->tris[k].verts[2] == GLmesh->triangles[i]->tris[j].verts[2]){
          break;
        }
      }
      if (k<T->num_tris) {
        j++;
      } else {
        tri_list_delete_index(GLmesh->triangles[i], j);
        vert_list_d_delete_index(GLmesh->normals[i], j);
      }
    }
    
    //add triangles in orthant which are new
    for (j=0; j<T->num_tris; j++) {
      for (k=0; k<GLmesh->triangles[i]->num_tris; k++) { 
        if (T->tris[j].verts[0] == GLmesh->triangles[i]->tris[k].verts[0]
           && T->tris[j].verts[1] == GLmesh->triangles[i]->tris[k].verts[1]
           && T->tris[j].verts[2] == GLmesh->triangles[i]->tris[k].verts[2]){
          break;
        }
      }
      if (k <  GLmesh->triangles[i]->num_tris) {
        //ok the triangle is in there, but we just need to make sure that we
        //keep is_scl_linear updated
        GLmesh->triangles[i]->tris[k].is_scl_linear = T->tris[j].is_scl_linear;
        continue;
      } else {
        tri_list_add_copy(GLmesh->triangles[i], &T->tris[j]);
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {
            point[k][l] = mpq_get_d(V->verts[T->tris[j].verts[k]].coord[l]); 
          }
          //we don't need to 
          point[k][0] = ((i&1)==1 ? -point[k][0] : point[k][0]);
          point[k][1] = (((i>>1)&1)==1 ? -point[k][1] : point[k][1]);
        }
        for (k=0; k<3; k++) {
          diff1[k] = point[1][k] - point[0][k];
          diff2[k] = point[2][k] - point[0][k];
        }
        dvector_cross(normal, diff1, diff2);
        //printf("I took the cross of [%f,%f,%f] and [%f,%f,%f] and got [%f,%f,%f]\n",
        //        diff1[0], diff1[1], diff1[2], 
        //        diff2[0], diff2[1], diff2[2],
        //        normal[0], normal[1], normal[2]);
        //if we only swapped one coordinate, we need to 
        //reverse orientation on the normals
        //if (i==1 || i==3) {
        //  normal[0] = -normal[0]; normal[1] = -normal[1]; //normal[2] = -normal[2];
        //}
        vert_list_d_add_copy(GLmesh->normals[i], normal);
      }
    }
  }

  draw_mesh();
  
  sem_post(&(E->read_data_sem));
  

}


/*****************************************************************************/
/* this function takes the basic input from the text boxes and starts it up  */
/*****************************************************************************/
void load_inputs_and_run(char* arg1, 
                         char* arg2, 
                         char* arg3, 
                         double tolerance,
                         int maxjun,
                         enum scallop_lp_solver solver,
                         GtkWidget* target_drawing_area) {
  char*** chains = (char***)malloc(3*sizeof(char**));
  int i;
  int j;
  int current_word;
  int current_total_word;
  int prev_j;
  int* chain_lens = (int*)malloc(3*sizeof(int));
  int* weights = NULL;
  int num_words = 0;
  char* args[] = { arg1, arg2, arg3 };
  for (i=0; i<3; i++) {
    chain_lens[i] = 0;
    prev_j = 0;
    j=1;
    while (1) {
      if (args[i][j] == ' ' || args[i][j] == '\0') {
        num_words++;
        chain_lens[i] ++;
        if (args[i][j] == '\0') {
          break;
        }
      }
      j++;
    }
  }
  
  if (VERBOSE) {
    printf("Called with: %s, %s, %s,\n", arg1, arg2, arg3);
    printf(" = %s, %s, %s, \n", args[0], args[1], args[2]);
    printf("I found chain lens: %d, %d, %d\n", chain_lens[0], chain_lens[1], chain_lens[2]);
    printf("For a total of %d words\n", num_words);
  }
  
  weights = (int*)malloc(num_words*sizeof(int));
  current_total_word = 0;
  
  //now actually put the chains in 
  for (i=0; i<3; i++) {
    current_word = 0;
    chains[i]  = (char**)malloc(chain_lens[i]*sizeof(char*));
    prev_j = 0;
    j = 1;
    while (1) {
      //printf("prev_j = %d, j = %d\n", prev_j, j);
      if (args[i][j] == ' ' || args[i][j] == '\0') {
        //printf("Found a breakpoint at %d -- %d\n", prev_j, j);
        chains[i][current_word] = (char*)malloc((j-prev_j+1)*sizeof(char));
        strncpy(chains[i][current_word], &(args[i][prev_j]), (j-prev_j));
        chains[i][current_word][j-prev_j] = '\0';
        weights[current_total_word] = 1;
        current_total_word++;
        current_word++;
        if (args[i][j] == ' ') {
          prev_j = j+1;
          j += 2;
        } else {
          break;
        }
      } else {
        j++;
      }
    }
  }
  
  if (VERBOSE) {
    printf("I got the following chains:\n");
    for (i=0; i<3; i++) {
      for (j=0; j<chain_lens[i]; j++) {
        printf("%s ", chains[i][j]);
      }
      printf("\n");
    }
  }
  
  execution* E = (execution*)malloc(sizeof(execution));
  
  //build the computation -- make initial scl computations, etc
  computation_init(E, chains, 
                      chain_lens, 
                      weights, 
                      num_words, 
                      tolerance, 
                      maxjun, 
                      solver,
                      target_drawing_area,
                      VERBOSE);
  
  if (EGlobal != NULL) {
    execution_free(EGlobal);
    free(EGlobal);
  }
  
  EGlobal = E;
  
  for (i=0; i<3; i++) {
    EGlobal->initial_arguments[i] = (char*)malloc((strlen(args[i])+1)*sizeof(char));
    strcpy(EGlobal->initial_arguments[i], args[i]);
  }
  
  //execution_print(E);
  
  //initialize the GLmesh
  GLmesh = (ball_mesh*)malloc(sizeof(ball_mesh));
  for (i=0; i<4; i++) {
    GLmesh->triangles[i] = (tri_list*)malloc(sizeof(tri_list));
    GLmesh->vertices[i] = (vert_list_d*)malloc(sizeof(vert_list_d));
    GLmesh->normals[i] = (vert_list_d*)malloc(sizeof(vert_list_d));
    GLmesh->triangles[i]->num_tris = 0;
    GLmesh->triangles[i]->tris = NULL;
    GLmesh->vertices[i]->num_verts = 0;
    GLmesh->vertices[i]->verts = NULL;
    GLmesh->normals[i]->num_verts = 0;
    GLmesh->normals[i]->verts = NULL;
  }
    
  sem_wait(&(E->read_data_sem));
  update_ball_picture_while_running(E);
  
  if (VERBOSE) {
    printf("Drew the initial ball\n");  
    printf("Done computation init -- about to start pthread\n"); 
  }
    
  //start the multi-threadedness
  pthread_t worker_thread;
  pthread_create(&worker_thread, NULL, run_execution, (void*)E);
  
  if (VERBOSE) 
    printf("started computation thread\n");
  
}

/*****************************************************************************/
/* find the rank (count the number of distinct letters)                      */
/*****************************************************************************/
int find_rank(char* e1, char* e2, char* e3) {
  char letters[20];
  int i,j,which;
  char* es[3] = {e1, e2, e3};
  int l[3] = { strlen(e1), strlen(e2), strlen(e3) };
  letters[0] = '\0';
  for (which=0; which<3; which++) {
    for (i=0; i<l[which]; i++) {
      j=0;
      while (letters[j] != '\0') {
        if (letters[j] == es[which][i]) {
          break;
        }
        j++;
      }
      if (letters[j] == '\0') {
        letters[j] = es[which][i];
        letters[j+1] = '\0';
      }
    }
  }
  return strlen(letters)/2;
}
  
  
/*****************************************************************************/
/* event handling functions:                                                 */
/*****************************************************************************/
static gboolean run_button_press(GtkWidget* widget,
                                 GdkEventButton* event,
                                 fieldList* fields) {
  char* e1;
  char* e2;
  char* e3;
  int rank;
  
  //if we're just pressing, don't do anything
  if (event->type == GDK_BUTTON_PRESS) {
    return FALSE;
  }
  
    
  //if something is running currently, don't do anything
  if (EGlobal != NULL) {
    sem_wait(&EGlobal->message_sem);
    if (EGlobal->status == 1) {
      sem_post(&EGlobal->message_sem);
      return TRUE;
    }
    sem_post(&EGlobal->message_sem);
  }
  
  //load the fields from the entries
  e1 = (char*)gtk_entry_get_text(fields->entry1);
  e2 = (char*)gtk_entry_get_text(fields->entry2);
  e3 = (char*)gtk_entry_get_text(fields->entry3);
  
  //find the rank
  rank = find_rank(e1, e2, e3);  
  if (VERBOSE) printf("I found the rank to be %d\n", rank);
  
  //check if either: we don't have any execution, or we've changed the inputs
  //so we have to re-init
  if (EGlobal == NULL || 
       (strcmp(EGlobal->initial_arguments[0], e1)!=0
        || strcmp(EGlobal->initial_arguments[1], e2)!=0
        || strcmp(EGlobal->initial_arguments[2], e3)!=0
        || (EGlobal->maxjun == 2*rank+1 && !gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(fields->m5_check)))
        || (EGlobal->maxjun == 2*rank && gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(fields->m5_check))) ) ) {
  
     if (strlen(e1) == 0
        || strlen(e2) == 0
        || strlen(e3) == 0) {
      printf("You forgot a chain or something\n");
      return TRUE;
    }
    if (VERBOSE) {
      printf("Looks like you want to start a new computation\n"); fflush(stdout);
    }
    double tolerance = atof((char*)gtk_entry_get_text(fields->tol_entry));
    load_inputs_and_run((char*)gtk_entry_get_text(fields->entry1),
                        (char*)gtk_entry_get_text(fields->entry2),
                        (char*)gtk_entry_get_text(fields->entry3),
                        tolerance,
                        (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(fields->m5_check)) ? 2*rank+1 : 2*rank),
                        EXLP,
                        fields->drawing);
                        

  } else { 
    //we already have an execution -- run it
    //make the tolerance whatever is in the box
    EGlobal->ball->tolerance = atof((char*)gtk_entry_get_text(fields->tol_entry));
    
    if (VERBOSE) {
      printf("Loading new tolerance %f\n", EGlobal->ball->tolerance);
      printf("from field with text %s\n", (char*)gtk_entry_get_text(fields->tol_entry));
    }
    pthread_t worker_thread;
    EGlobal->one_step = 0;
    EGlobal->status_message = 0;
    pthread_create(&worker_thread, NULL, run_execution, (void*)EGlobal);
    
    if (VERBOSE) 
      printf("started computation thread\n");  
  }
 
  return FALSE;
}       

/*****************************************************************************/
/* pause a currently-running execution                                       */
/*****************************************************************************/
static gboolean pause_button_press(GtkWidget* widget,
                                   GdkEventButton* event,
                                   fieldList* fields) {
  if (event->type == GDK_BUTTON_RELEASE && EGlobal != NULL) {
    sem_wait(&EGlobal->message_sem);
    EGlobal->status_message = 1;
    sem_post(&EGlobal->message_sem);
  }
  return FALSE;
}

/*****************************************************************************/
/* do a single execution step  (one linear program)                          */
/*****************************************************************************/
static gboolean step_button_press(GtkWidget* widget,
                                  GdkEventButton* event,
                                  fieldList* fields) {
  if (event->type == GDK_BUTTON_RELEASE && EGlobal != NULL) {
    sem_wait(&EGlobal->message_sem);
    if (EGlobal->status == 1) {
      sem_post(&EGlobal->message_sem);
      return TRUE;
    }
    sem_post(&EGlobal->message_sem);
    pthread_t worker_thread;
    EGlobal->one_step = 1;
    pthread_create(&worker_thread, NULL, run_execution, (void*)EGlobal);    
  }
  return FALSE;
}
    
/*****************************************************************************/
/* save the current ball as a pov file                                       */
/*****************************************************************************/    
static gboolean save_pov_button_press(GtkWidget* widget,
                                      GdkEventButton* event,
                                      fieldList* fields) {  
  int i,orth;
  double rotx, roty;
  char* file_text = (char*)gtk_entry_get_text(fields->file);
  if (event->type == GDK_BUTTON_PRESS || EGlobal == NULL) {
    return FALSE;
  } 
  
  //if we're running, don't mess with it
  sem_wait(&EGlobal->message_sem);
  if (EGlobal->status == 1) {
    sem_post(&EGlobal->message_sem);
    return FALSE;
  }
  sem_post(&EGlobal->message_sem);
  
  FILE* out_file = fopen(file_text, "w");
  if (out_file == NULL) {
    printf("Doesn't look like I can open that file\n");
    return FALSE;
  }
  fprintf(out_file, "#include \"colors.inc\"\n");
  fprintf(out_file, "background { color White }\n");
  fprintf(out_file, "camera {\n");
  rotx = (3.1415926535/180.0)*dstatus.constant_rotation_x;
  roty = (3.1415926535/180.0)*dstatus.constant_rotation_y;
  fprintf(out_file, "\tlocation <%f, %f, %f>\n",  8*sin(-rotx),
                                                  8*sin(-roty)*cos(rotx),
                                                  8*cos(roty)*cos(rotx));
  fprintf(out_file, "\tlook_at <0, 0, 0>\n");
  fprintf(out_file, "\tright<-1.33, 0, 0>\n");
  fprintf(out_file, "}\n");
  fprintf(out_file, "light_source { <0, 0, 7> color White}\n");
  fprintf(out_file, "light_source { <5, 0, -2> color White}\n");
  fprintf(out_file, "light_source { <-5, 0, -2> color White}\n");
  fprintf(out_file, "light_source { <10, 10, 10> color rgb<0.4, 0.4, 0.4>}\n");
  fprintf(out_file, "cylinder {\n  <0,0,0>, <3, 0, 0>, 0.01\n  pigment { Red }\n}\n");
  fprintf(out_file, "cylinder {\n  <0,0,0>, <0, 3, 0>, 0.01\n  pigment { Green }\n}\n");
  fprintf(out_file, "cylinder {\n  <0,0,0>, <0, 0, 3>, 0.01\n  pigment { Blue }\n}\n");
  
  //got through all the orthants
  for (orth=0; orth<8; orth++) {
    //now, list the vertices
    fprintf(out_file, "mesh2 {\n");
    fprintf(out_file, "\tvertex_vectors {\n");
    fprintf(out_file, "\t\t%d,\n", GLmesh->vertices[orth%4]->num_verts);
    for (i=0; i<GLmesh->vertices[orth%4]->num_verts; i++) {
      if (orth <4) {
        fprintf(out_file, "\t\t<%f,%f,%f>,\n", GLmesh->vertices[orth%4]->verts[i][0],
                                               GLmesh->vertices[orth%4]->verts[i][1],
                                               GLmesh->vertices[orth%4]->verts[i][2]);
      } else {
        fprintf(out_file, "\t\t<%f,%f,%f>,\n", -GLmesh->vertices[orth%4]->verts[i][0],
                                               -GLmesh->vertices[orth%4]->verts[i][1],
                                               -GLmesh->vertices[orth%4]->verts[i][2]);
      } 
    }
    fprintf(out_file, "\t}\n");
    //and the triangles
    fprintf(out_file, "\tface_indices {\n");
    fprintf(out_file, "\t\t%d,\n", GLmesh->triangles[orth%4]->num_tris);
    for (i=0; i<GLmesh->triangles[orth%4]->num_tris; i++) {
      fprintf(out_file, "\t\t<%d,%d,%d>,\n", GLmesh->triangles[orth%4]->tris[i].verts[0],
                                            GLmesh->triangles[orth%4]->tris[i].verts[1],
                                            GLmesh->triangles[orth%4]->tris[i].verts[2]);
                                            
    }
    fprintf(out_file, "\t}\n");
    fprintf(out_file, "pigment {rgb 0.8}\n");
    fprintf(out_file, "}\n");
  }
  fclose(out_file);
  return FALSE;
}

/*****************************************************************************/
/* save the ball as an eps file (not currently working)                      */
/*****************************************************************************/
static gboolean save_eps_button_press(GtkWidget* widget,
                                      GdkEventButton* event,
                                      fieldList* fields) {  
  return FALSE;
}
    
/*****************************************************************************/
/* this is called whenever we move the move (so we can rotate)               */
/*****************************************************************************/
static gboolean drawing_motion_notify(GtkWidget* area,
                                      GdkEventMotion* event,
                                      fieldList* fields) {
  GdkModifierType state;
  int x, y;
  if (EGlobal == NULL) {
    return TRUE;
  }
  if (event->is_hint) {
    gdk_window_get_pointer(event->window, &x, &y, &state);
  } else {
    x = event->x;
    y = event->y;
    state = event->state;
  }
  
  y = area->allocation.height - event->y;
  
  if (state && 
    (GDK_BUTTON1_MASK || GDK_BUTTON3_MASK)) {
    if (GDK_BUTTON1_MASK) {
      dstatus.current_rotation_x = (event->x - dstatus.button1_down_x)/10;
      dstatus.current_rotation_y = (area->allocation.height - event->y - dstatus.button1_down_y)/10;
    }
    draw_mesh();
    gdk_window_invalidate_rect(area->window, &area->allocation, FALSE);
	  gdk_window_process_updates(area->window, FALSE);
  }
  
  return FALSE;
}
    
/*****************************************************************************/
/* called whenever we click on the drawing area                              */
/*****************************************************************************/    
static gboolean drawing_mouse_click(GtkWidget* area,
                                    GdkEventButton* event,
                                    fieldList* fields) {
  if (EGlobal == NULL) {
    return TRUE;
  }
  if (event->button == 1) { //left button
    if (event->type == GDK_BUTTON_PRESS) {
      dstatus.button1_down_x = event->x;
      dstatus.button1_down_y = area->allocation.height - event->y;
    } else {
      dstatus.current_rotation_x = (event->x - dstatus.button1_down_x)/10;
      dstatus.constant_rotation_x += dstatus.current_rotation_x;
      dstatus.current_rotation_x = 0;
      
      dstatus.current_rotation_y = (area->allocation.height - event->y - dstatus.button1_down_y)/10;
      dstatus.constant_rotation_y += dstatus.current_rotation_y;
      dstatus.current_rotation_y = 0;
      //printf("Current rotations: %f, %f\n", dstatus.constant_rotation_x, dstatus.constant_rotation_y);
      draw_mesh();
      gdk_window_invalidate_rect(area->window, &area->allocation, FALSE);
	    gdk_window_process_updates(area->window, FALSE);
	  }
  }
  
  return TRUE;
}

/*****************************************************************************/
/* this gets called whenever an area of the window is invalidated            */
/*****************************************************************************/
static gboolean expose(GtkWidget* area, GdkEventExpose* event, gpointer data) {
  gdk_draw_drawable(area->window,
                    area->style->fg_gc[gtk_widget_get_state (area)],
                    pixmap,
                    event->area.x, event->area.y,
                    event->area.x, event->area.y,
                    event->area.width, event->area.height);
  return FALSE;
}

/*****************************************************************************/
/* called once when gtk starts up                                            */
/*****************************************************************************/
static gboolean configure(GtkWidget* area, 
                          GdkEventConfigure* event,
                          gpointer data) {
                          
  if (pixmap) {
    g_object_unref(pixmap);
  }
  pixmap = gdk_pixmap_new(area->window,
                          area->allocation.width,
                          area->allocation.height,
                          -1);
  GdkGLConfig* glconfig;   
  glconfig = gdk_gl_config_new_by_mode(GDK_GL_MODE_RGB |
                                       GDK_GL_MODE_DEPTH |
                                       GDK_GL_MODE_SINGLE); 
                                       
  //this returns the GdkGLPixmap, but we don't care, because we can get it 
  //out later
  GLPixmap = gdk_pixmap_set_gl_capability(pixmap,
                                           glconfig,
                                           NULL);
  if (GLPixmap == NULL) {
    printf("Couldn't initialize opengl pixmap\n");
  }                                          
  
  GdkGLDrawable* gldrawable = gdk_pixmap_get_gl_drawable(pixmap);
  
  GLContext = gdk_gl_context_new(gldrawable,
                                 NULL,
                                 FALSE,
                                 GDK_GL_RGBA_TYPE);                                        
  
  if (!gdk_gl_drawable_gl_begin(gldrawable, GLContext)) {
    printf("Couldn't start opengl drawing\n");
  }
  glLoadIdentity();
  glViewport(0,0,area->allocation.width, area->allocation.height);
  glOrtho (-4,4,-4,4,-4,4);
	glEnable (GL_BLEND);
	glEnable(GL_DEPTH_TEST);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glEnable(GL_NORMALIZE);
  
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //glEnable(GL_LIGHT1);
  //glEnable(GL_LIGHT2);
  //glEnable(GL_LIGHT3);
  //glEnable(GL_LIGHT4);
  glEnable(GL_LIGHT5);

  GLfloat lightSpec[] = {0.2, 0.2, 0.2, 1};
  GLfloat lightDiff[] = {0.7, 0.7, 0.7, 1};
  GLfloat lightAmbient[] = {0.5, 0.5, 0.5, 1};
  
  //put the lights in place
  GLfloat lightPos0[] = { 0, 0, 1, 0 };
  GLfloat lightPos1[] = { 1, 0, 0, 0 };
  GLfloat lightPos2[] = { 0, 1, 0, 0 };
  GLfloat lightPos3[] = { -1, 0, 0, 0 };
  GLfloat lightPos4[] = { 0, -1, 0, 0 };
  GLfloat lightPos5[] = { 0, 0, -1, 0 };
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
  glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);
  glLightfv(GL_LIGHT2, GL_POSITION, lightPos2);
  glLightfv(GL_LIGHT3, GL_POSITION, lightPos3);
  glLightfv(GL_LIGHT4, GL_POSITION, lightPos4);
  glLightfv(GL_LIGHT5, GL_POSITION, lightPos5);
  
  glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiff);  
  glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiff);  
  glLightfv(GL_LIGHT2, GL_DIFFUSE, lightDiff);  
  glLightfv(GL_LIGHT3, GL_DIFFUSE, lightDiff);  
  glLightfv(GL_LIGHT4, GL_DIFFUSE, lightDiff);
  glLightfv(GL_LIGHT5, GL_DIFFUSE, lightDiff);
  
  glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpec);
  glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpec);
  glLightfv(GL_LIGHT2, GL_SPECULAR, lightSpec);
  glLightfv(GL_LIGHT3, GL_SPECULAR, lightSpec);
  glLightfv(GL_LIGHT4, GL_SPECULAR, lightSpec);
  glLightfv(GL_LIGHT5, GL_SPECULAR, lightSpec);
  
  glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
  glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient);
  glLightfv(GL_LIGHT2, GL_AMBIENT, lightAmbient);
  glLightfv(GL_LIGHT3, GL_AMBIENT, lightAmbient);
  glLightfv(GL_LIGHT4, GL_AMBIENT, lightAmbient);
  glLightfv(GL_LIGHT5, GL_AMBIENT, lightAmbient);
  
        
  //glBegin (GL_LINES);
	//glColor3f (1., 0., 0.);
	//glVertex3f (0., 0., 0.);
	//glVertex3f (1., 0., 0.);
	//glEnd ();
  
  if(gdk_gl_drawable_is_double_buffered(gldrawable)) {
    gdk_gl_drawable_swap_buffers(gldrawable);
  } else {
	  glFlush();
	}
	gdk_gl_drawable_gl_end(gldrawable);
	
	if (VERBOSE)
	  printf("Configured opengl\n");
	
	dstatus.current_rotation_x = 0;
	dstatus.current_rotation_y = 0;
	
	gdk_window_invalidate_rect(area->window, &area->allocation, FALSE);
	gdk_window_process_updates(area->window, FALSE);
	
	return TRUE;
}
  
  
/*****************************************************************************/
/* this makes the X button in the upper right work                           */
/*****************************************************************************/  
static gboolean delete_event(GtkWidget* widget, GdkEvent* event, gpointer data) {
  return FALSE;
}
static void destroy(GtkWidget* widget, gpointer data) {
  gtk_main_quit();
}


int main(int argc, char* argv[]) {


  //main areas
  GtkWidget* window;
  GtkWidget* hBox;
  GtkWidget* control_box;
  GtkWidget* drawing_area;
  
  //buttons and stuff
  GtkWidget* chain1_entry;
  GtkEntryBuffer* chain1_text;  
  GtkWidget* chain2_entry;
  GtkEntryBuffer* chain2_text;
  GtkWidget* chain3_entry;
  GtkEntryBuffer* chain3_text;
  GtkWidget* run_pause_box;
  GtkWidget* run_button;
  GtkWidget* pause_button;
  GtkWidget* step_button;
  GtkWidget* tol_entry_box;
  GtkEntryBuffer* tol_text;
  GtkWidget* tol_label;
  GtkWidget* tolerance_entry;
  //GtkWidget* change_tol_button;
  GtkWidget* file_entry;
  GtkEntryBuffer* file_entry_text;
  GtkWidget* save_pov_button;
  GtkWidget* save_eps_button;
  GtkWidget* save_hbox;
  GtkWidget* dash_m5_check;
  
  gtk_init(&argc, &argv);
    
  gdk_gl_init(&argc, &argv);  //gtkglext
  
  //check if we're getting the -v option for verbosity
  if (argc>1) {
    if (strcmp(argv[1], "-v")==0) {
      VERBOSE = 1;
      printf("Looks like you want verbose output\n");
    }
  }
   
  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  
  //make the main drawing area
  drawing_area = gtk_drawing_area_new();
  gtk_widget_set_size_request(drawing_area, 800, 800); 
  
  //this holds everything
  hBox = gtk_hbox_new(FALSE, 0);
  
  //sidebar
  control_box = gtk_vbox_new(FALSE,0);
  
  //sidebar buttons, etc
  chain1_text = gtk_entry_buffer_new("abAABB ab", -1);
  chain2_text = gtk_entry_buffer_new("abAB", -1);
  chain3_text = gtk_entry_buffer_new("aBAbabAB", -1);
  chain1_entry = gtk_entry_new_with_buffer(chain1_text);
  chain2_entry = gtk_entry_new_with_buffer(chain2_text);
  chain3_entry = gtk_entry_new_with_buffer(chain3_text);
  run_pause_box = gtk_hbox_new(FALSE, 0);
  run_button = gtk_button_new_with_label("run");
  pause_button = gtk_button_new_with_label("pause");
  step_button = gtk_button_new_with_label("step");
  dash_m5_check = gtk_check_button_new_with_label("try harder (-m5)");
  
  tol_entry_box = gtk_hbox_new(FALSE, 0);
  tol_text = gtk_entry_buffer_new(NULL, -1);
  tol_label = gtk_label_new("Tolerance:");
  tolerance_entry = gtk_entry_new_with_buffer(tol_text);
  //change_tol_button = gtk_button_new_with_label("change tol.");
  file_entry_text = gtk_entry_buffer_new("", -1);
  file_entry = gtk_entry_new_with_buffer(file_entry_text);
  save_pov_button = gtk_button_new_with_label("save .pov");
  save_eps_button = gtk_button_new_with_label("save .eps");
  save_hbox = gtk_hbox_new(FALSE, 0);
  
  
  //build the sidebar
  gtk_box_pack_start(GTK_BOX(control_box), chain1_entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), chain2_entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), chain3_entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), dash_m5_check, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(tol_entry_box), tol_label, FALSE, FALSE, 0); 
  gtk_box_pack_start(GTK_BOX(tol_entry_box), tolerance_entry, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), tol_entry_box, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(run_pause_box), run_button, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(run_pause_box), pause_button, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(run_pause_box), step_button, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), run_pause_box, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), file_entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(save_hbox), save_pov_button, TRUE, TRUE, 0);
  //gtk_box_pack_start(GTK_BOX(save_hbox), save_eps_button, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), save_hbox, FALSE, FALSE, 0);
  
  //put them together
  gtk_box_pack_start(GTK_BOX(hBox), drawing_area, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(hBox), control_box, FALSE, FALSE, 0);
  gtk_widget_set_size_request(control_box, 200, 100);
  
  //put it in the window
  gtk_container_add(GTK_CONTAINER(window), hBox);
  
  g_signal_connect(GTK_OBJECT(window), 
                   "delete_event", 
                   G_CALLBACK(delete_event), 
                   NULL);
  g_signal_connect(GTK_OBJECT(window),
                   "destroy",
                   G_CALLBACK(destroy),
                   NULL);
                   
  //connect the drawing area signals
  g_signal_connect(GTK_OBJECT(drawing_area),
                   "configure-event",
                   G_CALLBACK(configure),
                   NULL);
  g_signal_connect(GTK_OBJECT(drawing_area),
                   "expose-event",
                   G_CALLBACK(expose),
                   NULL);
  g_signal_connect(GTK_OBJECT(drawing_area),
                   "button_press_event",
                   G_CALLBACK(drawing_mouse_click),
                   NULL);
  g_signal_connect(GTK_OBJECT(drawing_area),
                   "button_release_event",
                   G_CALLBACK(drawing_mouse_click),
                   NULL);
  g_signal_connect(GTK_OBJECT(drawing_area),
                   "motion_notify_event",
                   G_CALLBACK(drawing_motion_notify),
                   NULL);
  gtk_widget_add_events(drawing_area, GDK_CONFIGURE
                                    | GDK_BUTTON_PRESS_MASK 
                                    | GDK_BUTTON_RELEASE_MASK
                                    | GDK_POINTER_MOTION_MASK 
                                    | GDK_POINTER_MOTION_HINT_MASK);
  
           
  fieldList fields;
  fields.entry1 = GTK_ENTRY(chain1_entry);
  fields.entry2 = GTK_ENTRY(chain2_entry);
  fields.entry3 = GTK_ENTRY(chain3_entry);
  fields.tol_entry = GTK_ENTRY(tolerance_entry);
  fields.m5_check = dash_m5_check;
  fields.file = GTK_ENTRY(file_entry);
  fields.drawing = drawing_area;
  
  //run
  gtk_widget_add_events(run_button, GDK_BUTTON_RELEASE_MASK
                                  | GDK_BUTTON_PRESS_MASK
                                  | GDK_LEAVE_NOTIFY_MASK);
  gtk_widget_add_events(pause_button, GDK_BUTTON_RELEASE_MASK
                                  | GDK_BUTTON_PRESS_MASK
                                  | GDK_LEAVE_NOTIFY_MASK);
  gtk_widget_add_events(step_button, GDK_BUTTON_RELEASE_MASK
                                  | GDK_BUTTON_PRESS_MASK
                                  | GDK_LEAVE_NOTIFY_MASK);
  //buttons
  gtk_signal_connect(GTK_OBJECT(run_button),
                     "button_press_event",
                     G_CALLBACK(run_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(run_button),
                     "button_release_event",
                     G_CALLBACK(run_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(step_button),
                     "button_press_event",
                     G_CALLBACK(step_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(step_button),
                     "button_release_event",
                     G_CALLBACK(step_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(pause_button),
                     "button_press_event",
                     G_CALLBACK(pause_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(pause_button),
                     "button_release_event",
                     G_CALLBACK(pause_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(save_pov_button),
                     "button_press_event",
                     G_CALLBACK(save_pov_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(save_pov_button),
                     "button_release_event",
                     G_CALLBACK(save_pov_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(save_eps_button),
                     "button_press_event",
                     G_CALLBACK(save_eps_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(save_eps_button),
                     "button_release_event",
                     G_CALLBACK(save_eps_button_press),
                     &fields);
  
  
  
  //show everything
  gtk_widget_show_all(window);  
 
  //init the lp
  init_lp();
  
  
  gtk_main();
  
  return 0;
}








