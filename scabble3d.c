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
  GtkWidget* drawing;
} fieldList;

//these are the opengl globals that I think I need
GdkPixmap* pixmap = NULL;
GdkGLPixmap* GLPixmap = NULL;
GdkGLContext* GLContext = NULL;
ball_mesh* GLmesh = NULL;

//this is the global execution going on 
execution* EGlobal = NULL;


/*****************************************************************************/
/* this function draws the opengl scene                                      */
/*****************************************************************************/
void draw_mesh() {
  int i,j,orth;
  GdkGLDrawable* gldrawable = gdk_pixmap_get_gl_drawable(pixmap);
  
  if (!gdk_gl_drawable_gl_begin(gldrawable, GLContext)) {
    printf("Couldn't start opengl drawing\n");
  }
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  
  GLfloat blueish[] = {0.6, 0.6, 0.9, 1};
  GLfloat redish[] = {0.9, 0.6, 0.6, 1};
  GLfloat spec[] = {0.2, 0.2, 0.2, 1};
  
  glRotatef(15 , 1, 0, -1);
  glRotatef(-45, 0, 1, 0);
  
  glMaterialfv(GL_FRONT, GL_DIFFUSE, blueish);
  glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
  
  //go through the triangles and print them
  //I want to make this more efficient
  glBegin(GL_TRIANGLES);
  for (orth=0; orth<4; orth++) {
    for (i=0; i<GLmesh->triangles[orth]->num_tris; i++) {
      glNormal3f(GLmesh->normals[orth]->verts[i][0], 
                 GLmesh->normals[orth]->verts[i][1], 
                 GLmesh->normals[orth]->verts[i][2]);
      for (j=0; j<3; j++) {
        glVertex3f( GLmesh->vertices[orth]->verts[ GLmesh->triangles[orth]->tris[i].verts[j] ][0],
                    GLmesh->vertices[orth]->verts[ GLmesh->triangles[orth]->tris[i].verts[j] ][1],
                    GLmesh->vertices[orth]->verts[ GLmesh->triangles[orth]->tris[i].verts[j] ][2] );
      }    
    }
  }
  glEnd();
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
/* this is a handy cross product   (only for 3d)                             */
/*****************************************************************************/
void dvector_cross(double* dest, double* a, double* b) {
  dest[0] = a[1]*b[2] - a[2]*b[1];
  dest[1] = a[2]*b[0] - a[0]*b[2];
  dest[2] = a[0]*b[1] - a[1]*b[0];
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
  int offset = 0;
  
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
        continue;
      } else {
        tri_list_add_copy(GLmesh->triangles[i], &T->tris[j]);
        //add the normal
        for (k=0; k<3; k++) {
          for (l=0; l<3; l++) {
            point[k][l] = mpq_get_d(V->verts[T->tris[j].verts[k]].coord[l]); 
          }
          point[k][0] = ((i&1)==1 ? -point[k][0] : point[k][0]);
          point[k][1] = (((i>>1)&1)==1 ? -point[k][1] : point[k][1]);
        }
        for (k=0; k<3; k++) {
          diff1[k] = point[1][k] - point[0][k];
          diff2[k] = point[2][k] - point[0][k];
        }
        dvector_cross(normal, diff1, diff2);
        //if we only swapped one coordinate, we need to 
        //reverse orientation on the normals
        if (i==1 || i==3) {
          normal[0] = -normal[0]; normal[1] = -normal[1]; normal[2] = -normal[2];
        }
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
  
  printf("Called with: %s, %s, %s,\n", arg1, arg2, arg3);
  printf(" = %s, %s, %s, \n", args[0], args[1], args[2]);
  printf("I found chain lens: %d, %d, %d\n", chain_lens[0], chain_lens[1], chain_lens[2]);
  printf("For a total of %d words\n", num_words);
  
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
  
  printf("I got the following chains:\n");
  for (i=0; i<3; i++) {
    for (j=0; j<chain_lens[i]; j++) {
      printf("%s ", chains[i][j]);
    }
    printf("\n");
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
                      target_drawing_area);
  
  EGlobal = E;
  
  for (i=0; i<3; i++) {
    EGlobal->initial_arguments[i] = (char*)malloc((strlen(args[i])+1)*sizeof(char));
    strcpy(EGlobal->initial_arguments[i], args[i]);
  }
  
  execution_print(E);
  
  //initialize the GLmesh
  GLmesh = (ball_mesh*)malloc(sizeof(ball_mesh));
  for (i=0; i<4; i++) {
    GLmesh->triangles[i] = (tri_list*)malloc(sizeof(tri_list));
    GLmesh->vertices[i] = (vert_list_d*)malloc(sizeof(vert_list_d));
    GLmesh->normals[i] = (vert_list_d*)malloc(sizeof(vert_list_d*));
    GLmesh->triangles[i]->num_tris = 0;
    GLmesh->triangles[i]->tris = NULL;
    GLmesh->vertices[i]->num_verts = 0;
    GLmesh->vertices[i]->verts = NULL;
    GLmesh->normals[i]->num_verts = 0;
    GLmesh->normals[i]->verts = NULL;
  }
    
  sem_wait(&(E->read_data_sem));
  update_ball_picture_while_running(E);
  printf("Drew the initial ball\n");
    
  printf("Done computation init -- about to start pthread\n"); 
    
  //start the multi-threadedness
  pthread_t worker_thread;
  pthread_create(&worker_thread, NULL, run_execution, (void*)E);
  
  printf("started computation thread\n");
  
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
  
  if (event->type == GDK_BUTTON_RELEASE) {
    
    e1 = (char*)gtk_entry_get_text(fields->entry1);
    e2 = (char*)gtk_entry_get_text(fields->entry2);
    e3 = (char*)gtk_entry_get_text(fields->entry3);
    
    if (EGlobal == NULL || 
         (strcmp(EGlobal->initial_arguments[0], e1)!=0
          || strcmp(EGlobal->initial_arguments[1], e2)!=0
          || strcmp(EGlobal->initial_arguments[2], e3)!=0)) {
    
       if (strlen(e1) == 0
          || strlen(e2) == 0
          || strlen(e3) == 0) {
        printf("You forgot a chain or something\n");
        return TRUE;
      }
      double tolerance = atof((char*)gtk_entry_get_text(fields->tol_entry));
      load_inputs_and_run((char*)gtk_entry_get_text(fields->entry1),
                          (char*)gtk_entry_get_text(fields->entry2),
                          (char*)gtk_entry_get_text(fields->entry3),
                          tolerance,
                          4,
                          EXLP,
                          fields->drawing);
                          
  
    } else { //we already have an execution -- run it
      //make the tolerance whatever is in the box
      EGlobal->ball->tolerance = atof((char*)gtk_entry_get_text(fields->tol_entry));
      pthread_t worker_thread;
      pthread_create(&worker_thread, NULL, run_execution, (void*)EGlobal);
    
      printf("started computation thread\n");  
    }
  }
  return TRUE;
}       

static gboolean pause_button_press(GtkWidget* widget,
                                   GdkEventButton* event,
                                   fieldList* fields) {
  if (event->type == GDK_BUTTON_RELEASE) {
    sem_wait(&EGlobal->message_sem);
    EGlobal->status_message = 1;
    sem_post(&EGlobal->message_sem);
  }
  return FALSE;
}


static gboolean expose(GtkWidget* area, GdkEventExpose* event, gpointer data) {
  gdk_draw_drawable(area->window,
                    area->style->fg_gc[gtk_widget_get_state (area)],
                    pixmap,
                    event->area.x, event->area.y,
                    event->area.x, event->area.y,
                    event->area.width, event->area.height);
  return FALSE;
}

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
  glOrtho (-3,3,-3,3,-3,3);
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glEnable(GL_NORMALIZE);
  
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_LIGHT2);
  glEnable(GL_LIGHT3);
  glEnable(GL_LIGHT4);
  
  GLfloat lightPos0[] = { 0, 0, 1, 0 };
  GLfloat lightPos1[] = { 1, 0, 0, 0 };
  GLfloat lightPos2[] = { 0, 1, 0, 0 };
  GLfloat lightPos3[] = { -1, 0, 0, 0 };
  GLfloat lightPos4[] = { 0, -1, 0, 0 };
  
  glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);
  glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);
  glLightfv(GL_LIGHT2, GL_POSITION, lightPos2);
  glLightfv(GL_LIGHT3, GL_POSITION, lightPos3);
  glLightfv(GL_LIGHT4, GL_POSITION, lightPos4);
        
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
	
	printf("Configured opengl\n");
	
	gdk_window_invalidate_rect(area->window, &area->allocation, FALSE);
	gdk_window_process_updates(area->window, FALSE);
	
	return TRUE;
}
  
  
  
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
  GtkWidget* tol_entry_box;
  GtkEntryBuffer* tol_text;
  GtkWidget* tolerance_entry;
  GtkWidget* change_tol_button;
  
  gtk_init(&argc, &argv);
    
  gdk_gl_init(&argc, &argv);  //gtkglext
    
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
  tol_entry_box = gtk_hbox_new(FALSE, 0);
  tol_text = gtk_entry_buffer_new(NULL, -1);
  tolerance_entry = gtk_entry_new_with_buffer(tol_text);
  change_tol_button = gtk_button_new_with_label("change tol.");
  
  
  //build the sidebar
  gtk_box_pack_start(GTK_BOX(control_box), chain1_entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), chain2_entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), chain3_entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(run_pause_box), run_button, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(run_pause_box), pause_button, TRUE, TRUE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), run_pause_box, FALSE, FALSE, 0);
  //gtk_box_pack_start(GTK_BOX(tol_entry_box), tolerance_entry, FALSE, FALSE, 0);
  //gtk_widget_set_size_request(tol_entry_box, 20, 40);
  //gtk_box_pack_start(GTK_BOX(tol_entry_box), change_tol_button, FALSE, FALSE, 0);
  //gtk_box_pack_start(GTK_BOX(control_box), tol_entry_box, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), tolerance_entry, FALSE, FALSE, 0);
  gtk_box_pack_start(GTK_BOX(control_box), change_tol_button, FALSE, FALSE, 0);
  
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
           
  fieldList fields;
  fields.entry1 = GTK_ENTRY(chain1_entry);
  fields.entry2 = GTK_ENTRY(chain2_entry);
  fields.entry3 = GTK_ENTRY(chain3_entry);
  fields.tol_entry = GTK_ENTRY(tolerance_entry);
  fields.drawing = drawing_area;
  
  //run
  gtk_widget_add_events(run_button, GDK_BUTTON_RELEASE_MASK
                                  | GDK_BUTTON_PRESS_MASK
                                  | GDK_LEAVE_NOTIFY_MASK);
  gtk_signal_connect(GTK_OBJECT(run_button),
                     "button_press_event",
                     G_CALLBACK(run_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(run_button),
                     "button_release_event",
                     G_CALLBACK(run_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(pause_button),
                     "button_press_event",
                     G_CALLBACK(pause_button_press),
                     &fields);
  gtk_signal_connect(GTK_OBJECT(pause_button),
                     "button_release_event",
                     G_CALLBACK(pause_button_press),
                     &fields);
  
  
  
  //show everything
  gtk_widget_show_all(window);  
 
  //init the lp
  init_lp();
  
  
  gtk_main();
  
  return 0;
}








