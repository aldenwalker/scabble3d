#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include <gtk/gtk.h>
#include "scabble3d.h"
#include "lp.h"




/*****************************************************************************/
/* this function takes the basic input from the text boxes and starts it up  */
/*****************************************************************************/
void load_inputs_and_run(char* arg1, 
                         char* arg2, 
                         char* arg3, 
                         double tolerance,
                         enum scallop_lp_solver solver) {
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
    while (args[i][j] != '\0') {
      if (args[i][j] == ' ') {
        num_words++;
        chain_lens[i] ++;
      }
      j++;
    }
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
      if (args[i][j] == ' ' || args[i][j] == '\0') {
        chains[i][current_word] = (char*)malloc((j-prev_j+1)*sizeof(char));
        strncpy(chains[i][current_word], args[i]+prev_j, (j-prev_j));
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
      }
    }
  }
  
  execution E;
  
  //build the computation -- make initial scl computations, etc
  //init_computation(&E, chains, chain_lens, weights, num_words, tolerance, solver);
  
  //start the multi-threadedness
  //PTHREADs
}
  
  

/*****************************************************************************/
/* this function draws the polygon ball                                      */
/*****************************************************************************/
void update_ball_picture(ball_problem* ball) {
  

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
    
  window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
  
  //make the main drawing area
  drawing_area = gtk_drawing_area_new();
  gtk_widget_set_size_request(drawing_area, 400, 400);
  
  //this holds everything
  hBox = gtk_hbox_new(FALSE, 0);
  
  //sidebar
  control_box = gtk_vbox_new(FALSE,0);
  
  //sidebar buttons, etc
  chain1_text = gtk_entry_buffer_new(NULL, -1);
  chain2_text = gtk_entry_buffer_new(NULL, -1);
  chain3_text = gtk_entry_buffer_new(NULL, -1);
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
  gtk_widget_set_size_request(control_box, 100, 200);
  
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
                   
  //show everything
  gtk_widget_show_all(window);  
 
  
  gtk_main();
  
  return 0;
}
