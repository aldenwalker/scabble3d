
#ifndef __scabble3d__
#define __scabble3d__

#include "scl_backend.h"

typedef struct { 
  tri_list* triangles[4];
  vert_list_d* normals[4];
  vert_list_d* vertices[4];
} ball_mesh;

void update_ball_picture_while_running(execution* E);


#endif 
