/*
 *  scl_setup.c
 *  
 *
 *
 */

#include "scl_setup.h"


/*****************************************************************************/
/* these compute the arcs and polygons                                       */
/*****************************************************************************/
void generate_arcs(arc** arc_list, 
                   int* num_arcs,
                   char** word_list,
                   int word_list_len) {
  int word1, word2, i, j;
  arc temp_arc;
  int a, b;
  int* word_lens = (int*)malloc(word_list_len*sizeof(int));
  (*arc_list) = NULL;
  (*num_arcs) = 0;
  for (i=0; i<word_list_len; i++) {
    word_lens[i] = strlen(word_list[i]);
  }
  
  for (word1=0; word1<word_list_len; word1++) {
    for (i=0; i<word_lens[word1]; i++) {
      
      a = (int)word_list[word1][i];
      for (word2=word1; word2<word_list_len; word2++) {
        for (j=(word1==word2 ? i : 0);
             j<word_lens[word2]; 
             j++) {
          b = (int)word_list[word2][j];
          if ( (32+a-b)%64 == 0 ) {     //the letters are inverses
            (*arc_list) = (arc*)realloc((void*)(*arc_list),
                                        (*num_arcs+2)*sizeof(arc));
            temp_arc.first = i;
            temp_arc.last = j;
            temp_arc.first_word = word1;
            temp_arc.last_word = word2;
            (*arc_list)[*num_arcs] = temp_arc;
            temp_arc.first = j;
            temp_arc.last = i;
            temp_arc.first_word = word2;
            temp_arc.last_word = word1;
            (*arc_list)[*num_arcs+1] = temp_arc;
            (*num_arcs) += 2;
          }
        }
      }
      
    }
  }
  
  free(word_lens);
  
}
            
  



void generate_polygons(char** word_list,
                       int word_list_len,
                       arc* arc_list,
                       int num_arcs,
                       polygon** poly_list,
                       int* num_polys,
                       int maxjun) {
  int i,j,k,add,close,appeared,l;
  polygon testpoly;
  int* word_lens = (int*)malloc(word_list_len*sizeof(int));
  (*poly_list) = NULL;
  (*num_polys) = 0;
  for (i=0; i<word_list_len; i++) {
    word_lens[i] = strlen(word_list[i]);
  }
  
  testpoly.arc = (int*)malloc(maxjun*sizeof(int));
  testpoly.num_arcs = 0;
  
  
  for(i=0;i<arc_list_length-1;i++){  // i is initial index of hypothetical polygon
    testpoly.arc[0]=i;
    j=(i+1);
    size=1;
    while(size>0){
      
      add=1;    // haven't yet added
      close=1;  // haven't yet closed
      
      // does arc j glue up to end of last arc of testpoly?
      if(arc_list[j].first_word == arc_list[testpoly.arc[size-1]].last_word && 
         (arc_list[j].first - arc_list[testpoly.arc[size-1]].last-1)%word_lens[arc_list[j].first_word]==0){    
        appeared=1;
        for(k=0;k<size;k++){
          if(j==testpoly.arc[k]){    //has arc j already appeared in testpoly?
            appeared=0;
          }
        }
        
        if(appeared==1){          // if j has not appeared
          testpoly.arc[size]=j;    // then add it to testpoly
          add=0;                  // note that we have added an arc
          
          // does it close up?
          if(arc_list[j].last_word == arc_list[testpoly.arc[0]].first_word && 
             (arc_list[testpoly.arc[0]].first-arc_list[j].last-1)%word_lens[arc_list[j].last_word]==0){    
            testpoly.size = size+1;
            (*poly_list) = (polygon*)realloc((void*)(*poly_list),
                                             (*num_polys+1)*sizeof(polygon));
            (*poly_list)[*num_polys].num_arcs = size+1;
            (*poly_list)[*num_polys].arc = (int*)malloc((size+1)*sizeof(int));
            for (l=0; l<size+1; l++) {
              (*poly_list)[*num_polys].arc[l] = testpoly.arc[l];
            }
            (*num_polys)++;
            close=0;        // note that we have closed a polygon
          }
        }
      }
      
      
      if(add==0 && close==1){
        size++;
        j=i+1;
        if(size>=maxjun){
          j=arc_list_length;
        }
      } else {
        j++;    // increment j
      }
      while(j>=arc_list_length){
        size--;
        j=testpoly.arc[size]+1;
      }
    }
  }
  
}



/*****************************************************************************/
/* creates the basic constraint matrix (the restrictions from the arcs, plus */
/* the constraints saying that the boundary is ag + bh + ck (or just two)    */
/* this function works for 2d *and* 3d (actually, for any d)                 */
/*****************************************************************************/
void create_constraint_matrix(char*** chains, 
                              arc* arc_list,
                              int num_arcs,
                              polygon* polygon_list,
                              int num_polys,
                              int* weights,
                              RatMat* constraints,
                              int* equalityType) {
  
  int arc_list_length = arc_list.size();
  int polygon_list_length = polygon_list.size();
  int numChains = chains.size();
  int totalNumWords = 0;
  int offset;  // as we build the matrix, we use this for sanity
  int i,j,k,l;
  int n;
  int myWordNumber, firstWordNumber;
  //vector<int> equalityType; //-1 is <=, 0 is ==
  int coefFirstWord;
  int coefMyWord;
  
  for (i=0; i<numChains; i++) {
    totalNumWords += chains[i].size();
  }
  
  
  /*
   build the matrix of constraints for the unit ball:
   each column is a polygon, saying how much of it we want
   rows come from:
   - every arc, saying that that arc must appear the same number of times
   as positive and negative
   - the chains, saying that the boundary decomposes as ag + bh
   The last column is the b vector, as in Ax <= b is the polyhedron
   */
  RatMat_init(constraints,  arc_list_length/2 + totalNumWords - numChains, 
              polygon_list_length + 1);
  equalityType.resize(arc_list_length/2 + totalNumWords - numChains);
  
  RatMat_reset_int(constraints, 0);
  
  /*  this section removed since I believe glpk will do it for us in lp.cc
   //each polygon must appear >= 0 times:
   for (i=0; i<polygon_list_length; i++) {
   RatMat_set_int(constraints, i,i, -1,1);
   RatMat_set_int(constraints, i, polygon_list_length, 0,1);
   equalityType[i] = -1;
   }
   */
  
  offset = 0;
  
  //each arc must appear the same number of times coming and going
  for (i=0; i<arc_list_length/2; i++) {  //for each arc
    for (j=0; j<polygon_list_length; j++) {  //for each polygon
      //count the number of times the arc appears in the polygon, with signs
      n=0;
      for (k=0; k<polygon_list[j].size; k++) {  // for each arc in the polygon
        if(polygon_list[j].arc[k]==2*i){
          n=n+1;
        }
        if(polygon_list[j].arc[k]==2*i+1){  //the same arc, but backwards
          n=n-1;
        }
      }
      // n is now the coefficient in the offset+i row and jth column
      RatMat_set_int(constraints, offset+i, j, n, 1);
    }
    equalityType[offset+i] = 0;
    // note that these rows are == 0, which the last entry in the jth row already is
  }
  
  offset += arc_list_length/2;
  
  //each chain must appear as one -- that if we have aw + bv, they must appear
  //as caw + cbv = c(aw+bv) (if the weights are 1, then just the same # of times)
  
  firstWordNumber = 0; //note this is a running total through all words
  for (i=0; i<numChains; i++) { //for each chain
    for (j=1; j<(int)chains[i].size(); j++) {  //for each word past the first one
      myWordNumber = firstWordNumber + j;
      for (k=0; k<polygon_list_length; k++) { //for each polygon
        //the number of times our word appears, divided by its weight is the same
        //as the number of times that the first word appears, divided by its weight
        coefFirstWord = 0;
        coefMyWord = 0;
        for (l=0; l<polygon_list[k].size; l++) {
          if (arc_list[polygon_list[k].arc[l]].first_word == myWordNumber &&
              arc_list[polygon_list[k].arc[l]].first == 0 ) { //for each copy of the word, this is 1
            coefMyWord ++;
          }
          else if (arc_list[polygon_list[k].arc[l]].first_word == firstWordNumber &&
                   arc_list[polygon_list[k].arc[l]].first == 0) {
            coefFirstWord++;
          }
        }
        //the entry to put in is mycoef/myweight - firstcoef/firstweight
        RatMat_set_int(constraints, offset, k, 
                       weights[firstWordNumber]*coefMyWord - weights[myWordNumber]*coefFirstWord,
                       weights[myWordNumber]*weights[firstWordNumber]);
      }
      equalityType[offset] = 0;
      //again, the right hand side is 0, so we don't need to set it
      offset++;
    }
    firstWordNumber += chains[i].size();
  }  
  
}




