/*
 *  word.c
 *  
 *
 *  Created by Alden Walker on 1/28/11.
 *  Copyright 2011 Caltech. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "word.h"


char swapCaseChar(char c) {
  return (char)((int)c > 96 ? (int)c-32 : (int)c+32);
}


void swapCase(char* s) {
  int i;
  int sl = strlen(s);
  for (i=0; i<sl; i++) {
    s[i] = 	(char)(s[i] > 96 ? (int)s[i]-32 : (int)s[i]+32);
  }
}

void  invert(char* s) {
  int i;
  char temp;
  int sl = strlen(s);
  for (i=0; i<sl/2; i++) {
    temp = s[i];
    s[i] = s[sl-i-1];
    s[sl-i-1] = temp;
  }
  swapCase(s);
}

void red(char* s) {			// reduce by cancelling adjacent inverse letters
	int i,a,b;
	i=0;
  int j;
	while(i<=(int)strlen(s)){
		if(i>0){
			a = (int) s[i-1];
			b = (int) s[i];
			if((32+a-b)%64==0){
        //			cout << (*s) << " " << i << "\n";
				for (j=i-1; j<strlen(s)-1; j++) {
          s[j] = s[j+2];
        }
        //			cout << (*s) << "\n";
				i=i-2;
			}
		}
		i++;
	}
	return;
}


void cyc_red(char* s) {
  int sl = strlen(s);
  int i;
  while (1) {
    if (sl == 0) {
      return;
    }
    if (s[0] != swapCaseChar(s[sl-1])) {
      return;
    }
    //remove the first and last
    for (i=0; i<sl; i++) {
      s[i] = s[i+1];
    }
    s[sl-1] = '\0';
    sl -= 2;
  }
}


//string multiply_words(string& s1, string& s2) {
//  string ans = s1+s2;
//  red(ans);
//  return ans;
//}


















