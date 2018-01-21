#include <stdio.h>
#include <stdlib.h>

void free2d_float(float **p){
  free(p[0]);
  free(p);
}

float **alloc2d_float(int w, int h){
  int   i;
  float **p;
  p = (float **) malloc(h * sizeof(float *));
  if (p == NULL) {
    fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n",
    __FILE__,__LINE__);
    return(p);
  }

  p[0] = (float *) malloc(w * h * sizeof(float));
  if(p[0] == NULL){
    fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n",
    __FILE__,__LINE__);
    free(p);
    return(NULL);
  }

  for(i = 1; i < h; i++){
    p[i] = &(p[0][i*w]);
  }

  return(p);
}


void free2d_short(short **p){
  free(p[0]);
  free(p);
}

short **alloc2d_short(int w, int h){
  int   i;
  short **p;
  p = (short **) malloc(h * sizeof(short *));
  if (p == NULL) {
    fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n",
    __FILE__,__LINE__);
    return(p);
  }

  p[0] = (short *) malloc(w * h * sizeof(short));
  if(p[0] == NULL){
    fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n",
    __FILE__,__LINE__);
    free(p);
    return(NULL);
  }

  for(i = 1; i < h; i++){
    p[i] = &(p[0][i*w]);
  }

  return(p);
}


void free2d_char(char **p){
  free(p[0]);
  free(p);
}

char **alloc2d_char(int w, int h){
  int   i;
  char **p;
  p = (char **) malloc(h * sizeof(char *));
  if (p == NULL) {
    fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n",
    __FILE__,__LINE__);
    return(p);
  }

  p[0] = (char *) malloc(w * h * sizeof(char));
  if(p[0] == NULL){
    fprintf(stderr,"-E- %s line %d: Memory allocation failed.\n",
    __FILE__,__LINE__);
    free(p);
    return(NULL);
  }

  for(i = 1; i < h; i++){
    p[i] = &(p[0][i*w]);
  }

  return(p);
}


#ifdef TESTALLOC2D
void main(void) 
{
    char  **cp;
    float **fp;
    int   i,j;

    cp = alloc2d_char (3,2);
    fp = alloc2d_float(3,2);

    for (j=0; j<2; j++) for (i=0; i<3; i++) {
        fp[j][i] = i+j;
        cp[j][i] = i+j;
    }

    printf("%f %f %f %f\n",fp[0][0],fp[0][2],fp[1][0],fp[1][2]);
    printf("%d %d %d %d\n",cp[0][0],cp[0][2],cp[1][0],cp[1][2]);
}
#endif

