/* test.c: test program for STLS 
   .\test i - test example i, 1 <= i <= 24
   .\test   - test all examples            */ 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[])
{
  char fname[50], cmd[60];
  
  if (argc < 3) {
    printf("Usage: %s <m> <n> [<testnum> <stepm> <stepn> <m2> <n2>]\n", argv[0]);
    return 0;
  }
  
  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int testnum = argc > 3 ? atoi(argv[3]) : 10;
  int stepm = argc > 4 ? atoi(argv[4]) : 0;
  int stepn = argc > 5 ? atoi(argv[5]) : 1;
  //int weighted = argc > 6 ? atoi(argv[6]) : 0;
  int m2 = argc > 6 ? atoi(argv[6]) : 0;
  int n2 = argc > 7 ? atoi(argv[7]) : 0;
  
  int i, j;
  
  for (i = 1; i <= testnum; i++, m += stepm, n += stepn) {
    sprintf(fname, "s%d.txt", i); 

    
    FILE *F = fopen(fname, "wt");
    
    if (F != NULL) {
      fprintf(F, "%d %d %d %d %d\n", 1 + (n2 > 0), 1 + (m2 > 0), m + m2, m + m2 - 1, 0);
      fprintf(F, "%d", n);
      if (n2 > 0) {
        fprintf(F, " %d", n2);
      }
      fprintf(F, "\n");
      
      fprintf(F, "%d", m);
      if (m2 > 0) {
        fprintf(F, " %d", m2);
      }

      fprintf(F, "\n");
      
          
           
      sprintf(cmd, "/bin/cp -p \'p1.txt\' \'p%d.txt\'", i);
      system(cmd);   
      fclose(F);
    }
    
    if (m2 > 0) {
      m2 += stepm;
    }
    if (n2 > 0) {
      n2 += stepn;
    }
  }
  
  
 return 0;
}
