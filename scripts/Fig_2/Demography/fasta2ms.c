/*#################-SCRIPT TO CONVERT FASTA TO MS-FORMAT-########################*/
/*#################-WRITTEN BY PABLO DUCHEN-########################*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <math.h>

void readAlignment_simple(char *in_file,char **alignment,int *site,int *sample,char **names);
int readPositions(char *positions_file, int *positions);


/*#################-MAIN-########################*/
int main(int argc,char **argv) {
    int i,j;
    int largo=500000;
    int largo_names=100;
    int alto=100;
    int site=0;
    int sample=0;
    
    char **alignment = malloc(sizeof (char*) * alto);
    for (i=0; i<alto; i++) {
        alignment[i] = malloc(sizeof(char) * largo);
    }
        
    char **names = malloc(sizeof (char*) * alto);
    for (i=0; i<alto; i++) {
        names[i] = malloc(sizeof(char) * largo_names);
    }
    
    int **msdata = malloc(sizeof (int*) * alto);
    for (i=0; i<alto; i++) {
        msdata[i] = malloc(sizeof(int) * largo);
    }
    
    int *positions={0};
    positions = (int *)malloc(20000*sizeof(int));
    if (positions==NULL) {
        printf("Error allocating memory!\n");
        return 1;
    }
    

    
    readAlignment_simple(argv[1],alignment,&site,&sample,names);
//    for (j=0;j<sample;j++) {
//        printf("%s\n",names[j]);
//        printf("%s\n",alignment[j]);
//    }

//    printf("Site %d, sample %d\n",site,sample);


    int length_positions=readPositions(argv[2],positions);
//    printf("Test positions: %d\n",positions[1]);
//    printf("Length positions: %d\n",length_positions);
    
    int k=0;
    int check;
    for (i=0;i<length_positions;i++) {
        check = 0;
        for (j=0;j<=sample-2;j++) {
            if ( alignment[sample-1][positions[i]]=='N' ) break;
            else if ( alignment[sample-1][positions[i]]==alignment[j][positions[i]] ) msdata[j][k]=0;
            else msdata[j][k]=1;
            check = 1;
        }
        if (check==1) k++;
    }
    
    printf("666 666 666\n");
    printf("\n");
    printf("//\n");
    printf("segsites: %d\n",k);
    printf("positions:\n");
    for (j=0;j<=sample-2;j++) {
        for (i=0;i<k;i++) {
            printf("%d",msdata[j][i]);
        }
        printf("\n");
    }
    
    
    for (i=0;i<alto;i++) free(alignment[i]);
    free(alignment);
    for (i=0;i<alto;i++) free(names[i]);
    free(names);
    for (i=0;i<alto;i++) free(msdata[i]);
    free(msdata);
    free(positions);
    
    
    return(0);
}

/*########################-FUNCTIONS-########################*/

/*-------------Read fasta file (simpler version)-------------*/

void readAlignment_simple(char *in_file,char **alignment,int *site,int *sample,char **names)
{
  FILE*alignment_in;
  int w=0,z=0;
  int m=0,n=0;
  int count=0;
  char ch;

  alignment_in = fopen ( in_file,"r" );
  while ( (ch = fgetc(alignment_in)) != EOF ) {
    if (ch == '>') {
      names[m][n] = ch;
      n++;
      count++;
    }
    else if ((ch != '\n') && (count == 1)) {
      names[m][n] = ch;
      n++;
    }
    else if ((ch == '\n') && (count == 1)) {
      m++;
      n=0;
      count=0;
    }
    else if ((ch != '\n') && (count == 0)) {
      //data[w][z] = atoi(&ch);
      alignment[w][z] = ch;
      z++;
    }
    else if ((ch == '\n') && (count == 0)) {
      w++;
      *site = z;
      z=0;
    }
  }
  fclose(alignment_in);
  *sample = w;
}

/*-------------Read positions file-------------*/

int readPositions(char *positions_file, int *positions)
{

  FILE*file;
  char line[100];
  int j=0;

  file = fopen ( positions_file,"r" );
  j=0;
  while (fgets(line,100,file) != NULL) {
    positions[j]=atoi(line);
    j++;
  }
  
  fclose(file);
  return(j);

}
