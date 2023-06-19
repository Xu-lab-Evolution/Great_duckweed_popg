/*#################-SCRIPT TO CALCULATE SUMMARY STATISTICS FROM MS-FORMAT DATA-########################*/
/*#################-WRITTEN BY PABLO DUCHEN-########################*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

void readmsFileComp(char *in_file,int **msdata,int *site,int *sample);
int extractPop(int **msdata,int **msdataMod,int sites,int inicio,int fin,int *newsites);
int extract2Pop(int **msdata,int **msdataMod,int sites,int inicio1,int fin1,int inicio2,int fin2,int *newsites);

int segsites(int site);
double thetaW(int sample,int segs);
double thetaPi(int **msdata, int site, int sample);
double tajimasD(double Pi,double thetaW,int segsites,int sample);
int dvk(int **msdata,int site,int sample);
double zns(int **msdata,int site,int sample);
void getJFS(int **dataTwoPop,int sizePop1,int sizePop2,int sites,int **jfs);
void summarizeJFS(int **jfs,int sizePop1, int sizePop2,int *Wone,int *Wtwo,int *Wthree,int *Wfour);
int sfs(int **msdata,int site,int sample,int *site_freq_spec);
double thetaH(int sample,int *site_freq_spec);
double FayWusH(double Pi, double thetaH);
double FuLisD(int sample,int segsites,int singletons);
int sequenceDifferences(int *seq1, int *seq2, int site);
double Pi_between(int **msdata, int site, int sizePop1, int sizePop2);
double da(double Pi_btw,double Pi_pop1, double Pi_pop2,int size_pop1,int size_pop2);

/*#################-MAIN-########################*/

int main(int argc,char **argv) {
    int alto=atoi(argv[2]);
    int largo=atoi(argv[3]);
    int site,sample;
    int i,j,a;
 
    int **data = malloc(sizeof (int*) * alto);
    for (i=0; i<alto; i++) {
      data[i] = malloc(sizeof(int) * largo);
    }
    
    int **data2pop = malloc(sizeof (int*) * alto);
    for (i=0; i<alto; i++) {
      data2pop[i] = malloc(sizeof(int) * largo);
    }
    
    int **dataPop1 = malloc(sizeof (int*) * alto);
    for (i=0; i<alto; i++) {
      dataPop1[i] = malloc(sizeof(int) * largo);
    }
    
    int **dataPop2 = malloc(sizeof (int*) * alto);
    for (i=0; i<alto; i++) {
      dataPop2[i] = malloc(sizeof(int) * largo);
    }
    
    int **jfs = malloc(sizeof (int*) * alto);
    for (i=0; i<alto; i++) {
      jfs[i] = malloc(sizeof(int) * largo);
    }
    
    int *site_freq_spec={0};
    site_freq_spec = (int *)malloc(alto*sizeof(int));
    if (site_freq_spec==NULL) {
      printf("Error allocating memory!\n");
      return 1;
    }
    
    int *site_freq_spec_pop1={0};
    site_freq_spec_pop1 = (int *)malloc(alto*sizeof(int));
    if (site_freq_spec_pop1==NULL) {
      printf("Error allocating memory!\n");
      return 1;
    }

    int *site_freq_spec_pop2={0};
    site_freq_spec_pop2 = (int *)malloc(alto*sizeof(int));
    if (site_freq_spec_pop2==NULL) {
      printf("Error allocating memory!\n");
      return 1;
    }
    
    int new_site,new_site_pop1,new_site_pop2;
    int Wone=0;
    int Wtwo=0;
    int Wthree=0;
    int Wfour=0;

    int site_init_pop1=atoi(argv[4]);
    int site_end_pop1=atoi(argv[5]);
    int site_init_pop2=atoi(argv[6]);
    int site_end_pop2=atoi(argv[7]);

    int new_site_init_pop1=0;
    int new_site_end_pop1=site_end_pop1-new_site_init_pop1;
    int new_site_init_pop2=new_site_end_pop1+1;
    int new_site_end_pop2=new_site_init_pop2+site_end_pop2-site_init_pop2;

    char *tag_pop1=argv[8];
    char *tag_pop2=argv[9];

    int only_pop=atoi(argv[10]); //0 for all and 2 for pop 2.
    int header=atoi(argv[11]); //0 for no header, 1 to print only header, and 2 to print both.
    
    readmsFileComp(argv[1],data,&site,&sample);
    /*printf("\nData set: site = %d, sample = %d\n",site,sample);
    for (i=0;i<sample;i++) {
        for (j=0;j<site;j++) {
            printf("%d",data[i][j]);
        }
        printf("\n");
    }*/
    
    int new_sample=extract2Pop(data,data2pop,site,site_init_pop1,site_end_pop1,site_init_pop2,site_end_pop2,&new_site); //0 16 17 33
    /*printf("\nNew data set: site = %d, sample = %d\n",new_site,new_sample);
    for (i=0;i<new_sample;i++) {
        for (j=0;j<new_site;j++) {
            printf("%d",dataCurr[i][j]);
        }
        printf("\n");
    }*/
    
    int new_sample_pop1=extractPop(data2pop,dataPop1,new_site,new_site_init_pop1,new_site_end_pop1,&new_site_pop1);
    int new_sample_pop2=extractPop(data2pop,dataPop2,new_site,new_site_init_pop2,new_site_end_pop2,&new_site_pop2);

    if (header==1 || header==2) {
      if (only_pop==0) printf("S_tot\t");
      if (only_pop==0) printf("S_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("S_pop_%s\t",tag_pop2);
      
      if (only_pop==0) printf("thetaW_tot\t");
      if (only_pop==0) printf("thetaW_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("thetaW_pop_%s\t",tag_pop2);
      
      if (only_pop==0) printf("Pi_tot\t");
      if (only_pop==0) printf("Pi_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("Pi_pop_%s\t",tag_pop2);
      
      if (only_pop==0) printf("TajD_tot\t");
      if (only_pop==0) printf("TajD_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("TajD_pop_%s\t",tag_pop2);
      
      if (only_pop==0) printf("dvk_tot\t");
      if (only_pop==0) printf("dvk_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("dvk_pop_%s\t",tag_pop2);
      
      if (only_pop==0) printf("ZnS_tot\t");
      if (only_pop==0) printf("ZnS_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("ZnS_pop_%s\t",tag_pop2);
      
      if (only_pop==0 || only_pop==2 || only_pop==9) printf("W1_%s_%s\tW2_%s_%s\tW3_%s_%s\tW4_%s_%s\t",tag_pop1,tag_pop2,tag_pop1,tag_pop2,tag_pop1,tag_pop2,tag_pop1,tag_pop2);
      
      if (only_pop==0) printf("thetaH_tot\t");
      if (only_pop==0) printf("thetaH_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("thetaH_pop_%s\t",tag_pop2);
      
      if (only_pop==0) printf("FayWuH_tot\t");
      if (only_pop==0) printf("FayWuH_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("FayWuH_pop_%s\t",tag_pop2);
      
      if (only_pop==0) printf("FuLiD_tot\t");
      if (only_pop==0) printf("FuLiD_pop_%s\t",tag_pop1);
      if (only_pop==0 || only_pop==2) printf("FuLiD_pop_%s\t",tag_pop2);
      
      if (only_pop==0 || only_pop==2 || only_pop==9) printf("daNei_%s_%s\t",tag_pop1,tag_pop2);
      
      if (only_pop==0) for (a=1;a<new_sample;a++) printf("sfs_tot_%d\t",a);
      if (only_pop==0) for (a=1;a<new_sample_pop1;a++) printf("sfs_pop_%s_%d\t",tag_pop1,a);
      if (only_pop==0 || only_pop==2) for (a=1;a<new_sample_pop2;a++) printf("sfs_pop_%s_%d\t",tag_pop2,a);
      
      if (header==2) printf("\n");
    }

    if (header==0 || header==2) {
      double thetaW_2pop=thetaW(new_sample,new_site);
      double thetaW_pop1=thetaW(new_sample_pop1,new_site_pop1);  
      double thetaW_pop2=thetaW(new_sample_pop2,new_site_pop2);
      
      double thetaPi_2pop=thetaPi(data2pop,new_site,new_sample);
      double thetaPi_pop1=thetaPi(dataPop1,new_site_pop1,new_sample_pop1);
      double thetaPi_pop2=thetaPi(dataPop2,new_site_pop2,new_sample_pop2);
      
      double TajD_2pop=tajimasD(thetaPi_2pop,thetaW_2pop,new_site,new_sample);
      double TajD_pop1=tajimasD(thetaPi_pop1,thetaW_pop1,new_site_pop1,new_sample_pop1);
      double TajD_pop2=tajimasD(thetaPi_pop2,thetaW_pop2,new_site_pop2,new_sample_pop2);
      
      int dvk_2pop=dvk(data2pop,new_site,new_sample);
      int dvk_pop1=dvk(dataPop1,new_site_pop1,new_sample_pop2);
      int dvk_pop2=dvk(dataPop2,new_site_pop2,new_sample_pop2);
      
      double zns_2pop=zns(data2pop,new_site,new_sample);
      double zns_pop1=zns(dataPop1,new_site_pop1,new_sample_pop1);
      double zns_pop2=zns(dataPop2,new_site_pop2,new_sample_pop2);
      
      getJFS(data2pop,new_sample_pop1,new_sample_pop2,new_site,jfs);
      summarizeJFS(jfs,new_sample_pop1,new_sample_pop2,&Wone,&Wtwo,&Wthree,&Wfour);
      
      sfs(data2pop,new_site,new_sample,site_freq_spec);
      sfs(dataPop1,new_site_pop1,new_sample_pop1,site_freq_spec_pop1);
      sfs(dataPop2,new_site_pop2,new_sample_pop2,site_freq_spec_pop2);
      
      double thetaH_2pop=thetaH(new_sample,site_freq_spec);
      double thetaH_pop1=thetaH(new_sample_pop1,site_freq_spec_pop1);
      double thetaH_pop2=thetaH(new_sample_pop2,site_freq_spec_pop2);
      
      double FayWuH_2pop=FayWusH(thetaPi_2pop,thetaH_2pop);
      double FayWuH_pop1=FayWusH(thetaPi_pop1,thetaH_pop1);
      double FayWuH_pop2=FayWusH(thetaPi_pop2,thetaH_pop2);
      
      double FuLiD_2pop=FuLisD(new_sample,new_site,site_freq_spec[1]);
      double FuLiD_pop1=FuLisD(new_sample_pop1,new_site_pop1,site_freq_spec_pop1[1]);
      double FuLiD_pop2=FuLisD(new_sample_pop2,new_site_pop2,site_freq_spec_pop2[1]);
      
      double Pi_btw=Pi_between(data2pop,new_site,new_sample_pop1,new_sample_pop2);
      double daNei=da(Pi_btw,thetaPi_pop1,thetaPi_pop2,new_sample_pop1,new_sample_pop2);
      
      if (only_pop==0) printf("%d\t",segsites(new_site));
      if (only_pop==0) printf("%d\t",segsites(new_site_pop1));
      if (only_pop==0 || only_pop==2) printf("%d\t",segsites(new_site_pop2));
      
      if (only_pop==0) printf("%f\t",thetaW_2pop);
      if (only_pop==0) printf("%f\t",thetaW_pop1);
      if (only_pop==0 || only_pop==2) printf("%f\t",thetaW_pop2);
      
      if (only_pop==0) printf("%f\t",thetaPi_2pop); 
      if (only_pop==0) printf("%f\t",thetaPi_pop1);
      if (only_pop==0 || only_pop==2) printf("%f\t",thetaPi_pop2);
      
      if (only_pop==0) printf("%f\t",TajD_2pop);
      if (only_pop==0) printf("%f\t",TajD_pop1);
      if (only_pop==0 || only_pop==2) printf("%f\t",TajD_pop2);
      
      if (only_pop==0) printf("%d\t",dvk_2pop);  
      if (only_pop==0) printf("%d\t",dvk_pop1);
      if (only_pop==0 || only_pop==2) printf("%d\t",dvk_pop2);
      
      if (only_pop==0) printf("%f\t",zns_2pop);
      if (only_pop==0) printf("%f\t",zns_pop1);
      if (only_pop==0 || only_pop==2) printf("%f\t",zns_pop2);
      
      if (only_pop==0 || only_pop==2 || only_pop==9) printf("%d\t%d\t%d\t%d\t",Wone,Wtwo,Wthree,Wfour);
      
      if (only_pop==0) printf("%f\t",thetaH_2pop);    
      if (only_pop==0) printf("%f\t",thetaH_pop1);
      if (only_pop==0 || only_pop==2) printf("%f\t",thetaH_pop2);
      
      if (only_pop==0) printf("%f\t",FayWuH_2pop); 
      if (only_pop==0) printf("%f\t",FayWuH_pop1);
      if (only_pop==0 || only_pop==2) printf("%f\t",FayWuH_pop2);
      
      if (only_pop==0) printf("%f\t",FuLiD_2pop);
      if (only_pop==0) printf("%f\t",FuLiD_pop1);
      if (only_pop==0 || only_pop==2) printf("%f\t",FuLiD_pop2);
      
      if (only_pop==0 || only_pop==2 || only_pop==9) printf("%f\t",daNei);
      
      if (only_pop==0) for (a=1;a<new_sample;a++) printf("%d\t",site_freq_spec[a]);
      if (only_pop==0) for (a=1;a<new_sample_pop1;a++) printf("%d\t",site_freq_spec_pop1[a]);
      if (only_pop==0 || only_pop==2) for (a=1;a<new_sample_pop2;a++) printf("%d\t",site_freq_spec_pop2[a]);
      
      if (header==2) printf("\n");
    }

    for (i=0;i<alto;i++) free(data[i]);
    free(data);
    for (i=0;i<alto;i++) free(data2pop[i]);
    free(data2pop);
    for (i=0;i<alto;i++) free(dataPop1[i]);
    free(dataPop1);
    for (i=0;i<alto;i++) free(dataPop2[i]);
    free(dataPop2);
    for (i=0;i<alto;i++) free(jfs[i]);
    free(jfs);
    free(site_freq_spec);
    free(site_freq_spec_pop1);
    free(site_freq_spec_pop2);

    return(0);
}


/*=============FUNCTIONS=============*/

/*-------------Read ms file complete-------------*/

void readmsFileComp(char *in_file,int **msdata,int *site,int *sample)
{
  FILE*data_ms;
  int w=0,z=0;
  char ch;
  int count=0;

  //the following loop reads in data_sfs and converts the polytable to an array of integers to be saved in data
  data_ms = fopen ( in_file,"r" );
  while ( (ch = fgetc(data_ms)) != EOF ) {
    if (ch == 'p') count=1;
    else if ((ch != '\n') && (count==1)) continue;
    else if ((ch == '\n') && (count==1)) count=2;
    else if ((ch != '\n') && (count==2)) {
      //data[w][z] = atoi(&ch);
      msdata[w][z] = ch - '0';
      z++;
    }
    else if ((ch == '\n') && (count==2)) {
      w++;
      *site = z;
      z=0;
    }
  }
  fclose(data_ms);
  *sample = w;
  if (count==0) {
    *site=0;
    *sample=0;
  }

}

/*-------------Extract population from ms Data-------------*/

int extractPop(int **msdata,int **msdataMod,int sites,int inicio,int fin,int *newsites)
{
  int n,i,j,k,m,zero,one;

  n=0;
  for (j=0;j<sites;j++) {
    zero=0;
    one=0;
    for (k=inicio;k<=fin;k++) {
      if ((msdata[k][j] == 0) && (zero == 0)) zero=1;
      if ((msdata[k][j] == 1) && (one == 0)) one=1;
      if (zero+one==2) {
          m=0;
          for (i=inicio;i<=fin;i++) {
              msdataMod[m][n] = msdata[i][j];
              m++;
          }
          n++;
          break;
      }
    }
  }

  *newsites = n;

  return(fin-inicio+1);

}

/*-------------Extract 2 populations from ms Data-------------*/

int extract2Pop(int **msdata,int **msdataMod,int sites,int inicio1,int fin1,int inicio2,int fin2,int *newsites)
{

  int a,n,i,j,k,m,zero,one,inicio,fin;
  int ancho=fin2-inicio2+1+fin1-inicio1+1;

  int **msdata2pop = malloc(sizeof (int*) * (ancho*2));
  for (i=0; i<(ancho*2); i++) {
    msdata2pop[i] = malloc(sizeof(int) * sites);
  }

  for (j=0;j<sites;j++) {
    a=0;
    for (i=inicio1;i<=fin1;i++) {
      msdata2pop[a][j]=msdata[i][j];
      a++;
    }
    for (i=inicio2;i<=fin2;i++) {
      msdata2pop[a][j]=msdata[i][j];
      a++;
    }
  }

  /*fprintf(stderr,"inside, inicio1=%d, fin1=%d, inicio2=%d, fin2=%d, ancho=%d\n",inicio1,fin1,inicio2,fin2,ancho);
  for (i=0;i<20;i++) {
    for (j=0;j<sites;j++) fprintf(stderr,"%d",msdata2pop[i][j]);
    fprintf(stderr,"\n");
  }*/

  fin=ancho;
  n=0;
  for (j=0;j<sites;j++) {
    zero=0;
    one=0;
    for (k=0;k<fin;k++) {
      if ((msdata2pop[k][j] == 0) && (zero == 0)) zero=1;
      if ((msdata2pop[k][j] == 1) && (one == 0)) one=1;
      if (zero+one==2) {
          m=0;
          for (i=0;i<fin;i++) {
              msdataMod[m][n] = msdata2pop[i][j];
              m++;
          }
          n++;
          break;
      }
    }
  }

  *newsites = n;

  for (i=0;i<(ancho*2);i++) free(msdata2pop[i]);
  free(msdata2pop);

  return ancho;

}

/*-------------Number of segregating sites S-------------*/

int segsites(int site)
{

  return site;

}

/*-------------ThetaW-------------*/

double thetaW(int sample,int segs)
{
  int b,i;
  double theta_waterson=0;
  double sum_inverses=0;

  for (i=0;i<segs;i++) {
    for (b=1;b<=(sample-1);b++) {
      sum_inverses = sum_inverses + 1/(double)b;
    }
    theta_waterson = theta_waterson + 1/sum_inverses;
    sum_inverses=0;
  }


  return theta_waterson;
}

/*-------------ThetaPi-------------*/

double thetaPi(int **msdata, int site, int sample)
{
  double sum_heterozigosity=0, Pi=0;
  int n,m,count;
  
  for (n=0;n<site;n++) {
    count=0;
    for (m=0;m<sample;m++) {
      if (msdata[m][n]==1) count++;
    }
    sum_heterozigosity = sum_heterozigosity + (2*((double)count/(double)sample)*(1-((double)count/(double)sample)));
  }

  Pi = ((double)sample/((double)sample-1))*sum_heterozigosity;

  return Pi;
}

/*-------------Tajima's D-------------*/

double tajimasD(double Pi,double thetaW,int segsites,int sample)
{
  double sum_inverses=0,sum_squared_inverses=0;
  double a_one,a_two,b_one,b_two,c_one,c_two;
  double denom_TajD,TajD;
  int b;

  for (b=1;b<sample;b++) sum_inverses = sum_inverses + 1/(double)b;
  for (b=1;b<sample;b++) sum_squared_inverses = sum_squared_inverses + 1/pow((double)b,2);
  
  a_one = sum_inverses;
  a_two = sum_squared_inverses;
  b_one = ((double)sample+1)/(3*((double)sample-1));
  b_two = (2*(pow((double)sample,2)+(double)sample+3))/(9*(double)sample*((double)sample-1));
  c_one = b_one-(1/a_one);
  c_two = b_two-(((double)sample+2)/(a_one*(double)sample))+(a_two/pow(a_one,2));
  
  denom_TajD = sqrt((c_one/a_one)*(double)segsites+((c_two/(pow(a_one,2)+a_two))*(double)segsites*((double)segsites-1)));

  TajD = (Pi-thetaW)/denom_TajD;

  return TajD;
}

/*-------------DandVK-------------*/

int dvk(int **msdata,int site,int sample)
{
  double *number={0};
  number = (double *)malloc(sample*sizeof(double));
  if (number==NULL) {
    printf("Error allocating memory!\n");
    return 1;
  }

  double *haplotype;
  haplotype = (double *)malloc(sample*sizeof(double));
  if (haplotype==NULL) {
    printf("Error allocating memory!\n");
    return 1;
  }

  int i,j,b,count=0,test;
  double dandvk;

  for (i=0;i<sample;i++) {
    for (j=0;j<site;j++) {
      number[i] = number[i] + (msdata[i][j] * pow(2,site-j-1));
    }
  }

  haplotype[0] = number[0];
  for (i=1;i<sample;i++) {
    test = 0;
    for (b=0;b<count+1;b++) {
      if (number[i] != haplotype[b]) test++;
    }
    if (test > count) {
      count++;
      haplotype[count] = number[i];
    }
  }

  free(number);
  free(haplotype);
  /*printf("Numbers:\n");
  for (i=0;i<sample;i++) printf("%0.f\n",number[i]);
  printf("\n\nHaplotypes:\n");
  for (i=0;i<sample;i++) printf("%0.f\n",haplotype[i]);
  printf("\n");*/

  dandvk = count+1;
  return dandvk;
}

/*-------------Kelly's ZnS-------------*/
double zns(int **msdata,int site,int sample)
{
  double Pij,pi,pj,D,delta,sum=0,denom_Pij,ZnS;
  int count=0,count_i=0,count_j=0,count_ij=0,site_control=0;

  int *number={0};
  number = (int *)malloc(sample*sizeof(int));
  if (number==NULL) {
    printf("Error allocating memory!\n");
    return 1;
  }

  int a,b,c,d;
  int truesites;

  for (d=0;d<site;d++) {
    //fprintf(stderr,"site: %d\n",d);
    for (c=d+1;c<site;c++) {
      for (b=0;b<sample;b++) number[b] = 0;
      count_i = count_j = count_ij = 0;

      for (a=0;a<sample;a++) {
          if (msdata[a][d] == 1) count_i++;
      }

      if ((count_i==0) || (count_i==sample)) {
          if (c==d+1) site_control++;
          continue;//check if it is correct to do it like this
      }

      pi=(double)count_i/(double)sample;
      //fprintf(stderr,"pi=%f, count_i=%d, samplesize=%d, comparison %d with %d\n",pi,count_i,samplesizes[d],d,c);
      
      for (a=0;a<sample;a++) {
          if (msdata[a][c] == 1) count_j++;
      }
      if ((count_j==0) || (count_j==sample)) continue;//check if it is correct to do it like this
      pj=(double)count_j/(double)sample;
      //fprintf(stderr,"pj=%f, count_j=%d, samplesize=%d, comparison %d with %d\n\n",pj,count_j,samplesizes[c],d,c);

      for (b=0;b<sample;b++) number[b] = number[b] + (msdata[b][d]*10 + msdata[b][c]);
      for (a=0;a<sample;a++) {
          if (number[a] == 11) count_ij++;
      }
//      if (samplesizes[d] < samplesizes[c]) denom_Pij = samplesizes[d];
//      else if (samplesizes[d] > samplesizes[c]) denom_Pij = samplesizes[c];
//      else if (samplesizes[d] == samplesizes[c]) denom_Pij = samplesizes[d];
      denom_Pij = sample;
      Pij=(double)count_ij/denom_Pij;

      D = Pij-(pi*pj);
      delta = pow(D,2)/(pi*(1-pi)*pj*(1-pj));
      //fprintf(stderr,"delta=%f, D=%f, pi=%f, pj=%f\n",delta,D,pi,pj);
      //if (delta > 0.5) {
      //printf ("counti=%d, countj=%d, countij=%d, denom=%f, D=%f, pi=%f, pj=%f, delta=%f, d=%d, c=%d\n",count_i,count_j,count_ij,denom_Pij,D,pi,pj,delta,d+1,c+1);
    //printf ("%d\n",d+1);
      //}
      sum = sum + delta;
      
      //for (b=0;b<sample;b++) printf("number: %d\n",number[b]);
      //printf("pi=%f\npj=%f\nPij=%f\n",pi,pj,Pij);
      //printf("D=%f\ndelta=%f\nsum=%f\nsite comparison:%d with %d\n",D,delta,sum,d,c);
    }
  }

  free(number);
  //fprintf(stderr,"site control = %d, site = %d\n",site_control,site);
  truesites=site-site_control;

  ZnS = 2*pow((truesites*(truesites-1)),-1)*sum;
  //ZnS=(2/(6*(6-1)))*sum; //it doesn't work, why????!!!
  //fprintf(stderr,"ZnS = %f\n",ZnS);
  return ZnS;
}

/*-------------Get the joint site frequency spectrum-------------*/

void getJFS(int **dataTwoPop,int sizePop1,int sizePop2,int sites,int **jfs)
{

  int i,j;
  int countPop1=0,countPop2=0;
  int sizeTotal=sizePop1+sizePop2;

  for (j=0;j<sites;j++) {
    countPop1=0;
    countPop2=0;
    for (i=0;i<sizePop1;i++) {
      if (dataTwoPop[i][j]==1) countPop1++;
    }
    for (i=sizePop1;i<sizeTotal;i++) {
      if (dataTwoPop[i][j]==1) countPop2++;
    }
    jfs[countPop1][countPop2]++;
  }
  
}

/*-------------Summarize the joint site frequency spectrum-------------*/

void summarizeJFS(int **jfs,int sizePop1, int sizePop2,int *Wone,int *Wtwo,int *Wthree,int *Wfour)
{
  int i,j;

  *Wone=0;
  for (i=1;i<=sizePop1-1;i++) *Wone = *Wone + jfs[i][0];
  for (i=1;i<=sizePop1-1;i++) *Wone = *Wone + jfs[i][sizePop2];

  *Wtwo=0;
  for (i=1;i<=sizePop2-1;i++) *Wtwo = *Wtwo + jfs[0][i];
  for (i=1;i<=sizePop2-1;i++) *Wtwo = *Wtwo + jfs[sizePop1][i];

  *Wthree=jfs[sizePop1][0]+jfs[0][sizePop2];

  *Wfour=0;
  for (i=1;i<=sizePop2-1;i++) {
    for (j=1;j<=sizePop1-1;j++) *Wfour = *Wfour + jfs[j][i];
  }

}

/*-------------Get the site frequency spectrum-------------*/

int sfs(int **msdata,int site,int sample,int *site_freq_spec)
{
  //double site_freq_spec[500]={0};
  int *raw_data;
  raw_data = (int *)malloc(site*sizeof(int));
  if (raw_data==NULL) {
    printf("Error allocating memory!\n");
    return 1;
  }

  int n,t,p,s,count,m,max_counter;
  
  for (n=0;n<site;n++) { //counts the number of derived states per site and stores them in the vector: raw_data
    count=0;
    for (m=0;m<sample;m++) {
      if (msdata[m][n]==1) count++;
      raw_data[n]=count;
    }
  }
  max_counter=0;
  for(p=0; p<site; p++){
    int state = raw_data[p];
    if(state) max_counter++;
    site_freq_spec[ state ] ++;
  }

  free(raw_data);
  //for (s=0;s<sample;s++) norm_freq_spec[s]=(site_freq_spec[s]/max_counter);//normalizes the site frequency spectrum
  return 1;

}

/*-------------ThetaH-------------*/

double thetaH(int sample,int *site_freq_spec)
{
  int i;
  double theta_H=0;

  for (i=1;i<sample;i++) {
    theta_H = theta_H + (2*site_freq_spec[i]*pow(i,2))/(sample*(sample-1));
  }

  return theta_H;

}

/*-------------Fay_and_Wu's H-------------*/

double FayWusH(double Pi, double thetaH)
{

  double FayAndWusH;

  FayAndWusH = Pi - thetaH;

  return FayAndWusH;

}

/*-------------Fu_and_Li's D*-------------*/

double FuLisD(int sample,int segsites,int singletons)
{

  double D;
  double a=0,ap=0,b=0,c=0,d=0,u=0,v=0;
  int i;

  for (i=1;i<=sample-1;i++) a = a + (1/(double)i);
  for (i=1;i<=sample;i++) ap = ap + (1/(double)i);
  for (i=1;i<=sample-1;i++) b = b + (1/pow(i,2));
  c = 2*( ( ((double)sample*a) - 2*((double)sample-1) ) / ( ((double)sample-1) * ((double)sample-2) ) );
  d = c + (((double)sample-2)/pow(((double)sample-1),2)) + ((2/((double)sample-1))*((3/2)-(((2*ap)-3)/((double)sample-2))-(1/(double)sample)));
  v = ((pow(((double)sample/((double)sample-1)),2)*b)+(pow(a,2)*d)-(2*(((double)sample*a*(a+1))/(pow(((double)sample-1),2)))))/(pow(a,2)+b);
  u = (((double)sample/((double)sample-1))*(a-(((double)sample)/((double)sample-1))))-v;

  D = (((double)sample/((double)sample-1))*(double)segsites-(a*(double)singletons))/(sqrt((u*(double)segsites)+(v*pow((double)segsites,2))));

  return D;

}

/*-------------Number of sequence differences (ms type)-------------*/
int sequenceDifferences(int *seq1, int *seq2, int site)
{
    int i;
    int count=0;
    
    for (i=0;i<site;i++) {
        if (seq1[i]!=seq2[i]) count++;
    }
    
    return count;
    
}

/*-------------Pi between two populations-------------*/

double Pi_between(int **msdata, int site, int sizePop1, int sizePop2)
{
    int m,n;
    int numDiff=0,numComparisons=0;
    int sizeTotal=sizePop1+sizePop2;
    
    for (n=0;n<sizePop1;n++) {
        for (m=sizePop1;m<sizeTotal;m++) {
            numDiff = numDiff + sequenceDifferences(msdata[n],msdata[m],site);
            numComparisons++;
        }
    }
    
    return ((double)numDiff/(double)numComparisons);
    
}
/*-------------Distance of Nei-------------*/

double da(double Pi_btw,double Pi_pop1, double Pi_pop2,int size_pop1,int size_pop2)
{

  return Pi_btw - (((double)size_pop1/((double)size_pop1+(double)size_pop2))*Pi_pop1 + ((double)size_pop2/((double)size_pop1+(double)size_pop2))*Pi_pop2);

}
