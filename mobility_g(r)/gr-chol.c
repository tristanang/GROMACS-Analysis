#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NDIM 3
#define DR 0.02
#define bins 200
#define PI 3.141592654
#define DPPC 0
#define CHOL 1
#define dt 0.02
#define anint(x) ((x >= .5) ? (1.0) : (x <= -.5) ? (-1.0) : (0.0))
#define dt 0.02

int main(int argc, char **argv)
{

  double **config, **x, *L, **y, *r,*lipid_gr, Li[NDIM], *chol_gr, *cross_gr;
  double r2, lipid_norm, zavg, zavg_total, chol_norm, cross_norm;

  int i, k, t, *index, n, tt, start, count, j;
  int N, Nlipids, Nlipidbeads, Nconf, Nperlipid, Nchol, Ncholbeads, Nperchol;
  int Nblock, nlog, nlin, step;
  FILE *ifp, *fp;
  char junk[200], outfile[100];
  int counttop = 0;
  int countbottom = 0;

char command[100];
  
  if (argc < 4) {
    printf("Usage: trajectory  number of configurations");
    printf("# of Configurations");
    exit(0);
  }
  
  sscanf(argv[2], "%d", &Nconf);
  sscanf(argv[3], "%d", &nlog);

  sprintf(command, "awk '/CHOL/ {sum+=$2} END {print sum}' %s", argv[4]);
  if ((fp = popen(command, "r")) == NULL) {
    printf("Could not get ending config\n");
    exit(0);
  }
  fscanf(fp, "%d", &Nchol);
  pclose(fp);

  /* open file */
  if ((ifp = fopen(argv[1], "r"))==NULL){
    printf("\nError opening file %s;  Program aborted!\n\n", argv[1]);
    exit(1);
  }
  
  //  L = (double*) calloc(NDIM, sizeof(double));

  /*number of atoms and frame*/
  fscanf(ifp, "%d %d", &N, &step);
  
  for(k=0; k<NDIM; k++)
  {
    fscanf(ifp, "%lf", &Li[k]);
  }
  
  Nperlipid = 12;
  Nperchol = 8;
  Nlipids = (N - 8*Nchol)/12;

  Ncholbeads = Nchol*Nperchol;
  Nlipidbeads = Nlipids*Nperlipid;

  printf("%d\n", Nlipids);

  fclose(ifp);

  /*memory allocating*/

  index = (int*) calloc(N, sizeof(int));


  config = (double**) calloc(N, sizeof(double*));
  for (i=0; i < N; i++){
    config[i] = (double*) calloc(NDIM, sizeof(double));
  }

  x = (double**) calloc(N, sizeof(double*)); //coordinates for all lipid beads, for each time-->Nconf
  for (i=0; i < N; i++){
    x[i] = (double*) calloc(NDIM, sizeof(double));
  }
  
  y = (double**) calloc(N, sizeof(double*)); //coordinates for all lipid beads, for each time-->Nconf
  for (i=0; i < N; i++){
    y[i] = (double*) calloc(NDIM, sizeof(double));
  }

  L = (double*) calloc(NDIM, sizeof(double));
  // step = (double*)calloc(Nconf, sizeof (double));

  /*allocate lipid_gr*/
  lipid_gr = (double*)calloc((int)(Li[0]/DR), sizeof (double));
  chol_gr = (double*)calloc((int)(Li[0]/DR), sizeof (double));
  cross_gr = (double*)calloc((int)(Li[0]/DR), sizeof (double));

  r = (double*) calloc(NDIM, sizeof(double));
   
  
  /*getting coordinates*/
  
  if ((ifp = fopen(argv[1], "r"))==NULL){
    printf("\nError opening file %s;  Program aborted!\n\n", argv[1]);
    exit(1);
    } 
  // while (feof(ifp) == 0){
  
  for (t=0; t<Nconf; t++) {

    fgets(junk, 200, ifp);

    /*box sizes*/    
    
    for(k=0; k<NDIM; k++)
      {
	fscanf (ifp, "%lf", &L[k]);
      }

      

    /*coordinates*/
    
    for(i=0; i<N; i++)
      {
	fscanf (ifp, "%d %lf %lf %lf", &index[i], &config[i][0], &config[i][1], &config[i][2]);
      }

 
    // can i instead put if loop here
    
 /*coordinates for each timestep*/
     for(i=0; i<Nlipidbeads/2; i++){
      for (k=0; k<NDIM; k++){
	x[i][k]=config[i][k];
      }
    }

    for(i=0; i<Ncholbeads/2; i++){
      for (k=0; k<NDIM; k++){
	y[i][k]=config[i+Nlipidbeads/2][k];
      }
    }

    for(i=0; i<Nlipidbeads/2; i++){
      for (k=0; k<NDIM; k++){
	x[i+Nlipidbeads/2][k]=config[i+Nlipidbeads/2+Ncholbeads/2][k];
      }
    }

    for(i=0; i<Ncholbeads/2; i++){
      for (k=0; k<NDIM; k++){
	y[i+Ncholbeads/2][k]=config[i+Nlipidbeads+Ncholbeads/2][k];
      }
    }
    
    //printf("TOTAL = %d", count);
    
    fgets(junk, 200, ifp);
    
    if (t%nlog==0) {
      
      zavg_total = 0;
      for (i=0; i<Nlipidbeads; i+=12){
	zavg_total += x[i][2];
      }
      zavg = zavg_total/(double)Nlipids;
      
	counttop=0;
/*calculating gr*/
   for (i=0; i<Nlipids-1; i++){
     
     for (j=i+1; j<Nlipids;j++){
       
       if (x[i*12][2] > zavg && x[j*12][2] > zavg)
       // are we in the same leaf
		{
	 
	 for (k=0; k<2; k++){
	   
	   r[k]= x[i*12][k]- x[j*12][k];
	   
	    if(r[k] > (0.5*L[k]))
	      {
		r[k] = r[k] - L[k];
	      }
	    
	    else if(r[k] < (-0.5*L[k]))
	      {
		r[k] = r[k] + L[k];
	      }  
	 }
	 
	 r2 = sqrt(r[0]*r[0] + r[1]*r[1]);
	  //	printf("%lf\n", r2);
	 lipid_gr[(int)(r2/DR)]++;
	counttop++;

   
       }
     }
   }

	counttop=0;
 for (i=0; i<Nchol-1; i++){                 //loop over chol
       
     for (j=i+1; j<Nchol;j++){             //loop over lipids
       
       if (y[i*8][2] > zavg && y[j*8][2] > zavg)
       // are we in the same leaf
		{
	 
	 for (k=0; k<2; k++){
	   
	   r[k]= y[i*8][k]- y[j*8][k];
	   
	    if(r[k] > (0.5*L[k]))
	      {
		r[k] = r[k] - L[k];
	      }
	    
	    else if(r[k] < (-0.5*L[k]))
	      {
		r[k] = r[k] + L[k];
	      }  
	 }
	 
	 r2 = sqrt(r[0]*r[0] + r[1]*r[1]);
	  //	printf("%lf\n", r2);
	 chol_gr[(int)(r2/DR)]++;
	counttop++;
		}
	 }
 }
  
for (i=0; i<Nchol; i++){                 //loop over chol
       
     for (j=0; j<Nlipids;j++){             //loop over lipids
       
       if (x[i*12][2] > zavg && y[j*8][2] > zavg)
       // are we in the same leaf
		{
	 
	 for (k=0; k<2; k++){
	   
	   r[k]= x[i*12][k]- y[j*8][k];
	   
	    if(r[k] > (0.5*L[k]))
	      {
		r[k] = r[k] - L[k];
	      }
	    
	    else if(r[k] < (-0.5*L[k]))
	      {
		r[k] = r[k] + L[k];
	      }  
	 }
	 
	 r2 = sqrt(r[0]*r[0] + r[1]*r[1]);
	  //	printf("%lf\n", r2);
	 cross_gr[(int)(r2/DR)]++;
		}
	 }
 }  
    
   
    
  


    
  }
  }
  

  lipid_norm = (double)((Nlipids/2.) * ((Nlipids/2.)-1))*PI*DR;
  lipid_norm = (double)(Li[0])*(double)(Li[1])/((lipid_norm)*(double)(Nconf/nlog));
  chol_norm = (double)((Nchol/2.) * ((Nchol/2.)-1))*PI*DR;
  chol_norm = (double)(Li[0])*(double)(Li[1])/((chol_norm)*(double)(Nconf/nlog));
  cross_norm = (double)((Nlipids/2.) * (Nchol/2.))*PI*DR;
  cross_norm = (double)(Li[0])*(double)(Li[1])/((cross_norm)*(double)(Nconf/nlog));
  
  sprintf(outfile, "gr-lipid.dat");
  fp = fopen(outfile, "w");
  for (i=0; i<(int)(Li[0]/(2.0*DR)); i++){
    
    fprintf(fp, "%lf\t%lf\n", ((double)i+0.5)*DR,  lipid_norm*lipid_gr[i]/(((double)i+0.5)*DR));  
    
  }

  sprintf(outfile, "gr-chol.dat");
  fp = fopen(outfile, "w");
  for (i=0; i<(int)(Li[0]/(2.0*DR)); i++){
    
    fprintf(fp, "%lf\t%lf\n", ((double)i+0.5)*DR,  chol_norm*chol_gr[i]/(((double)i+0.5)*DR));  
    
  }
  
  sprintf(outfile, "gr-cross.dat");
  fp = fopen(outfile, "w");
  for (i=0; i<(int)(Li[0]/(2.0*DR)); i++){
    
    fprintf(fp, "%lf\t%lf\n", ((double)i+0.5)*DR,  cross_norm*cross_gr[i]/(((double)i+0.5)*DR));  
    
  }

  fclose(fp);
   
  }  
  

