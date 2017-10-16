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


int main(int argc, char **argv)
{

	double **config, **x, *L, **y, *r,*lipid_gr, Li[NDIM], *chol_gr, *cross_gr,**comlipids, **comchol;
	double r2, lipid_norm, zavg, zavg_total, chol_norm, cross_norm,Lave[NDIM];

	int i, k, t, *index, n, tt, start, count, j;
	int N, Nlipids, Nlipidbeads, Nconf, Nperlipid, Nchol, Ncholbeads, Nperchol;
	int Nblock, nlog, nlin, step;
	int lipid_pair, chol_pair, cross_pair;
	FILE *ifp, *fp;
	char junk[200], outfile[100];

	char command[100];

	if (argc < 4) {
		printf("Usage: trajectory  number of configurations");
		printf("# of Configurations");
		exit(0);
	}

	lipid_pair=chol_pair=cross_pair=0;

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

	x = (double**) calloc(Nlipidbeads, sizeof(double*)); //coordinates for all lipid beads, for each time-->Nconf
	for (i=0; i < Nlipidbeads; i++){
		x[i] = (double*) calloc(NDIM, sizeof(double));
	}

	comlipids = (double**) calloc(Nlipids, sizeof(double*)); //coordinates for all lipid beads, for each time-->Nconf
	for (i=0; i < Nlipids; i++){
		comlipids[i] = (double*) calloc(NDIM, sizeof(double));
	}

	y = (double**) calloc(Ncholbeads, sizeof(double*)); //coordinates for all lipid beads, for each time-->Nconf
	for (i=0; i < Ncholbeads; i++){
		y[i] = (double*) calloc(NDIM, sizeof(double));
	}

	comchol = (double**) calloc(Nchol, sizeof(double*)); //coordinates for all lipid beads, for each time-->Nconf
	for (i=0; i < Nchol; i++){
		comchol[i] = (double*) calloc(NDIM, sizeof(double));
	}

	L = (double*) calloc(NDIM, sizeof(double));
	// step = (double*)calloc(Nconf, sizeof (double));

	/*allocate lipid_gr*/
	lipid_gr = (double*)calloc((int)(Li[0]/DR), sizeof (double));
	chol_gr = (double*)calloc((int)(Li[0]/DR), sizeof (double));
	cross_gr = (double*)calloc((int)(Li[0]/DR), sizeof (double));

	r = (double*) calloc(NDIM, sizeof(double));
	Lave[0]=Lave[2]=Lave[1]=0;


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
		for (i=0; i<Nlipids; i++){
			comlipids[i][0]=comlipids[i][1]=comlipids[i][2]=0;
		}

		for (i=0; i<Nchol; i++){
			comchol[i][0]=comchol[i][1]=comchol[i][2]=0;
		}
		//printf("TOTAL = %d", count);

		fgets(junk, 200, ifp);

		if (t%nlog==0) {
			Lave[0] += L[0];
			Lave[1] += L[1];
			Lave[2] += L[2];

			zavg_total = 0;
			for (i=0; i<Nlipidbeads; i++){
				zavg_total += x[i][2];
			}

			for (i=0; i<Ncholbeads; i++){
				zavg_total += y[i][2];
			}

			zavg = zavg_total/(double)N;
			/*calculating coms*/

			for (i=0; i<Nlipids; i++){

				for (j=0; j<Nperlipid; j++){
					comlipids[i][0] += x[(i*Nperlipid)+j][0];
					comlipids[i][1] += x[(i*Nperlipid)+j][1];
					comlipids[i][2] += x[(i*Nperlipid)+j][2];
				}
			}
			for (i=0; i<Nlipids; i++){
				for (k=0; k<NDIM; k++){
					comlipids[i][k]=comlipids[i][k]/Nperlipid;
				}
			}

			for (i=0; i<Nchol; i++){

				for (j=0; j<Nperchol; j++){
					comchol[i][0] += y[(i*Nperchol)+j][0];
					comchol[i][1] += y[(i*Nperchol)+j][1];
					comchol[i][2] += y[(i*Nperchol)+j][2];
				}
			}
			for (i=0; i<Nchol; i++){
				for (k=0; k<NDIM; k++){
					comchol[i][k]=comchol[i][k]/Nperchol;
				}
			}


			for (i=0; i<Nlipids-1; i++){

				for (j=i+1; j<Nlipids;j++){

					if (comlipids[i][2] > zavg && comlipids[j][2] > zavg){
						// are we in the same leaf

						for (k=0; k<2; k++){

							r[k]= comlipids[i][k] - comlipids[j][k];

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
						lipid_pair++;


						}

					if (comlipids[i][2] < zavg && comlipids[j][2] < zavg)
					{

						for (k=0; k<2; k++){

							r[k]= comlipids[i][k] - comlipids[j][k];

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
						lipid_pair++;

					}
				}

			}

			for (i=0; i<Nchol-1; i++){

				for (j=i+1; j<Nchol;j++){

					if (comchol[i][2] > zavg && comchol[j][2] > zavg)
					{

						for (k=0; k<2; k++){

							r[k]= comchol[i][k]- comchol[j][k];

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
						chol_pair++;
					}

					if (comchol[i][2] < zavg && comchol[j][2] < zavg)
					{

						for (k=0; k<2; k++){

							r[k]= comchol[i][k]- comchol[j][k];

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
						chol_pair++;
					}

				}
			}

			for (j=0; j<Nchol; j++){

				for (i=0; i<Nlipids;i++){

					if (comlipids[i][2] > zavg && comchol[j][2] > zavg)
					{

						for (k=0; k<2; k++){

							r[k]= comlipids[i][k]- comchol[j][k];

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
						cross_pair++;
					}

					if (comlipids[i][2] < zavg && comchol[j][2] < zavg)
					{

						for (k=0; k<2; k++){

							r[k]= comlipids[i][k]- comchol[j][k];

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
						cross_pair++;
					}

				}
			}







		}
	}


	Lave[0] = Lave[0]/(Nconf/nlog);
	Lave[1] = Lave[1]/(Nconf/nlog);
	Lave[2] = Lave[2]/(Nconf/nlog);
	lipid_norm = (double)(lipid_pair)*PI*DR/2.;
	lipid_norm = (double)(Lave[0])*(double)(Lave[1])/(4*(lipid_norm));

	chol_norm = (double)(chol_pair)*PI*DR/2.;
	chol_norm = (double)(Lave[0])*(double)(Lave[1])/(4*(chol_norm));

	cross_norm = (double)(cross_pair)*PI*DR/2.;
	cross_norm = (double)(Lave[0])*(double)(Lave[1])/(4*(cross_norm));


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


