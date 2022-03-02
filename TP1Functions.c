#include "TP1Functions.h"
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include<stdio.h>
#include<ilcplex/cplex.h>

int read_TP1_instance(FILE*fin,dataSet* dsptr)
{
	int rval = 0;

	//Taille des boites
	int V;
	//Nombre d'objets
	int n;
	fscanf(fin,"%d,%d\n",&n,&V);
	dsptr->V = V;
	dsptr->n = n;
	dsptr->size = (int*)malloc(sizeof(int)*n);

	int i;
	for( i = 0 ; i < n ; i++)
		fscanf(fin,"%d\n",&(dsptr->size[i]));

	// fprintf(stderr,"Instance file read, each bin is %d long and there is %d items of lengths:\n",
	// 		V,n);
	// for( i = 0 ; i < n ; i++)
	// 	fprintf(stderr,"%d\t",dsptr->size[i]);
	fprintf(stderr,"\n");


	return rval;
}

int TP1_solve_exact(dataSet* dsptr)
{
	int rval = 0;

	IP_problem* ip_prob_ptr = &(dsptr->master);
	ip_prob_ptr->env = NULL;
	ip_prob_ptr->lp = NULL;
	ip_prob_ptr->env = CPXopenCPLEX (&rval);
	if(rval) fprintf(stderr,"ERROR WHILE CALLING CPXopenCPLEX\n");
	if ( ip_prob_ptr->env == NULL ) 
	{
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (ip_prob_ptr->env, rval, errmsg);
		fprintf (stderr, "%s", errmsg);
		exit(0);	
	}

	//We create the MIP problem
	ip_prob_ptr->lp = CPXcreateprob (ip_prob_ptr->env, &rval, "TP1");
	if(rval) fprintf(stderr,"ERROR WHILE CALLING CPXcreateprob\n");

	rval = CPXsetintparam (ip_prob_ptr->env, CPX_PARAM_DATACHECK, CPX_ON); 
	rval = CPXsetintparam (ip_prob_ptr->env, CPX_PARAM_SCRIND, CPX_ON);

	int n = dsptr->n;
	int*size = dsptr->size;
	int V = dsptr->V;
	int nv = n+ n*n;

	//We fill our arrays
	//Memory
	ip_prob_ptr->nv = nv;
        ip_prob_ptr->x = (double*)malloc(sizeof(double)*nv);
        ip_prob_ptr->cost = (double*)malloc(sizeof(double)*nv);
        ip_prob_ptr->c_type = (char*)malloc(sizeof(char)*nv);
        ip_prob_ptr->up_bound = (double*)malloc(sizeof(double)*nv);
        ip_prob_ptr->low_bound = (double*)malloc(sizeof(double)*nv);
	ip_prob_ptr->var_name = (char**)malloc(sizeof(char*)*nv);

	int i,j,id = 0;
	//Structures keeping the index of each variable
	int*id_y_i = (int*)malloc(sizeof(int)*n);
	int**id_x_ij = (int**)malloc(sizeof(int*)*n);
	for( i = 0 ; i < n ; i++)
		id_x_ij[i] = (int*)malloc(sizeof(int)*n);

	//First the variables yi (bin #i used or not)
	for( i = 0 ; i < n ; i++)
	{
		//We keep the id
		id_y_i[i] = id;

		//We generate the variable attributes
		ip_prob_ptr->x[id] = 0;
		ip_prob_ptr->cost[id] = 1;
		ip_prob_ptr->c_type[id] = 'B';
		ip_prob_ptr->up_bound[id] = 1;
		ip_prob_ptr->low_bound[id] = 0;
		ip_prob_ptr->var_name[id] = (char*)malloc(sizeof(char)*1024);
	        snprintf(       ip_prob_ptr->var_name[id],
        	                1024,
                	        "y_i%d",
                        	i);
		id++;
	}


	/********************************************/
	/******* FILL HERE FOR VARIABLES x_ij *******/
	/********************************************/



	rval = CPXnewcols( ip_prob_ptr->env, ip_prob_ptr->lp, 
			nv, 
			ip_prob_ptr->cost, 
			ip_prob_ptr->low_bound,
			ip_prob_ptr->up_bound,
			ip_prob_ptr->c_type,
			ip_prob_ptr->var_name);
	if(rval)
		fprintf(stderr,"CPXnewcols returned errcode %d\n",rval);



	//Constraints part
        ip_prob_ptr->rhs = (double*)malloc(sizeof(double));
        ip_prob_ptr->sense = (char*)malloc(sizeof(char));
        ip_prob_ptr->rmatbeg = (int*)malloc(sizeof(int));
	ip_prob_ptr->nz = n+1;


        ip_prob_ptr->rmatind = (int*)malloc(sizeof(int)*nv);
        ip_prob_ptr->rmatval = (double*)malloc(sizeof(double)*nv);
	ip_prob_ptr->const_name = (char**)malloc(sizeof(char*));
	ip_prob_ptr->const_name[0] = (char*)malloc(sizeof(char)*1024);

	//We fill what we can 
	ip_prob_ptr->rmatbeg[0] = 0;

	//We generate and add each constraint to the model
	//Bin capacity constraints
	ip_prob_ptr->rhs[0] = 0;
	ip_prob_ptr->sense[0] = 'L';
	for( i = 0 ; i < n ; i++)
	{
		//Constraint name
	        snprintf(       ip_prob_ptr->const_name[0],
        	                1024,
                	        "capacity_bin_i%d",
                        	i);
		id=0;
		//variable y_i coefficient
	        ip_prob_ptr->rmatind[id] = id_y_i[i];
        	ip_prob_ptr->rmatval[id] =  -V;
		id++;
		//variables x_ij coefficients
		for( j = 0 ; j < n ; j++)
		{
		        ip_prob_ptr->rmatind[id] = id_x_ij[i][j];
        		ip_prob_ptr->rmatval[id] =  size[j];
			id++;
		}
		rval = CPXaddrows( ip_prob_ptr->env, ip_prob_ptr->lp, 
			0,//No new column
			1,//One new row
			n+1,//Number of nonzero coefficients
			ip_prob_ptr->rhs, 
			ip_prob_ptr->sense, 
			ip_prob_ptr->rmatbeg, 
			ip_prob_ptr->rmatind, 
			ip_prob_ptr->rmatval,
			NULL,//No new column
			ip_prob_ptr->const_name );
		if(rval)
			fprintf(stderr,"CPXaddrows returned errcode %d\n",rval);

	}


	/****************FILL HERE*******************/
	/*****Each item must be in exactly one bin***/
	/********************************************/


	//We write the problem for debugging purposes, can be commented afterwards
	rval = CPXwriteprob (ip_prob_ptr->env, ip_prob_ptr->lp, "bin_packing.lp", NULL);
	if(rval)
		fprintf(stderr,"CPXwriteprob returned errcode %d\n",rval);

	//We solve the model
	rval = CPXmipopt (ip_prob_ptr->env, ip_prob_ptr->lp);
	if(rval)
		fprintf(stderr,"CPXmipopt returned errcode %d\n",rval);

	rval = CPXsolwrite( ip_prob_ptr->env, ip_prob_ptr->lp, "bin_packing.sol" );
	if(rval)
		fprintf(stderr,"CPXsolwrite returned errcode %d\n",rval);

	//We get the objective value
	rval = CPXgetobjval( ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval) );
	if(rval)
		fprintf(stderr,"CPXgetobjval returned errcode %d\n",rval);

	//We get the best solution found 
	rval = CPXgetobjval( ip_prob_ptr->env, ip_prob_ptr->lp, &(ip_prob_ptr->objval) );
	rval = CPXgetx( ip_prob_ptr->env, ip_prob_ptr->lp, ip_prob_ptr->x, 0, nv-1 );
	if(rval)
		fprintf(stderr,"CPXgetx returned errcode %d\n",rval);

	//We display the solution
	double tolerance = 0.0001;
	int remaining;
	for( i = 0 ; i < n ; i++) 
	{
		id = id_y_i[i];
		if(ip_prob_ptr->x[id] <= 1-tolerance)
			continue;
		remaining = V;
		fprintf(stderr,"Bin #%d is used with volume %d and contains:\n",i,V);
		for( j = 0 ; j < n ; j++)
		{
			id = id_x_ij[i][j];
			if(ip_prob_ptr->x[id] <= 1-tolerance)
				continue;
			remaining -= size[j];
			fprintf(stderr,"\tItem #%d of volume %d (remaining: %d)\n",j,size[j],remaining);
		}
	}









	return rval;
}

int TP1_solve_heuristic(dataSet* dsptr)
{
	printf("------heuristic------\n");
	printf("\n");

	// nombre objets
	int n = dsptr->n;

	// taille boîtes
	int V = dsptr->V;

	// tableau de poids des objets
	int *array = malloc(n * sizeof(int));

	// structure de boîtes
	typedef struct bin {
		int *array;
		int size;
		int nbelt;
	} bin;

	bin* bins = malloc(n * sizeof(bins));

	for (int i=0; i<n; i++) {
		bins[i].array = malloc(n * sizeof(int));
		bins[i].size = 0;
		bins[i].nbelt = 0;
	} 
    
	// fonction pour qsort (a-b: croissant, b-a: décroissant)
	int compare (const void * a, const void * b)
	{
		return ( *(int*)b - *(int*)a );
	}

	// copie des valeurs dans notre tableau
	for (int i=0; i<n; i++) {
		array[i] =  dsptr->size[i];
	}
	
	// tri le tableau dans l'ordre décroissant
	qsort (array, n, sizeof(int), compare);

	printf("nombre objets: %d \n", n);
	printf("tailles boîtes: %d \n", V);
	printf("poids objets: ");

	// affichage poids
	for (int i=0; i<n; i++) {
		printf("%d ", array[i]);
	}

	// stocke le nombre de boîtes nécessaires
	int nbBin = 0;

	printf("\n");
	printf("\n");

	printf("------first fit decreasing------\n");
	printf("\n");

	// parcourir les valeurs du tableau de poids
	// mettre ces valeurs dans tableau de boîtes

	// parcours les valeurs une par une
	for (int i = 0; i < n; i++) {

		// parcours les struct de boîtes une par une
		for (int j=0; j<n; j++) {

			// si on peut ajouter une valeur dans une boîte
			if (array[i] + bins[j].size <= V) {

				// on ajoute la valeur dans le tableau de la boîte
				bins[j].array[bins[j].nbelt] = array[i];

				// on ajoute la valeur au total de la boîte
				bins[j].size += array[i];

				// on incrémente le nombre d'éléments pour 
				// pouvoir ajouter d'autres valeurs en itérant 
				// à partir de cette valeur
				bins[j].nbelt ++;

				// mets à jour le nombre de boîtes
				if (j >= nbBin) {
					nbBin = j+1;
				}

				// sort de la boucle pour itérer à la prochaine
				// valeur et reparcourir les struct de boîtes
				break;
			}
		} 
    }

	// Affichage du résultats

	// s'arrête au nombre de boîte que l'on a
	for (int i=0; i<nbBin; i++) {
		printf("boîte n° %d (taille: %d): ", i+1, bins[i].size);

		// s'arrête au nombre d'élément pour chaque
		// tableau de chaque structure de boîte
		for (int j=0; j<bins[i].nbelt; j++) {
			printf("%d ", bins[i].array[j]); 
		}
		printf("\n");
	} 

	printf("\n");
	printf("------next fit decreasing------\n");
	printf("\n");

	// nouvelles boîtes
	// allocation n*n peut-être trop grande
	bin* bins2 = malloc((n*n) * sizeof(bins2));

	for (int i=0; i<n; i++) {
		bins2[i].array = malloc(n * sizeof(int));
		bins2[i].size = 0;
		bins2[i].nbelt = 0;
	} 

	// on reset notre nombre de boîtes
	nbBin = 0;

	// parcourir les valeurs du tableau de poids
	// mettre ces valeurs dans tableau de boîtes

	// parcours les valeurs une par une
	for (int i = 0; i < n; i++) {

		// parcours les struct de boîtes
		// cette fois on part de la boîte précédente
		for (int j=nbBin; j<n; j++) {

			// si on peut ajouter une valeur dans une boîte
			if (array[i] + bins2[j].size <= V) {

				// on ajoute la valeur dans le tableau de la boîte
				bins2[j].array[bins2[j].nbelt] = array[i];

				// on ajoute la valeur au total de la boîte
				bins2[j].size += array[i];

				// on incrémente le nombre d'éléments pour 
				// pouvoir ajouter d'autres valeurs en itérant 
				// à partir de cette valeur
				bins2[j].nbelt ++;

				// mets à jour le nombre de boîtes
				if (j >= nbBin) {
					nbBin = j;
				}

				// sort de la boucle pour itérer à la prochaine
				// valeur et reparcourir les struct de boîtes
				break;
			}
		} 
    }

	// on rajoute la valeur manquante car on commence à 0
	nbBin += 1;

	// Affichage du résultats

	// s'arrête au nombre de boîte que l'on a
	for (int i=0; i<nbBin; i++) {
		printf("boîte n° %d (taille: %d): ", i+1, bins2[i].size);

		// s'arrête au nombre d'élément pour chaque
		// tableau de chaque structure de boîte
		for (int j=0; j<bins2[i].nbelt; j++) {
			printf("%d ", bins2[i].array[j]); 
		}
		printf("\n");
	} 

	printf("\n");
	printf("------best fit decreasing------\n");
	printf("\n");

	// nouvelles boîtes
	bin* bins3 = malloc(n * sizeof(bins3));

	for (int i=0; i<n; i++) {
		bins3[i].array = malloc(n * sizeof(int));
		bins3[i].size = 0;
		bins3[i].nbelt = 0;
	} 

	// on reset notre nombre de boîtes
	nbBin = 0;

	// parcourir les valeurs du tableau de poids
	// mettre ces valeurs dans tableau de boîtes

	// parcours les valeurs une par une
	for (int i = 0; i < n; i++) {

		// meilleur boîte
		int bestBin = 0;

		// meilleur différence valeur/tailleBoite
		int bestFit = V;

		// parcours les struct de boîtes une par une
		for (int j=0; j<n; j++) {

			// si on peut ajouter une valeur dans une boîte
			if (array[i] + bins3[j].size <= V) {

				// si on trouve un bestFit
				if (bins3[j].size - array[i] <= bestFit) {
					bestBin = j;
				} 
			}
		} 

		// on ajoute la valeur dans le tableau de la boîte
		bins3[bestBin].array[bins3[bestBin].nbelt] = array[i];

		// on ajoute la valeur au total de la boîte
		bins3[bestBin].size += array[i];

		// on incrémente le nombre d'éléments pour 
		// pouvoir ajouter d'autres valeurs en itérant 
		// à partir de cette valeur
		bins3[bestBin].nbelt ++;

		// mets à jour le nombre de boîtes
		if (bestBin >= nbBin) {
			nbBin = bestBin+1;
		}
    }

	// Affichage du résultats

	// Pour enlever les boîtes vides.
	int cpt = 0;

	// s'arrête au nombre de boîte que l'on a
	for (int i=0; i<nbBin; i++) {

		if (bins3[i].size != 0) {

			printf("boîte n° %d (taille: %d): ", i+1-cpt, bins3[i].size);

			// s'arrête au nombre d'élément pour chaque
			// tableau de chaque structure de boîte
			for (int j=0; j<bins3[i].nbelt; j++) {
				printf("%d ", bins3[i].array[j]); 
			}
			printf("\n");

		} else {
			cpt++;
		} 
	}

	return 0;
}


