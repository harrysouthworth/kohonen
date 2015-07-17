#include "S.h"
#include "sconnect.h"
#include "Kohonen.h"
#include <math.h>

#define EPSILON 1e-4

/* Construct the coordinate matrix
**    centers		: centered coordinates' matrix of the grid
**    n				: nrow(coord)
**    x_template	: prototype 'x' centered on 0.0
**    y_template	: prototype 'y' centered on 0.0
**    out			: coordinate matrix of the grid
**    p				: polygon size ( 4 for rectangular, 7 for hexagonal )

** Return coordinates' matrix of the grid
*/
void get_coordinates(double* centers, long* n, double* x_template, double* y_template, double* out, long* p)
{
int index = 0;
long iterx, itery;
for(iterx = 0 ; iterx < *n; iterx ++){
	for(itery = 0 ; itery < *p; itery ++){
		out[index * 2] = x_template[ itery ] + centers[ iterx * 2 ];
		out[index * 2 + 1] = y_template[ itery ] + centers[ iterx * 2 + 1 ];
		index++;
		}
	out[ index * 2 ] = -1;
	out[ index * 2 + 1 ] = -1;
	index++;
	}
return;
}


/* Construct the sum of squares vector for each individual regarding to its class
**    data			: data matrix
**    nrowdata		: nrow(data)
**    ncol			: ncol(data)
**    weights		: weights matrix
**    nrowweights	: nrow(nrowweights)
**    groups		: class vector for each individual
**    ssvalues		: Vector where SS are stored

** Return the sum of squares vector
*/
void get_individual_ss(double* data, long* nrowdata, long* ncol, double* weights, long *nrowweights, long *groups, double * ssvalues)
{
double tempsums;
long nbcol = *ncol;
long iterx, itery;
for(iterx = 0 ; iterx < *nrowdata; iterx ++)
	{
	for(itery = 0 ; itery < nbcol ; itery ++)
		{
		tempsums = data[iterx * nbcol + itery] - weights[(groups[iterx] - 1) * nbcol + itery];
		ssvalues[iterx] = ssvalues[iterx] + (tempsums * tempsums);
		}
	}
return;
}


/* Construct the intra sum of squares vector for each cluster
**    weights	: weights matrix
**    nrow		: nrow(weights)
**    ncol		: ncol(weights)
**    g			: vecteur centre
**    ssvalues	: Vector where intra SS are stored

** Return the intra sum of squares vector
*/
void get_inter_ss(double* weights, long* nrow, long* ncol, double* g, double * ssvalues)
{
double tempsums;
long nbcol = *ncol;
long iterx, itery;
for( iterx = 0 ; iterx < *nrow; iterx ++)
	{
	for( itery = 0 ; itery < nbcol ; itery ++)
		{
		tempsums = weights[iterx * nbcol + itery] - g[itery];
		ssvalues[iterx] = ssvalues[iterx] + (tempsums * tempsums);
		}
	}
return;
}


/* Construct a distance matrix between each row and all weights
**    data		: data matrix
**    nrow		: nrow(data)
**    ncol		: ncol(data)
**    weights	: weights matrix
**    nclass	: nrow(weights)
**    result	: distance matrix

** Return the distance matrix
*/
void rows_dist_weights(double* data, long* nrow, long* ncol, double* weights, long* nclass, double *result)
{
double tempsums;
long iterind, iterclass;
long itervar;

for(iterind = 0 ; iterind < *nrow ; iterind ++)
	{
	for(iterclass = 0 ; iterclass < *nclass ; iterclass ++)
		{
		for(itervar = 0 ; itervar < *ncol ; itervar ++)
			{
			tempsums = data[iterind * *ncol + itervar] - weights[iterclass * *ncol + itervar];
			result[iterind * *nclass + iterclass] = result[iterind * *nclass + iterclass] + (tempsums * tempsums);
			}
		}
	}
return;
}


void arrangedist(long* idx, long* idy, long* classx, long* classy, double* data, double* weights, 
				long* n, long* ncol, double* res)
{
long j, id;
long id1, id2, cl1, cl2;
double tempDist, tempSum;

for(id = 0; id < *n; id ++)
	{
	id1 = idx[id];
	id2 = idy[id];
	cl1 = classx[id];
	cl2 = classy[id];
	tempDist= 0.0;
	
	for(j = 0; j < *ncol; j ++)
		{
		tempSum = data[id1 * *ncol + j] - data[id2 * *ncol + j];
		tempDist = tempDist + (tempSum * tempSum);
		}
	res[id * 2] = sqrt(tempDist);
	tempDist = 0.0;
	for(j = 0; j < *ncol; j ++)
		{
		tempSum = weights[cl1 * *ncol + j] - weights[cl2 * *ncol + j];
		tempDist = tempDist + (tempSum * tempSum);
		}
	res[id * 2 + 1] = sqrt(tempDist);
	}
return;
}



/* Construct a weights matrix
**    dataset			: data matrix
**    nrow				: nrow(data)
**    ncol				: ncol(data)
**    weights			: weights matrix
**    nclass			: nrow(weights)
**    alpha				: alpha vector of length ntime
**    rayon				: rayon vector of length ntime
**    ntime				: number of iterations (learning steps)
**    dmcc				: distance matrix between all coordinates of centers on the grid
**    order				: row's indexes of data to present for each learning step
**    dist_unit			: a double that represent the dist between coordinates of 2 adjacent centers
**    neightbourhood	: type of neightbourhood ( 0 = uniform ; 1 = gaussian )
**    weights_eii		: weights (denominator) vector for calculating extended intra-inertia
**    winners_of_steps	: indexes of each winner for each learning step
**    temp_ss			: sum of squares vector used for calculation of extended intra-inertia
**    index_iter_eii	: indexes of each iteration where extended intra-inertia have to be calculated
**    ss_intra			: extended intra-inertia vector for energy plot

** Return the updated weights matrix
*/
void KohonenC(double* dataset, long* nrow, long* ncol, double* weights, long* nclass, 
			  double* alpha, long* rayon, long* ntime, double *dmcc, long *order, double *dist_unit, 
			  long *neightbourhood, double *weights_eii, long *winners_of_steps, double *temp_ss, long *index_iter_eii, double *ss_intra)
{
// number of lines
int n = (int)*nrow;
// number of cols
int np = (int)*ncol;
// number of classes
int nc = (int)*nclass ;

// index of the data line to present for a learning step
int id_current_line;
// current position in ss_intra
int id_ss_intra = 0; 
// current position in the weights matrix
int id_cluster;
// current position in the data matrix
int id_data_line;

double temp_diff_value;
double temp_min_distance = 0.0;
double temp_current_distance;
double temp_h_value;
double temp_ref_max_distance = 0.0;
int id_winner  = -1;
long iteration;

int j;
for(iteration = 0 ; iteration < *ntime ; iteration ++)
	{
	id_current_line = (int)order[ iteration ] ; // get id of line to present to the network
	temp_min_distance = 0.0 ; // get id of line to present to the network

	/*
	* search the winner 
	*/
	for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
		temp_current_distance = 0.0;
		for(j = 0 ; j < np ; j ++){
			temp_diff_value = dataset[( np * id_current_line + j)] - weights[(np*id_cluster + j)];
			temp_current_distance = temp_current_distance + (temp_diff_value * temp_diff_value);
			}
		temp_current_distance = sqrt(temp_current_distance);

		if(id_cluster==0 || temp_current_distance <= temp_min_distance){
			temp_min_distance = temp_current_distance;
			id_winner = id_cluster;            
			}
		}
	/*
	* set the maximum valid distance on the grid for an update
	*/
	temp_ref_max_distance = *dist_unit * rayon[ iteration ] + EPSILON;

	/*
	*
	* browse each cluster, if winner or in the neightbourhood, update codes
	*/
	for( id_cluster = 0 ; id_cluster < nc ; id_cluster ++){

		if(dmcc[id_winner * nc + id_cluster] < temp_ref_max_distance){
			if((int)rayon[iteration] > 0 && *neightbourhood == 1){
				temp_h_value = dmcc[ id_winner * nc + id_cluster];
				temp_h_value = temp_h_value/( 2 * temp_ref_max_distance * temp_ref_max_distance );
				temp_h_value = exp( - temp_h_value );
				}
			else
				temp_h_value = 1;

			for(j = 0 ; j < np ; j ++)
				{
				temp_diff_value = dataset[ ( np * id_current_line ) + j ] - weights[( id_cluster * np ) + j];
				weights[(id_cluster * np) + j] += temp_diff_value * alpha[iteration] * temp_h_value;
				}
			}
		}


		
	/*
	* THIS PART IS EVALUATED IF ( index_iter_eii != -1 )
	* AND the learning step iteration == a value of index_iter_eii
	*/
	if( index_iter_eii[ id_ss_intra ] == iteration ) {
		/*
		* initialization
		*/
		for ( id_cluster = 0 ; id_cluster < nc ; id_cluster ++) {
			temp_ss[id_cluster] = 0.0;
			weights_eii[id_cluster] = 0.0;
			}

		/*
		* predict classes with data in winners_of_steps
		*/
		for(id_data_line = 0 ; id_data_line < n ; id_data_line ++){
			for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
				temp_current_distance = 0.0;

				for(j = 0 ; j < np ; j ++){
					temp_diff_value = dataset[(np*id_data_line + j)] - weights[(np*id_cluster + j)];
					temp_current_distance = temp_current_distance + (temp_diff_value * temp_diff_value);
					}

				temp_current_distance = sqrt(temp_current_distance);

				if(id_cluster==0 || temp_current_distance < temp_min_distance + EPSILON ){
					temp_min_distance = temp_current_distance;
					id_winner = id_cluster;                               
					}
				}
			winners_of_steps[id_data_line] = id_winner;
			}

		/*
		* prepare calculation of ss_intra
		*/
		for ( id_data_line = 0 ; id_data_line < n ; id_data_line ++){
			for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
				
				if (dmcc[winners_of_steps[id_data_line] * nc + id_cluster] <= temp_ref_max_distance + EPSILON){
					
					temp_current_distance = 0.0;
					
					for(j = 0 ; j < np ; j ++){
						temp_diff_value = dataset[(np*id_data_line + j)] - weights[(np*id_cluster + j)];
						temp_current_distance = temp_current_distance + (temp_diff_value * temp_diff_value); 
						}
					
					if((int)rayon[iteration] > 0 && *neightbourhood==1){
						temp_h_value = dmcc[winners_of_steps[id_data_line] * nc + id_cluster];
						temp_h_value = temp_h_value/(2 * temp_ref_max_distance * temp_ref_max_distance);
						temp_h_value = exp(-temp_h_value);
						}
					else
						temp_h_value = 1;

					temp_ss[id_cluster] = temp_ss[id_cluster] + temp_h_value * temp_current_distance;
					weights_eii[id_cluster] = weights_eii[id_cluster] + temp_h_value;
					}
				}
			}
		/*
		* calculation of ss_intra
		*/
		for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++)
			{
			if(weights_eii[id_cluster] != 0)
				ss_intra[id_ss_intra] = ss_intra[id_ss_intra] + temp_ss[id_cluster] / weights_eii[id_cluster];
			}
		ss_intra[id_ss_intra] = ss_intra[id_ss_intra] / nc;
		id_ss_intra = id_ss_intra + 1;
		}

	}


	
	// not really clever, there should be a function for calculation of ss_intra
	for ( id_cluster = 0 ; id_cluster < nc ; id_cluster ++) {
		temp_ss[id_cluster] = 0.0;
		weights_eii[id_cluster] = 0.0;
		}

	/*
	* predict classes with data in winners_of_steps
	*/
	for(id_data_line = 0 ; id_data_line < n ; id_data_line ++){
		for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
			temp_current_distance = 0.0;

			for(j = 0 ; j < np ; j ++){
				temp_diff_value = dataset[(np*id_data_line + j)] - weights[(np*id_cluster + j)];
				temp_current_distance = temp_current_distance + (temp_diff_value * temp_diff_value);
				}

			temp_current_distance = sqrt(temp_current_distance);

			if(id_cluster==0 || temp_current_distance < temp_min_distance + EPSILON ){
				temp_min_distance = temp_current_distance;
				id_winner = id_cluster;                               
				}
			}
		winners_of_steps[id_data_line] = id_winner;
		}

	/*
	* prepare calculation of ss_intra
	*/
	for ( id_data_line = 0 ; id_data_line < n ; id_data_line ++){
		for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
			
			if (dmcc[winners_of_steps[id_data_line] * nc + id_cluster] <= temp_ref_max_distance + EPSILON){
				
				temp_current_distance = 0.0;
				
				for(j = 0 ; j < np ; j ++){
					temp_diff_value = dataset[(np*id_data_line + j)] - weights[(np*id_cluster + j)];
					temp_current_distance = temp_current_distance + (temp_diff_value * temp_diff_value); 
					}
				
				if((int)rayon[iteration] > 0 && *neightbourhood==1){
					temp_h_value = dmcc[winners_of_steps[id_data_line] * nc + id_cluster];
					temp_h_value = temp_h_value/(2 * temp_ref_max_distance * temp_ref_max_distance);
					temp_h_value = exp(-temp_h_value);
					}
				else
					temp_h_value = 1;

				temp_ss[id_cluster] = temp_ss[id_cluster] + temp_h_value * temp_current_distance;
				weights_eii[id_cluster] = weights_eii[id_cluster] + temp_h_value;
				}
			}
		}
	/*
	* calculation of ss_intra
	*/
	for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++)
		{
		if(weights_eii[id_cluster] != 0)
			ss_intra[id_ss_intra] = ss_intra[id_ss_intra] + temp_ss[id_cluster] / weights_eii[id_cluster];
		}
	ss_intra[id_ss_intra] = ss_intra[id_ss_intra] / nc;


return;
}



/* Construct a weights matrix
**    dataset			: data matrix
**    nrow				: nrow(data)
**    ncol				: ncol(data)
**    weights			: weights matrix
**    nclass			: nrow(weights)
**    rayon				: rayon vector of length ntime
**    ntime				: number of iterations (learning steps)
**    dmcc				: distance matrix between all coordinates of centers on the grid
**    dist_unit			: a double that represent the dist between coordinates of 2 adjacent centers
**    size				: size of each class
**    winners_of_steps	: indexes of each winner for each learning step
**    neightbourhood	: type of neightbourhood ( 0 = uniform ; 1 = gaussian )
**    weights_eii		: weights (denominator) vector for calculating extended intra-inertia
**    temp_ss			: sum of squares vector used for calculation of extended intra-inertia
**    ss_intra			: extended intra-inertia vector for energy plot

** Return the updated weights matrix
*/
void KohonenCBatch(double* dataset, long* nrow, long* ncol, double* weights, long* nclass, 
				long* rayon, long* ntime, double *dmcc, double *dist_unit, 
				long *size, long * winners_of_steps, long *neightbourhood, double *weights_eii, double *temp_ss, double *ss_intra)
{
// number of lines
int n = (int)*nrow;
// number of cols
int np = (int)*ncol;
// number of classes
int nc = (int)*nclass ;

// n (denominator) vector for calculating mean of each group
double extendsize;
// index of the winner for a learning step
int id_winner = -1;
// current position in the weights matrix
int id_cluster;
// current position in the data matrix
int id_data_line;

int j;


double temp_current_distance;
double temp_diff_value;
double temp_min_distance = 0.0;
double temp_ref_max_distance = 0.0;
double temp_h_value;
long iteration;

for(iteration = 0 ; iteration < *ntime ; iteration ++){
	temp_ref_max_distance = *dist_unit * rayon[iteration];
	for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
		size[id_cluster] = 0;
		}

	for(id_data_line = 0 ; id_data_line < n ; id_data_line ++){
		for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
			temp_current_distance = 0.0;
			for(j = 0 ; j < np ; j ++){
				temp_diff_value = dataset[(np *id_data_line + j)] - weights[(np *id_cluster + j)];
				temp_current_distance = temp_current_distance + (temp_diff_value * temp_diff_value);
				}
			temp_current_distance = sqrt(temp_current_distance);

			if(id_cluster==0 || temp_current_distance < temp_min_distance + EPSILON ){
				temp_min_distance = temp_current_distance;
				id_winner = id_cluster;                            
				}
			}
		winners_of_steps[id_data_line] = id_winner;
		size[id_winner] = size[id_winner] + 1;
		}

	/*
	* browse each cluster and update it if necessary
	*/
	for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
		if(size[id_cluster] != 0){
			// reset extendsize for this cluster
			extendsize = 0;
			
			// set values to 0 for all clusters that have to be modified
			for( j = 0 ; j < np ; j ++ )
				weights[ id_cluster * np + j ] = 0;

			// for each values in the neighbourhood of the current cluster
			for(id_data_line = 0 ; id_data_line < n ; id_data_line ++){
				if (dmcc[winners_of_steps[id_data_line] * nc + id_cluster] < temp_ref_max_distance + EPSILON ){
					if( ( int ) rayon [ iteration ] > 0 && *neightbourhood == 1 ){
						temp_h_value = dmcc [ winners_of_steps [ id_data_line ] * nc + id_cluster ];
						temp_h_value *= temp_h_value;
						temp_h_value = temp_h_value / ( 2 * temp_ref_max_distance * temp_ref_max_distance );
						temp_h_value = exp ( - temp_h_value );
						}
					else
						temp_h_value = 1;

					//somme une ligne du dataset dans le prototype
					for( j = 0 ; j < np ; j ++ ){
						weights[ id_cluster * np + j ] += dataset [ id_data_line * np + j ] * temp_h_value;
						}
					//incremente la ponderation du prototype
					extendsize = extendsize + temp_h_value;
					}
				}
			for ( j = 0 ; j < np ; j ++ ){
				weights [ id_cluster * np + j ] = weights [ id_cluster * np + j ] / extendsize;
				}
			}
		}

	/*
	* initialization
	*/
	for ( id_cluster = 0 ; id_cluster < nc ; id_cluster ++) {
		temp_ss[id_cluster] = 0.0;
		weights_eii[id_cluster] = 0.0;
		}

	/*
	* predict classes with data in winners_of_steps
	*/
	for(id_data_line = 0 ; id_data_line < n ; id_data_line ++){
		for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
			temp_current_distance = 0.0;

			for(j = 0 ; j < np ; j ++){
				temp_diff_value = dataset[(np*id_data_line + j)] - weights[(np*id_cluster + j)];
				temp_current_distance = temp_current_distance + (temp_diff_value * temp_diff_value);
				}

			temp_current_distance = sqrt(temp_current_distance);

			if(id_cluster==0 || temp_current_distance < temp_min_distance + EPSILON ){
				temp_min_distance = temp_current_distance;
				id_winner = id_cluster;                               
				}
			}
		winners_of_steps[id_data_line] = id_winner;
		}

	/*
	* prepare calculation of ss_intra
	*/
	for ( id_data_line = 0 ; id_data_line < n ; id_data_line ++){
		for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
			
			if (dmcc[winners_of_steps[id_data_line] * nc + id_cluster] <= temp_ref_max_distance + EPSILON ){
				
				temp_current_distance = 0.0;
				
				for(j = 0 ; j < np ; j ++){
					temp_diff_value = dataset[(np*id_data_line + j)] - weights[(np*id_cluster + j)];
					temp_current_distance = temp_current_distance + (temp_diff_value * temp_diff_value); 
					}
				
				if((int)rayon[iteration] > 0 && *neightbourhood == 1 ){
					temp_h_value = dmcc[winners_of_steps[id_data_line] * nc + id_cluster];
					temp_h_value *= temp_h_value;
					temp_h_value = temp_h_value/(2 * temp_ref_max_distance * temp_ref_max_distance );
					temp_h_value = exp(-temp_h_value);
					}
				else
					temp_h_value = 1;

				temp_ss[id_cluster] = temp_ss[id_cluster] + temp_h_value * temp_current_distance;
				weights_eii[id_cluster] = weights_eii[id_cluster] + temp_h_value;
				}
			}
		}
	/*
	* calculation of ss_intra
	*/
	for (id_cluster = 0 ; id_cluster < nc ; id_cluster ++){
		if(weights_eii[id_cluster] != 0)
			ss_intra[iteration] = ss_intra[iteration] + temp_ss[id_cluster] / weights_eii[id_cluster];
		}
	ss_intra[iteration] = ss_intra[iteration] / nc;
	
	}

return;
}
