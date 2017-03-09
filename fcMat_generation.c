/*******************************************************************************
 *
 * This routine imports roi data from csv file format and calculates the
 * Pearsons Correlation coefficient for all roi combinations. These values are
 * written to a NxN csv file (where N is the number of rois).
 *
 * Compile with: cc -o fcMat_generation fcMat_generation.c -lm
 * Execute with: ./csv_read filename.csv
 *
 * Written by Thomas Campbell Arnold, Florida State University 2016
 * arnold@psy.fsu.edu
 *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int main(argc,argv)
int argc;
char *argv[];
{
	FILE *f,*fp;
	char *csv,*roi_string,*tr_string;
	char delim;
	
	// check for proper number of commandline arguements 
	if (argc != 4) {
		fprintf(stderr, "\nUsage: %s  <csv file> <number of rois> <number of TRs>\n",argv[0]);
        exit(1);
    }
    
    //set filename equal to csv
    csv = argv[1]; 
	roi_string = argv[2];
	tr_string = argv[3];
	int ROIs = atoi(roi_string);
	int TRs = atoi(tr_string);
	fprintf(stderr, "\nxsum = %d\n",TRs);
	fprintf(stderr, "\nxsum = %d\n",ROIs);
	
	int i, j, val;
	float *data[ROIs];
	for (i=0; i<ROIs; i++){
	data[i] = (float *)malloc(TRs * sizeof(float));
	}
	float temp;
	
	// start timer
	time_t tstart, tend;
	tstart = time(0);
	
	// open csv file
	f = fopen(csv,"r");
	
	// used to index matrix, start at zero
	i=0;
	j=0;
	int n=0;
		
	// read in the data to matrix
	while(val = fscanf(f,"%f%c",&temp,&delim)){
		if(i <= ROIs-1){
			data[i][j] = temp;
			if(i>=ROIs-1){
				i=0;
				j++;
				if(j>=TRs){break;}
			}else{
				i++;
			}
		}else{
			break;
		}
	}
	
	// close csv file to write
	fclose(f);
	
	/***************************************************************************
	 *
	 * This section beneath is for calculating the correlation matrix, also
	 * known as the functional connectivity map (fcMat). Each roi is correlated
	 * with each other roi producing an NxN matrix where N is the number of
	 * rois.
	 *
	***************************************************************************/
	
	int k;
	float x[TRs], y[TRs], xy[TRs], xsquare[TRs], ysquare[TRs];
	float xsum, ysum, xysum, xsqr_sum, ysqr_sum;
    long double num, deno;
	long double coeff;
    
    // open file for writing fcMap to
    fp = fopen("fcMap.csv","w");
    // loop through each combination of i,j from 1 to N
    for(k=0; k<ROIs; k++){
		for(j=0; j<ROIs; j++){
			
			xsum = ysum = xysum = xsqr_sum = ysqr_sum = 0;
			
			// set X and Y
			for(n=0; n<TRs; n++){
				x[n]=data[k][n];
			}
				 
			for(n=0; n<TRs; n++){
				y[n]=data[j][n];
			}
			
			/* find the needed data to manipulate correlation coeff */
			for (i = 0; i < TRs; i++) {
                xy[i] = x[i] * y[i];
                xsquare[i] = x[i] * x[i];
                ysquare[i] = y[i] * y[i];
                xsum += x[i];
                ysum = ysum + y[i];
                xysum = xysum + xy[i];
                xsqr_sum = xsqr_sum + xsquare[i];
                ysqr_sum = ysqr_sum + ysquare[i];
            }
			
            // calculate numerator and denominator 
            // 1 formerly 1.0 for float type
            num = 1 * ((TRs * xysum) - (xsum * ysum));
            deno = 1 * ((TRs * xsqr_sum - xsum * xsum) * (TRs * ysqr_sum - ysum * ysum));
            
            // calcualte correlation coefficient [ num/sqrt(deno) ]
            coeff = num / sqrt(deno);
            
            if (k==1 && j==0){
                	fprintf(stderr, "\nxsum = %f\n",xsum);
                	fprintf(stderr, "\nysum = %f\n",ysum);
                	fprintf(stderr, "\nxysum = %f\n",xysum);
                	fprintf(stderr, "\nxsqr_sum = %f\n",xsqr_sum);
                	fprintf(stderr, "\nysqr_sum = %f\n",ysqr_sum);
                	fprintf(stderr, "\nnum = %Lf\n",num);
                	fprintf(stderr, "\ndeno = %Lf\n",deno);
                	fprintf(stderr, "\ncoeff = %Lf\n",coeff);
            }
                
            // put corr coef in fcMap
            fprintf(fp,"%Lf,",coeff);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	
	// call back time parameter
	tend = time(0);
	fprintf(stderr, "\nIt took %.3f second(s) to run this code.\n",difftime(tend, tstart));
	
    return 0;
}

