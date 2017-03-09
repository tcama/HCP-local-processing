/*********************************************************************
 *
 * This code is designed to carry out several processing steps
 * necessary before the functional connectivity matrices can
 * be calculated.
 * 
 * 1. update_goodvoxels
 * 		The HCP data come with a file called "goodvoxels.nii" which
 * 		is a binary brain mask with additional voxels removed based		
 * 		on local coefficient of variation (COV) (inverse of the
 *		signal-to-noise ratio [SNR]). This file is useful because it 
 *		removes voxels with a low SNR in a localized area, thus 
 * 		avoiding the overpenalization of regions with a lower SNR 
 *		(such as near the sinuses) that would occur when using a 
 *		fixed global. This provides similar voxel removal rates in 
 *		areas of signal dropout as well as in areas of normal BOLD 
 *		signal. HOWEVER, some voxels are still included that have 
 *		negative values for certain TRs. Thus this function also 
 *		searches for remaining negative values and provides an 
 *		updated mask with a "0" indicating the voxels to be excluded.
 * 
 *		Output: goodvoxels_new.csv
 *		Return: goodv_new, pointer to a binary int array (same as output)
 *
 * 2. check_snr
 * 		This function is included for checking the SNR and activity
 * 		of individual voxels. It prints out two files, one with the SNR 
 *      and one with the mean of each voxels timeseries. These two can 
 *		be used to calculate each voxels variance or standard deviation.
 *		It also returns the number of TRs in the data file.
 * 
 * 		Output: mean_data.csv & snr_data.csv
 *		Return: hdr.dim[4], int of the 4th dimension (number of TRs)
 * 
 * 3. roi_extract
 * 		This is the most crucial function in this program. It reads in
 * 		the atlas and fmri data and uses this information to construct
 * 		the ROI timeseries. The end result will be a csv file with
 * 		each column corresponding to a single roi. (i.e. timepoints by
 * 		roi matrix).
 * 
 * 		Output: roi_data.csv
 *		Return: max, int maximum value from roi atlas. This should be the
 *				number of rois. If it isn't, renumber the atlas in MATLAB
 * 
 * 4. meancenter_col
 * 		Perfroms mean centering and matrix transposition. This function
 * 		relys upon the returned values from check_snr and roi_extract to
 * 		setup the data matrix (data[max][TRs]).
 * 
 * 5. voxel_summary
 * 		Prints out voxel removal information to a csv file. Examine in
 * 		Excel.
 * 
 * To compile:
 * You need to put a copy of the nifti1.h header file in this directory
 * before compiling the code. It can be obtained from the NIFTI homepage
 * http://nifti.nimh.nih.gov/ or from the niftilib SourceForge site 
 * http://niftilib.sourceforge.net/
 * 
 * complilation command:
 * cc roi_preprocess -o roi_preprocess.cpp
 * 
 * Run:
 * ./roi_preprocess rfmri.nii atlas.nii goodvoxels.nii
 *
 * where rfmri.nii is the resting-state scan from HCP, atlas.nii is a MNI
 * normalized atlas of your ROIs (coded as 1 through R, with R equal to
 * the total number of regions), and goodvoxels.nii is the goodvoxels
 * image obtained from the HCP website. This last part can be removed if
 * desirable.
 *
 * This code was built using blocks of code obtained from Kate Fissel,
 * University of Pittsburgh. Her original comments are preserved at the
 * bottom of this file.
 *
 * Written by Thomas Campbell Arnold, Florida State University 2016
 * arnold@psy.fsu.edu
 *
 *********************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "nifti1.h"

typedef float MY_DATATYPE;

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352

using namespace std;

// PROTOTYPES
int * update_goodvoxels(char *data_file, char *goodv_file);
int check_snr(char *data_file);
int roi_extract(char *data_file, char *atlas_file, int *goodv_update);
int meancenter_col(int max, int trs);
int voxel_summary(int max, char* data_file, char *atlas_file, char *goodv_file, int *goodv_update);

int main(int argc, char *argv[]) 
{
	
	char *data_file, *atlas_file, *goodv_file;
	
	// check for proper number of commandline arguements 
	if (argc != 4) {
		fprintf(stderr, "\nUsage: %s  <header file> <data file> <atlas.nii>\n",argv[0]);
        exit(1);
    }
	
	// get filenames from arguement
	data_file = argv[1];
	atlas_file = argv[2];
	goodv_file = argv[3];

	// set a time parameter
	time_t tstart, tend;
	tstart = time(0);
	
	// call functions
	int *goodv_update = update_goodvoxels(data_file, goodv_file);
	int trs = check_snr(data_file);
	fprintf(stderr,"\nnumber of TRs = %d\n",trs);
	int max = roi_extract(data_file, atlas_file, goodv_update);
	fprintf(stderr,"\nNumber of Brain Regions = %d\n",max);
	meancenter_col(max,trs);
	fprintf(stderr,"\nfinished mean centering timeseries\n");
	voxel_summary(max, data_file, atlas_file, goodv_file, goodv_update);
	fprintf(stderr,"\nThe code has finished running. Check output.\n");

	// call back time parameter
	tend = time(0);
	fprintf(stderr, "\nIt took %.3f second(s) to run this code.\n",difftime(tend, tstart));
	
	exit(0);
}

/***************************** FIRST FUNCTION *****************************/

int * update_goodvoxels(char *data_file, char *goodv_file)
{
	nifti_1_header hdr;
	nifti_1_header goodv_hdr;
	FILE *f;
	int ret, i, t, neg_data, neg_goodv;
	MY_DATATYPE *data=NULL, *goodv_data=NULL;
	
	// initialize negative counts at 0
	neg_data = 0;
	neg_goodv = 0;
	
	// FIRST LOAD IN GOODVOXELS.NII INFORMATION
	
	// open hdr file
	f = fopen(goodv_file,"r");
	if (f == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",goodv_file);
        exit(1);
    }

    // read hdr information
    ret = fread(&goodv_hdr, MIN_HEADER_SIZE, 1, f);
    if (ret != 1) {
        fprintf(stderr, "\nError reading header file %s\n",goodv_file);
        exit(1);
    }
    fclose(f);
    
    // print hdr information to screen
    fprintf(stderr, "\n%s header information:",goodv_file);
    fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",goodv_hdr.dim[1],goodv_hdr.dim[2],goodv_hdr.dim[3],goodv_hdr.dim[4]);
    fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",goodv_hdr.datatype,goodv_hdr.bitpix);
    fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",goodv_hdr.scl_slope,goodv_hdr.scl_inter);
    fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(goodv_hdr.vox_offset));
    fprintf(stderr, "\n");
    
    // open data file jump to vox_offset before beginning to read
    f = fopen(goodv_file,"r");
    if (f == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",goodv_file);
        exit(1);
	}
	ret = fseek(f, (long)(goodv_hdr.vox_offset), SEEK_SET);
	
	// allocate space to hold the data
	goodv_data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]);
	if (goodv_data == NULL) {
		fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
		exit(1);
	}
	
	ret = fread(goodv_data, sizeof(MY_DATATYPE), goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3], f);
	if (ret != goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]) {
		fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",goodv_data,ret);
		exit(1);
	}
	
	// close data file
	fclose(f);
	
	// BEGIN LOADING RFMRI DATA
	
	// open hdr file
	f = fopen(data_file,"r");
	if (f == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",data_file);
        exit(1);
    }

    // read hdr information
    ret = fread(&hdr, MIN_HEADER_SIZE, 1, f);
    if (ret != 1) {
        fprintf(stderr, "\nError reading header file %s\n",data_file);
        exit(1);
    }
    fclose(f);
    
    // print hdr information to screen
    fprintf(stderr, "\n%s header information:",data_file);
    fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
    fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
    fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
    fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
    fprintf(stderr, "\n");
	
	// open data file jump to vox_offset before beginning to read
    f = fopen(data_file,"r");
    if (f == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",data_file);
        exit(1);
	}
	ret = fseek(f, (long)(hdr.vox_offset), SEEK_SET);
	
	// loop through each of the TRs (1:hdr.dim[4])
	for(t=1; t<=hdr.dim[4]; t++){
		
		// allocate space to hold the data
		data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
		if (data == NULL) {
			fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
			exit(1);
		}
		
		ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], f);
		if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
			fprintf(stderr, "\nError reading volume %d from %s (%d)\n",t,data_file,ret);
			exit(1);
		}
		
		// add data from new TR to the sum
		for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
		{
			// if data value is below zero add to count
			if (data[i] < 0){
				neg_data++;
			}
			
			// change goodv_data if it is a basal voxel
			if (goodv_data[i] == 0){
				//REMOVED BASAL AFFECT if (basal_data[i] == 1){
					goodv_data[i] = 0;
				//}
			}
			
			// multiply data by goodvoxels
			data[i] = data[i] * goodv_data[i];
			
			// if data value is below zero add to count (after goodvoxels)
			if (data[i] < 0){
				neg_goodv++;
				goodv_data[i] = 0;
			}
		}
		
	}
	
	// close data file
	fclose(f);

	// print results
	fprintf(stderr, "\nTotal instances of negative TRs in raw data = %d\n",neg_data);
	fprintf(stderr, "\nInstances of negative TRs after goodvoxels correction = %d\n",neg_goodv);
	
	// open file to print out new goodvoxels if needed
	if (neg_goodv > 0){
		f = fopen("goodvoxels_new.csv","w");
		for (i=0; i<goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]; i++){
				fprintf(f,"%d,",goodv_data[i]);
		}
		fclose(f);
	}
	
	// convert to int
	static int goodv_new[1000000]; // This hard coding should be changed to deal with dynamically sized files
	//static int goodv_new[goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]]; 
	
	memset(goodv_new,0,goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]*sizeof(int));
	for (i=0; i<goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]; i++)
	{
		goodv_new[i] = goodv_data[i];
	}
	

	delete[] data;
	delete[] goodv_data;
	return goodv_new;
};
 
/***************************** SECOND FUNCTION *****************************/
 
int check_snr(char *data_file)
{
	nifti_1_header hdr;
	FILE *f, *f2;
	int ret, i, t;
	MY_DATATYPE *data=NULL, *mean_data=NULL, *var_data=NULL;
	
	// open hdr file
	f = fopen(data_file,"r");
	if (f == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",data_file);
        exit(1);
    }

    // read hdr information
    ret = fread(&hdr, MIN_HEADER_SIZE, 1, f);
    if (ret != 1) {
        fprintf(stderr, "\nError reading header file %s\n",data_file);
        exit(1);
    }
    fclose(f);
    
    // open data file jump to vox_offset before beginning to read
    f = fopen(data_file,"r");
    if (f == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",data_file);
        exit(1);
	}
	ret = fseek(f, (long)(hdr.vox_offset), SEEK_SET);
	
	// allocate memory for holding mean
	mean_data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
	memset(mean_data,0,hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*sizeof(MY_DATATYPE));

	// loop through each of the TRs (1:hdr.dim[4])
	for(t=1; t<=hdr.dim[4]; t++){

		// allocate space to hold the data
		data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
		if (data == NULL) {
			fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
			exit(1);
		}
		
		ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], f);
		if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
			fprintf(stderr, "\nError reading volume %d from %s (%d)\n",t,data_file,ret);
			exit(1);
		}
		
		// add data from new TR to the sum
		for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
		{
			mean_data[i] += data[i];
		}
		
	}

	// divide the sum of all TR intensities by total number of TTs
	for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
	{
		mean_data[i] /= hdr.dim[4];
	}
	

	// close data file
	fclose(f);
	
	// BEGIN CALCULATION OF VARIANCE
	
	// open data file jump to vox_offset before beginning to read
    f = fopen(data_file,"r");
    if (f == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",data_file);
        exit(1);
	}
	ret = fseek(f, (long)(hdr.vox_offset), SEEK_SET);

	// allocate memory for holding mean
	var_data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
	memset(var_data,0,hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*sizeof(MY_DATATYPE));

	// loop through each of the TRs (1:hdr.dim[4])
	for(t=1; t<=hdr.dim[4]; t++){
		
		// allocate space to hold the data
		data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
		if (data == NULL) {
			fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
			exit(1);
		}
		
		ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], f);
		if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
			fprintf(stderr, "\nError reading volume %d from %s (%d)\n",t,data_file,ret);
			exit(1);
		}
		
		// add data from new TR to the sum
		for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
		{
			var_data[i] += (data[i] - mean_data[i])*(data[i] - mean_data[i]);
		}
		
	}
	
	// divide the sum of all TR intensities by total number of TTs
	for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
	{
		var_data[i] /= hdr.dim[4];
	}
	
	// close data file
	fclose(f);

	// variables to count voxels w/ mean > 1 and calculate mean SNR
	int n = 0;
	float tot = 0;
	
	// open file for writing
	f = fopen("mean_data.csv", "w");
	f2 = fopen("snr_data.csv", "w");
	
	// calculate SNR: divide mean by sqrt of variance (std)
	for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
	{
		if (mean_data[i] > 1){
			fprintf(f,"%f,",mean_data[i]); // print mean before calculating snr
			mean_data[i] /= sqrt(var_data[i]);
			fprintf(f2,"%f,",mean_data[i]); // print snr (note: mean_data is altered)
			n++;
			tot += mean_data[i];
		}else{
			mean_data[i] = 0;
			fprintf(f,"%f,",mean_data[i]); // NOTE: This does not print actual mean or SNR, it prints 0.
			fprintf(f2,"%f,",mean_data[i]); //      Mean/SNR could be in below 1 or even negative.
		}
		
	}
	
	// close output files
	fclose(f);
	fclose(f2);
	
	// print out voxel count and mean SNR to screen
	tot /= n;
	fprintf(stderr, "\nTotal voxels before goodvoxels = %d\n\nmean SNR = %f\n",n,tot);
	

	delete[] data;
	delete[] mean_data;
	delete[] var_data;
	return hdr.dim[4];
};
 
/***************************** THIRD FUNCTION *****************************/
 
int roi_extract(char *data_file, char *atlas_file, int *goodv_update){
	nifti_1_header hdr;
	nifti_1_header aal_hdr;
	FILE *fp,*f,*aal;
	int ret,i,t,max;
	MY_DATATYPE *data=NULL, *mean_intensity=NULL;
	double *roi_data=NULL;
	int *roi_array=NULL;
	
	/**********************************************************************
	 *
	 * Process Atlas
	 *
	**********************************************************************/
	
	// open and read header
	aal = fopen(atlas_file,"r");
	if (aal == NULL) {
	        fprintf(stderr, "\nError opening header file %s\n",atlas_file);
	        exit(1);
	}
	
	ret = fread(&aal_hdr, MIN_HEADER_SIZE, 1, aal);
	if (ret != 1) {
	        fprintf(stderr, "\nError reading header file %s\n",atlas_file);
	        exit(1);
	}
	fclose(aal);
	
	// print a little header information
	fprintf(stderr, "\n%s header information:",atlas_file);
	fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",aal_hdr.dim[1],aal_hdr.dim[2],aal_hdr.dim[3],aal_hdr.dim[4]);
	fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",aal_hdr.datatype,aal_hdr.bitpix);
	fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",aal_hdr.scl_slope,aal_hdr.scl_inter);
	fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(aal_hdr.vox_offset));
	fprintf(stderr, "\n");
	
	// check if there are multiple timepoints in atlas
	if (aal_hdr.dim[4] >= 2){
			fprintf(stderr, "\nError - multiple timepoints in atlas file: %s\n",atlas_file);
	        exit(1);
	}
	
	// open the datafile, jump to data offset
	aal = fopen(atlas_file,"r");
	if (aal == NULL) {
	        fprintf(stderr, "\nError opening data file %s\n",atlas_file);
	        exit(1);
	}
	ret = fseek(aal, (long)(aal_hdr.vox_offset), SEEK_SET);
	
	// allocate buffer and read first 3D volume from data file
	roi_data = (double *) malloc(sizeof(double) * aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]);
	if (roi_data == NULL) {
	        fprintf(stderr, "\nError allocating data buffer for %s\n",atlas_file);
	        exit(1);
	}
	ret = fread(roi_data, sizeof(double), aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3], aal);
	if (ret != aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]) {
	        fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",atlas_file,ret);
	        exit(1);
	}
	
	// the scaling section is commented out as it is unnecessary, should it be
	// needed it can simply be added back in.
	/********** scale the data buffer 
	if (aal_hdr.scl_slope != 0) {
	        for (i=0; i<aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]; i++)
	                roi_data[i] = (roi_data[i] * aal_hdr.scl_slope) + aal_hdr.scl_inter;
	}*/
	
	// get max value of roi
	max = roi_data[0]; // declare max is first value in data array
	for (i=0; i<aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]; i++)
	{
	 if(max < roi_data[i]){
	 	max = roi_data[i];
	 }
	}
	
	int n = 0;
	// convert from float to int
	// multiply by goodvoxels to remove voxels from atlas
	roi_array = (int *) malloc(sizeof(int) * aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]);
	memset(roi_array,0,aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]*sizeof(int));


	for (i=0; i<aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]; i++)
	{
		roi_array[i] = roi_data[i] * goodv_update[i];
		if (goodv_update[i] == 1){
			n++;
		}
	}
	
	fprintf(stderr, "\nTotal voxels after goodvoxels update = %d\n",n);
	
	fclose(aal);
	
	/**********************************************************************
	 *
	 * Process rfMRI data
	 *
	**********************************************************************/
	
	// open and read header
	fp = fopen(data_file,"r");
	if (fp == NULL) {
	        fprintf(stderr, "\nError opening header file %s\n",data_file);
	        exit(1);
	}
	
	ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
	if (ret != 1) {
	        fprintf(stderr, "\nError reading header file %s\n",data_file);
	        exit(1);
	}
	fclose(fp);
	
	mean_intensity = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
	memset(mean_intensity,0,hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*sizeof(MY_DATATYPE));

	f = fopen("roi_data.csv", "w"); // open CSV file
	fp = fopen(data_file,"r");
	for(t=1; t<=hdr.dim[4]; t++)
	{
	// open the datafile, jump to data offset
	if (fp == NULL) {
	        fprintf(stderr, "\nError opening data file %s\n",data_file);
	        exit(1);
	}
	if (t == 1){
		ret = fseek(fp, NII_HEADER_SIZE, SEEK_SET);
	}
	
	// allocate buffer and read first 3D volume from data file
	data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
	if (data == NULL) {
	        fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
	        exit(1);
	}
	ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
	if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
	        fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",data_file,ret);
	        exit(1);
	}
	
	// the scaling section is commented out as it is unnecessary, should it be
	// needed it can simply be added back in.
	/********** scale the data buffer  
	if (hdr.scl_slope != 0) {
	        for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
	                data[i] = (data[i] * hdr.scl_slope) + hdr.scl_inter;
	}*/

	 // roi count 
	 int roi_count[max+1];
	 memset(roi_count,0,(max+1)*sizeof(int));


	// print mean roi timeseries of data
	float total[max+1];
	memset(total,0,(max+1)*sizeof(float));
	for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++){
		if(data[i] >= 1){
	        total[roi_array[i]] = total[roi_array[i]] + data[i];
	        roi_count[roi_array[i]]++;
	    }
	}


	
	// divide by number of voxels included in roi
	for (i=1; i<max+1; i++){
	        total[i] /= roi_count[i];
	}
	
	
	// write CSV file
	
	for (i=1; i<max+1; i++)
	{
		if (i<max){
			fprintf(f,"%f,",total[i]);
		}else{
			fprintf(f,"%f",total[i]);
		}
	}
	
	fprintf(f,"\n");
	
	}
	
	 fclose(fp);
	 fclose(f);
	 
	delete[] data;
	delete[] mean_intensity;
	delete[] roi_data;
	delete[] roi_array;
	return max;

};
 
/***************************** FOURTH FUNCTION *****************************/
 
int meancenter_col(int max, int trs){

	FILE *f,*fp;
	char delim;
	int i, j, n, val;
	float temp;

	// allocate memory for data
	float *data[max];
	for (i=0; i<max; i++){
	data[i] = (float *)malloc(trs * sizeof(float));
	}

	// open csv file
	f = fopen("roi_data.csv","r");

	// used to index matrix, start at zero
	i=0;
	j=0;
	n=0;
	
	// read in data and put into matrix
	while(val = fscanf(f,"%f,",&temp)){
		if(i <= max-1){
			// is int value, put in array
			data[i][j] = temp;
			//fprintf(stderr,"\ni=%i j=%i\n",i,j); // used to check proper i,j output
			if(i>=max-1){
				i=0;
				j++;
				if(j>=trs){break;}
			}else{
				i++;
			}
		}else{
			break;
		}
	}

	// remove the mean
	float tp_sum;
	for(i=0; i<max; i++){
		tp_sum = 0;
		for(j=0; j<trs; j++){
			tp_sum += data[i][j];
		}
		tp_sum /= trs;
		for(j=0; j<trs; j++){
			data[i][j] -= tp_sum;
		}
	}
	
	// open csv file to write
	fp = fopen("mc.csv","w");			
	for(j=0; j<trs; j++){
		for(i=0; i<max; i++){
			if (i<max-1){
				fprintf(fp,"%f,",data[i][j]);
			}else{
				fprintf(fp,"%f",data[i][j]);
			}
			
		}
		fprintf(fp,"\n",data[i][j]);
	}
	
	
	// close csv file to write
	fclose(f);
	fclose(fp);
	
    return 0;
};


int voxel_summary(int max, char* data_file, char *atlas_file, char *goodv_file, int *goodv_update){
	
	nifti_1_header goodv_hdr;
	nifti_1_header aal_hdr;
	nifti_1_header hdr;
	FILE *f,*aal;
	int ret,i,j;
	MY_DATATYPE *data=NULL, *goodv_data=NULL;
	double *roi_data=NULL;
	int *goodv_array=NULL, *roi_array=NULL;
	
	// load in atlas information
	//open and read header
	aal = fopen(atlas_file,"r");
	if (aal == NULL) {
	        fprintf(stderr, "\nError opening header file %s\n",atlas_file);
	        exit(1);
	}
	
	ret = fread(&aal_hdr, MIN_HEADER_SIZE, 1, aal);
	if (ret != 1) {
	        fprintf(stderr, "\nError reading header file %s\n",atlas_file);
	        exit(1);
	}
	fclose(aal);
	
	// open the datafile, jump to data offset
	aal = fopen(atlas_file,"r");
	if (aal == NULL) {
	        fprintf(stderr, "\nError opening data file %s\n",atlas_file);
	        exit(1);
	}
	ret = fseek(aal, (long)(aal_hdr.vox_offset), SEEK_SET);
	
	// allocate buffer and read first 3D volume from data file
	roi_data = (double *) malloc(sizeof(double) * aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]);
	if (roi_data == NULL) {
	        fprintf(stderr, "\nError allocating data buffer for %s\n",atlas_file);
	        exit(1);
	}
	ret = fread(roi_data, sizeof(double), aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3], aal);
	if (ret != aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]) {
	        fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",atlas_file,ret);
	        exit(1);
	}
	
	// close data file
	fclose(aal);
	
	// convert from float to int
	roi_array = (int *) malloc(sizeof(int) * aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]);
	memset(roi_array,0,aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]*sizeof(int));

	for (i=0; i<aal_hdr.dim[1]*aal_hdr.dim[2]*aal_hdr.dim[3]; i++)
	{
		roi_array[i] = roi_data[i];
	}
	
	// load in goodvoxels.nii information
	
	// open hdr file
	f = fopen(goodv_file,"r");
	if (f == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",goodv_file);
        exit(1);
    }

    // read hdr information
    ret = fread(&goodv_hdr, MIN_HEADER_SIZE, 1, f);
    if (ret != 1) {
        fprintf(stderr, "\nError reading header file %s\n",goodv_file);
        exit(1);
    }
    fclose(f);
    
    // open data file jump to vox_offset before beginning to read
    f = fopen(goodv_file,"r");
    if (f == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",goodv_file);
        exit(1);
	}
	ret = fseek(f, (long)(goodv_hdr.vox_offset), SEEK_SET);
	
	// allocate space to hold the data
	goodv_data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]);
	if (goodv_data == NULL) {
		fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
		exit(1);
	}
	
	ret = fread(goodv_data, sizeof(MY_DATATYPE), goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3], f);
	if (ret != goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]) {
		fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",goodv_data,ret);
		exit(1);
	}
	
	// close data file
	fclose(f);
	
	// convert from float to int
	goodv_array = (int *) malloc(sizeof(int) * goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]);
	memset(goodv_array,0,goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]*sizeof(int));

	for (i=0; i<goodv_hdr.dim[1]*goodv_hdr.dim[2]*goodv_hdr.dim[3]; i++)
	{
		goodv_array[i] = goodv_data[i];
	}
	
	// load in data information, just first slice

	// open hdr file
	f = fopen(data_file,"r");
	if (f == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",data_file);
        exit(1);
    }

    // read hdr information
    ret = fread(&hdr, MIN_HEADER_SIZE, 1, f);
    if (ret != 1) {
        fprintf(stderr, "\nError reading header file %s\n",data_file);
        exit(1);
    }
    fclose(f);
    
    // open data file jump to vox_offset before beginning to read
    f = fopen(data_file,"r");
    if (f == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",data_file);
        exit(1);
	}
	ret = fseek(f, (long)(hdr.vox_offset), SEEK_SET);
	
	// allocate space to hold the data
	data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
	if (data == NULL) {
		fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
		exit(1);
	}
	
	ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], f);
	if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
		fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",data,ret);
		exit(1);
	}
	
	// close data file
	fclose(f);
	
	// open file for writing
	f = fopen("voxels_summary.csv","w");
	
	// count voxels for each condition
	int roi_count[max+1];
	memset(roi_count,0,(max+1)*sizeof(int));
	int goodv_count[max+1];
	memset(goodv_count,0,(max+1)*sizeof(int));
	int goodv_new_count[max+1];
	memset(goodv_new_count,0,(max+1)*sizeof(int));
	
	for (j=0; j<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; j++)
	{
		
		if(data[j] >= 1)
		{
			roi_count[roi_array[j]]++;
			goodv_count[roi_array[j]*goodv_array[j]]++;
			goodv_new_count[roi_array[j]*goodv_update[j]]++;
		};
		
	};
	
	for (i=0; i<max+1; i++)
	{
		fprintf(f,"%i,%i,%i,%i,%f,%f\n",i,roi_count[i],goodv_count[i],goodv_new_count[i],((goodv_count[i]*1.0)/(roi_count[i]*1.0))*100,((goodv_new_count[i]*1.0)/(roi_count[i]*1.0))*100);
	}
	
	// close output file
	fclose(f);
	
	return 0;
};

/*********************************************************************
 *
 * Helpful comments from parent package. Note the original code was designed
 * to be compiled in c, not c++. The links near the bottom can supply you
 * with the original code, which is brief. 
 *
 * Very simple code snippets to read/write nifti1 files
 * This code is placed in the public domain.
 *
 * If you are the type who doesn't want to use a file format unless
 * you can write your own i/o code in less than 30minutes, this
 * example is for you.
 *
 * This code does not deal with wrong-endian data, compressed data,
 * the new qform/sform orientation codes, parsing filenames, volume-
 * wise or timecourse-wise data access or any of a million other very useful
 * things that are in the niftilib i/o reference libraries.
 * We encourage people to use the niftilib reference library and send
 * feedback/suggestions, see http://niftilib.sourceforge.net/
 * But, if that is too much to tackle and you just want to jump in, this
 * code is a starting point.
 * This code was written for maximum readability, not for the greatest
 * coding style.
 *
 *
 * If you are already a little familiar with reading/writing Analyze
 * files of some flavor, and maybe even have some of your own code, here
 * are the most important things to be aware of in transitioning to nifti1:
 *
 * 1. nii vs .hdr/.img
 *      nifti1 datasets can be stored either in .hdr/.img pairs of files
 *      or in 1 .nii file.  In a .nii file the data will start at the byte
 *      specified by the vox_offset field, which will be 352 if no extensions
 *      have been added.  And, nifti1 really does like that magic field set
 *      to "n+1" for .nii and "ni1" for .img/.hdr
 *
 * 2. scaling
 *      nifti1 datasets can contain a scaling factor.  You need to check the
 *      scl_slope field and if that isn't 0, scale your data by 
 *      Y * scl_slope  + scl_inter
 *
 * 3. extensions
 *      nifti1 datasets can have some "extension data" stuffed after the 
 *      regular header.  You can just ignore it, but, be aware that a
 *      .hdr file may be longer than 348 bytes, and, in a .nii file
 *      you can't just jump to byte 352, you need to use the vox_offset
 *      field to get the start of the image data.
 *
 * 4. new datatypes
 *      nifti1 added a few new datatypes that were not in the Analyze 7.5
 *      format from which nifti1 is derived.  If you're just working with
 *      your own data this is not an issue but if you get a foreign nifti1
 *      file, be aware of exotic datatypes like DT_COMPLEX256 and mundane
 *      things like DT_UINT16.
 *
 * 5. other stuff
 *     nifti1 really does like the dim[0] field set to the number of
 *     dimensions of the dataset.  Other Analyze flavors might not
 *     have been so scrupulous about that.
 *     nifti1 has a bunch of other new fields such as intent codes,
 *     qform/sform, etc. but, if you just want to get your hands on
 *     the data blob you can ignore these.  Example use of these fields
 *     is in the niftilib reference libraries.
 *
 *
 *
 * To compile:
 * You need to put a copy of the nifti1.h header file in this directory.
 * It can be obtained from the NIFTI homepage  http://nifti.nimh.nih.gov/
 * or from the niftilib SourceForge site http://niftilib.sourceforge.net/
 * 
 * cc -o nifti1_read_write nifti1_read_write.c
 * 
 * 
 * To run:
 * nifti1_read_write -w abc.nii abc.nii
 * nifti1_read_write -r abc.nii abc.nii
 * 
 * 
 * The read method is hardcoded to read float32 data.  To change
 * to your datatype, just change the line:
 * typedef float MY_DATATYPE;
 *
 * The write method is hardcoded to write float32 data.  To change
 * to your datatype, change the line:
 * typedef float MY_DATATYPE;
 * and change the lines:
 * hdr.datatype = NIFTI_TYPE_FLOAT32;
 * hdr.bitpix = 32;
 *
 *
 * Written by Kate Fissell, University of Pittsburgh, May 2005.
 *
 *********************************************************************/