#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

struct DATA
{
	int count;
	char SAMPLE[20];
	char TARGET[20];
	double CT_VAL;
	double DCT;
	double DDCT;
	double FOLD_CHANGE;
};

int getN_lines(char *FILENAME)
{
	FILE *pCSV; pCSV = fopen(FILENAME,"r");
	int count=0;

	if(pCSV == NULL)
	{
		printf("ERROR: Could not open file.\n");
		exit(1);
	}
	for (char c = getc(pCSV); c!=EOF; c=getc(pCSV))
	{
		if (c == '\n')
			++count;
	}

	return count;

}

int * GET_COL_LOCATION(char *FILENAME)
{
	FILE * pCSV; 
	pCSV = fopen(FILENAME,"r");
	
	if(pCSV==NULL)
	{
		printf("ERROR: No file found!\n");
		exit(1);
	}

	char line[1024]; char * pline;
	fgets(line,1024,pCSV);
	int sample_N=0, target_N=0, CT_n=0, c=0;
	pline = strtok(line,",");
	while(pline!=NULL)
	{
		if(strncmp(pline,"Sample Name",11) == 0)
			sample_N = c;
		if(strncmp(pline,"Target Name",11) == 0)
			target_N = c;
		if(strncmp(pline,"CT",2) == 0)
			CT_n = c;
		++c;
		pline = strtok(NULL,",");
	}
	int check = sample_N + target_N + CT_n;
	if(check < 3)
	{
		printf("ERROR: First line of CSV file does not contain Sample/Target Names or CT values\n");
		exit(1);
	}
	
	fclose(pCSV);
	int vals[3] = {sample_N, target_N, CT_n};
	int * locations = vals;
	return locations;
}

bool IS_IN(int len, char LIST[len][40], char * VALUE)
{
	int t1, t2, j=len;
	bool found = false;
	for (int i=0; i<len;i++)
	{
		t1 = strncmp(VALUE,LIST[i],40);
		t2 = strncmp(VALUE,LIST[j],40); j--;
		if(t1==0 || t2==0)
		{
			found = true;
			break;
		}
	}
	return found;
}

double get_DCT(int i, int control_location,int struct_len, struct DATA DATA[struct_len])
{
	if(strcmp(DATA[i].SAMPLE,DATA[control_location].SAMPLE)==0)  //if the sample is the control
	{
		if(strcmp(DATA[i].TARGET, DATA[control_location].TARGET)==0) // if the tarfet is control
		{
			DATA[i].DCT = 0; 
		}
		else
		{
			DATA[i].DCT = DATA[i].CT_VAL - DATA[control_location].CT_VAL; //if sample is control, but target isnt control, subtract the value from gene control
		}
	}
	else //if the sample is not the control. find where its corresponding control is and subtract
	{
		for(int j=0;j<struct_len;j++) 
		{
			if(strcmp(DATA[i].SAMPLE, DATA[j].SAMPLE)==0 && strcmp(DATA[j].TARGET, DATA[control_location].TARGET)==0)
				DATA[i].DCT = DATA[i].CT_VAL - DATA[j].CT_VAL;
		}
	}
	return DATA[i].DCT;
}

double get_DDCT(int i, int control_location, int struct_len, struct DATA DATA[struct_len])
{
	if(strcmp(DATA[i].SAMPLE,DATA[control_location].SAMPLE)==0) //if sample is the control, the difference from control is 0
	{
		DATA[i].DDCT = 0;
	}
	else if(strcmp(DATA[i].TARGET,DATA[control_location].TARGET)==0) //if the gene name is the control gene.. It shouldnt be calculated and it's made to = 0
	{
		DATA[i].DDCT = 0;
	}
	else //if its not the control, find the value of the target gene in the control samples and subtract it from the current value
	{
		for(int j=0;j<struct_len;j++)
		{
			if(strcmp(DATA[j].TARGET,DATA[i].TARGET)==0 && strcmp(DATA[j].SAMPLE,DATA[control_location].SAMPLE)==0)
			{
				DATA[i].DDCT = DATA[i].DCT - DATA[j].DCT;
				break;
			}

		}
	}
	return DATA[i].DDCT;
}

int main(int argc, char *argv[])
{
	/*Argument handling*/
	if (argv[1] == NULL)
	{
		puts("ERROR: No arguments supplied!");
		puts("Usage: ./qPCR <filename.csv> or use ./qPCR --help");
		exit(1);
	}
	else if (strncmp(argv[1],"--help",6) == 0 || strncmp(argv[1],"-h",2) == 0)
	{
		puts("Usage: ./qPCR <filename.csv>");
		exit(0);
	}
	else
	{
		printf("Opening %s...\n\n",argv[1]);
	}

	/*Getting filename from argument and opening file*/	
	char * FILENAME = argv[1];
	int N_LINES = getN_lines(FILENAME);
	int * COL_LOCATION = GET_COL_LOCATION(FILENAME);  //this determines where the rows are that contain the data I'm interested in calculating.
	int Ns=COL_LOCATION[0], Nt=COL_LOCATION[1], nCT=COL_LOCATION[2]; //indicates the colouln numbers for Samplenames, Targetnames, and thier CT values
	
	/*opening file and parsing file to extract data*/
	FILE * pCSV; pCSV = fopen(FILENAME,"r"); 
	char TARGET_NAMES[N_LINES][20], SAMPLE_NAMES[N_LINES][20], line[1024];
	double CT_VALS[N_LINES];
	line[1024];
	
	for(int l=0;l < N_LINES;l++)
	{
		fgets(line,sizeof(line),pCSV);
		int i=0, N_samples = 0, N_targets=0; 
		char * pline = strtok(line,",");
		
		while(pline!=NULL)
		{
			if(i == Ns)
				strncpy(SAMPLE_NAMES[l],pline,19);
			else if(i == Nt)
				strncpy(TARGET_NAMES[l],pline,19);
			else if(i == nCT)
				CT_VALS[l] = strtod(pline,&pline);
			i++; pline = strtok(NULL,",");
		}
	}

	/*Making an array of all the different conditions found in our dataset*/
	char CONDITIONS[N_LINES][40];
	for(int i=0; i<N_LINES; i++)
	{
		strncpy(CONDITIONS[i],SAMPLE_NAMES[i],19);
		strncat(CONDITIONS[i],TARGET_NAMES[i],19);
	}

	/*Start making our struct containing the data we want to extract*/
	int n_cond = 0;
	struct DATA DATA_ARR[N_LINES-1];

	/*lines 181-222 gets rid od duplicates, calulcates replicates per sample, and gets the mean of value of each found condition*/
	for(int i=0; i < N_LINES; i++)
	{
		if(i == 0)
			printf("Sample# %-14s%-14sCT\n",SAMPLE_NAMES[i],TARGET_NAMES[i]);
		else if(CONDITIONS[i] != NULL)
		{
			bool REPLICATE_EXISTS; 
			REPLICATE_EXISTS = IS_IN(i-1,CONDITIONS,CONDITIONS[i]);
			if (!REPLICATE_EXISTS)
			{
				printf("%-8d%-14s%-14s%-8lf",i,SAMPLE_NAMES[i],TARGET_NAMES[i],CT_VALS[i]);
				puts(" ---> New condition added"); 
				DATA_ARR[n_cond].count=1;
				DATA_ARR[n_cond].CT_VAL = CT_VALS[i];
				strncpy(DATA_ARR[n_cond].SAMPLE,SAMPLE_NAMES[i],19);
				strncpy(DATA_ARR[n_cond].TARGET,TARGET_NAMES[i],19);
				n_cond++;
			}
			else
			{
				printf("%-8d%-14s%-14s%-8lf",i,SAMPLE_NAMES[i],TARGET_NAMES[i],CT_VALS[i]);
				puts(" ---> Adding value to previous replicate"); 
				int j = 0;
				while(j<n_cond)
				{	
					char CUR_COND[40];
					strncpy(CUR_COND,DATA_ARR[j].SAMPLE,19);
					strncat(CUR_COND,DATA_ARR[j].TARGET,19);
					if(strcmp(CUR_COND,CONDITIONS[i]) == 0)
					{
						DATA_ARR[j].count++;
						double AVE_CT = (DATA_ARR[j].CT_VAL * DATA_ARR[j].count) + CT_VALS[i];
						AVE_CT = AVE_CT/(DATA_ARR[j].count+1);
						DATA_ARR[j].CT_VAL = AVE_CT;
						break;
					}
					j++;
				}
			}
		}
	}
	/*printing for user feedback*/
	puts("\nAveraged data:\n");
	puts("#   SAMPLE_NAME    TARGET_NAME #Reps   AVE_CT"); 
	for(int i=0; i<n_cond; i++)
	{
		printf("%-8d%-14s%-10s%-5d%5lf\n",i,DATA_ARR[i].SAMPLE,DATA_ARR[i].TARGET,DATA_ARR[i].count,DATA_ARR[i].CT_VAL);
	}
	/*Gets the control conditions to normalize data to. First assumes it's data from first column*/
	bool is_first_control; int control_location=0;
	printf("\nAre %s and %s the control conditions? (y/n) >> ",DATA_ARR[control_location].SAMPLE,DATA_ARR[control_location].TARGET);
	char a[3]; fgets(a, sizeof(a), stdin);
	if(strncmp(a,"y",1) != 0)
	{
		printf("\nEnter the number corresponsding to a control condition >> ");
		scanf("%d",&control_location);
		if(control_location > n_cond-1)
		{
			printf("ERROR: Not in range!\n");
			exit(1);
		}
		printf("You entered %s and %s as control conditions\n",DATA_ARR[control_location].SAMPLE,DATA_ARR[control_location].TARGET);
	}
	/*Calulating the difference between each target from the control target gene for each sample to get deltaCT
	Then calculating the difference of each samples value to the values determined by the control gene lastly determining fond change from controls*/
	
	puts("\n#   SAMPLE_NAME    TARGET_NAME #Reps   \tAVE_CT\t\tDCT\t\tDDCT\t\tFOLD_CHANGE"); 
	for(int i=0; i<n_cond; i++)
	{
		DATA_ARR[i].DCT = get_DCT(i,control_location,n_cond,DATA_ARR);
		printf("%-8d%-14s%-10s%-5d\t%lf\t%lf",i,DATA_ARR[i].SAMPLE,DATA_ARR[i].TARGET,DATA_ARR[i].count,DATA_ARR[i].CT_VAL,DATA_ARR[i].DCT);
		DATA_ARR[i].DDCT = get_DDCT(i,control_location,n_cond,DATA_ARR);
		printf("\t%lf",DATA_ARR[i].DDCT);
		DATA_ARR[i].FOLD_CHANGE = pow((double)2,(DATA_ARR[i].DDCT*-1));
        printf("\t%lf\n",DATA_ARR[i].FOLD_CHANGE);
	}
	puts("\nSaving data...");
	FILE * pRESULTS;
	pRESULTS = fopen("Results.csv","w");
	fprintf(pRESULTS,"NUMBER,SAMPLE_NAME,TARGET_NAME,AVE_CT,DCT,DDCT,FOLD_CHANGE\n");
	for(int i=0; i<n_cond; i++)
		fprintf(pRESULTS,"%d,%s,%s,%lf,%lf,%lf,%lf\n",i,DATA_ARR[i].SAMPLE,DATA_ARR[i].TARGET,DATA_ARR[i].CT_VAL,DATA_ARR[i].DCT,DATA_ARR[i].DDCT,DATA_ARR[i].FOLD_CHANGE);
	fclose(pCSV); fclose(pRESULTS);
	puts("Done!");
	
	return 0;
}
