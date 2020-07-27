#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int main(void) 
{
    char buf[4096], token[20],filename[256];
    int row_count = 0, field_count, in_double_quotes = 0, token_pos = 0;
    int i,s=0,t=0,c=0;
    printf("Enter CSV file to anaylize: ");
    scanf("%s",filename);
    FILE *fp = fopen(filename, "r");

    char SAMPLE[100][256],SAMPLELIST[100][256];
    char TARGET[100][256],TARGETLIST[100][256];

    double CT[100];

    if (!fp) 
    {
        printf("Can't open file");
        return 0;
    }

    printf("\nProcessing CSV file.. ");

    while (fgets(buf, 4096, fp)) 
    {
        row_count++;

    	if (row_count == 1) 
		{
    	    continue;
    	}

        field_count = 0;
        i = 0;

        do 
		{
            token[token_pos++] = buf[i];

            if (!in_double_quotes && (buf[i] == ',' || buf[i] == '\n')) 
	   		{
                token[token_pos - 1] = 0;
                token_pos = 0;
				field_count++;
				if (field_count == 1)
				{
					strcpy(SAMPLE[s],token);
					strcpy(SAMPLELIST[s],token);
					char token[256];
					s++;
				}
				if (field_count == 2)
				{
					strcpy(TARGET[t],token);
					strcpy(TARGETLIST[t],token);
					char token[256];
					t++;
				}
				if (field_count == 3)
				{
					CT[c] = strtod(token,NULL);
					char token[256];
					c++;
				}
		
        	}
			if (buf[i] == '"' && buf[i + 1] != '"') 
			{
				token_pos--;
				in_double_quotes = !in_double_quotes;
			}

			if (buf[i] == '"' && buf[i + 1] == '"')
			{
				i++;
			}

		}
		while (buf[++i]);
    }


    printf("\nCalculating Average CT values...");

	
    int len=0,n=0;
    for(i=0;i<100;i++)
    {
	    if(CT[i] != 0)
		    len++;
    }	

    int cmp,j,k,nT,nS;
    n=0;
    i=0;
    cmp=0;
    nT=len;

    for(i=0;i<len;i++)
	    for(j=i+1;j<nT;)
	    {
		    cmp = strcmp(TARGETLIST[i],TARGETLIST[j]);
		    if(cmp == 0)
		    {
			    for(k=j;k<nT-1;k++)
				    strcpy(TARGETLIST[k],TARGETLIST[k+1]);
			    nT--;
		    }
		    else
			    j++;
	    }
	
    j=0;
	k=0;
    i=0;
    cmp=0;
    nS=len;

    for(i=0;i<len;i++)
	    for(j=i+1;j<nS;)
	    {
		    cmp = strcmp(SAMPLELIST[i],SAMPLELIST[j]);
		    if(cmp == 0)
		    {
			    for(k=j;k<nS-1;k++)
				    strcpy(SAMPLELIST[k],SAMPLELIST[k+1]);
			    nS--;
		    }
		    else
			    j++;
	    }

    s=0;
    t=0;
    i=0;
    cmp=0;
    int compT,compS,n_reps=1;
    char targsamp[100][512];
    
    for(n=0;n<len-1;n++)
    {
	    strcpy(targsamp[n],TARGET[n]);
	    strcat(targsamp[n],SAMPLE[n]);
	    cmp = strcmp(targsamp[1],targsamp[n]);
	    if(cmp == 0)
	    {
		    n_reps++;
	    }
	    
    }   

    float AVECTs[len],val,sum=0,mean;
    char TARG[100][256],SAMP[100][256];
    int v=-1;
	i=0;

    for(n=0;n<len;n++)
    {
	    for(t=0;t<nT;t++)
	    {
		    for(s=0;s<nS;s++)
		    {
			    compT = strcmp(TARGET[n],TARGETLIST[t]);
			    compS = strcmp(SAMPLE[n],SAMPLELIST[s]);
			    if(compT == 0 & compS == 0)
			    {
				    val = CT[n];
				    sum = val+sum;
				    mean = sum/n_reps;
				    i++;
			    } else if (i >= n_reps)
				{
					i = 0;
					sum=0;
					v++;
					AVECTs[v] = mean;
					strcpy(TARG[v],TARGET[n]);
					strcpy(SAMP[v],SAMPLE[n]);
				}
		    }
			
	    }
		AVECTs[v+1] = mean;
		strcpy(TARG[v+1],TARGET[n]);
		strcpy(SAMP[v+1],SAMPLE[n]);
    }

	int len_ave_data;
	len_ave_data = v+2;

    printf("\nGetting DCT values...");
    char control_s[256],control_t[256];
    float DCT[len],ave_cntl_ct;
    compT = 0;
    compS = 0;
    strcpy(control_s,SAMP[0]);
    strcpy(control_t,TARG[0]);

    printf("Using %s and %s as sample and gene control...",control_s,control_t);
	
    for(i=0;i<len_ave_data;i++)
    {
		compT = strcmp(control_t,TARG[i]);
		if(compT == 0)
		{
			ave_cntl_ct = AVECTs[i];
			DCT[i] = 0;
		} 
		else
		{
			DCT[i] = AVECTs[i] - ave_cntl_ct;
		}
    }

    printf("\nGetting DDCT values...");
    s=0;
    t=0;

	//getting DDCT values

    float DCTvals[len],DDCT[len],DDCTs[len], FOLD_CHANGE[len];

    for(i=0;i<len;i++)
    {
	    if(DCT[i] != 0)
	    {
		    for(t=0;t<n;t++)
		    {
			  compT = strcmp(TARGETLIST[t],TARG[i]);
			  if(compT == 0)
			  {
			  	compS = strcmp(control_s,SAMP[i]);
				if(compS == 0)
				{
					DCTvals[t] = DCT[i];
				}
				else
				{
					DDCT[i] = DCT[i] - DCTvals[t];
				}
			  }				  
		    }
	    }
	    else
	    {
		    DDCT[i] = 0;
	    }
    }

    printf("\nCalculating fold change...\n");    

    printf("\nSample\t Target\t Mean\t\t DCT\t\t DDCT\t\t FOLD_CHANGE\t \n\n");
    for(i=0;i<len_ave_data;i++)
    {
	    FOLD_CHANGE[i] = pow(2,(-1*DDCT[i]));
	    printf("%s\t %s\t %lf\t %lf\t %lf\t %lf\t \n",SAMP[i],TARG[i],AVECTs[i],DCT[i],DDCT[i],FOLD_CHANGE[i]);
    }
	
    fclose(fp);

	char savedfile[256];

	printf("\nEnter filename to save CSV: ");
	scanf("%s", savedfile);
	printf("\nSaving %s",savedfile);
	FILE *sp;
	sp = fopen(savedfile,"w+");
	fprintf(sp,"SAMPLE, TARGET, CT, DCT, DDCT, FOLD_CHANGE\n");
	for(i=0;i<len_ave_data;i++)
	{
		fprintf(sp,"%s,%s,%lf,%lf,%lf,%lf", SAMP[i],TARG[i],AVECTs[i],DCT[i],DDCT[i],FOLD_CHANGE[i]);
		fprintf(sp,"\n");
	}

	fclose(sp);

	printf("\n\nDone! Press CTRL+C to quit!\n");

	scanf("%s");

    return 0;
}
