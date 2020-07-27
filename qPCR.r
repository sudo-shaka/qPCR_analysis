#!/usr/bin/Rscript

file <- ('ENTER NAME OF CSV FILE HERE!!!')


data <- read.csv(file)
data <- cbind.data.frame(data$Sample.Name,data$Target.Name, data$CT)# extract important colulmns

# get names of genes/targets
n_uniq <- as.integer(rapply(data[1:2],function(x)length(unique(x))))
sample_names <- list(); gene_names <- list()
for(x in 1:n_uniq[1]){sample_names[x] <- as.character(unlist(unique(data[1]))[x])}
for(y in 1:n_uniq[2]){gene_names[y] <- as.character(unlist(unique(data[[2]]))[y])}
sample_names <-unlist(sample_names); gene_names <-unlist(gene_names)

# make new data frame for AVE_CT
ave_data <- as.data.frame(unique(data[1:2]))
len_data <- unname(colSums(!is.na(data)))[1]
len_ave_data <- unname(colSums(!is.na(unique(data[1:2])[1])))
ave_data <- cbind.data.frame(ave_data,1:len_ave_data)
ave_data[[3]][1:as.integer(colSums(!is.na(unique(data[1:2]))))[1]] <- 0
colnames(ave_data) <- c("Sample.Name","Target.Name","Average.CT"); rownames(ave_data) <- NULL


#Calculating Average CT and putting into data.frame
replicates <- len_data/len_ave_data; value <- rep(NA,len_data)

print("Calculating Average CT values per sample...")

for(i in 1:len_data)
{
	for(j in 1:length(gene_names))
	{
		for(k in 1:length(sample_names))
		{
			if(as.character(data[[1]][i]) == sample_names[k] && as.character(data[[2]][i]) == gene_names[j])
			{
				value[i] <- data[[3]][i]
				ave_data[[3]][j+(length(gene_names)*(k-1))] <- mean(as.numeric(value),trim=0,na.rm=T)
			}
	  		else if(sum(!is.na(value)) >= replicates)
      			{
        			rm(value)
				value <- rep(NA,len_data)
      			}
		}
	}
}

warnings("Error with getting Ave_CT values")
DCT_data <- ave_data; DCT_data <- cbind.data.frame(ave_data,1:len_ave_data)
colnames(DCT_data) <- c("Sample.Name","Target.Name","AVE.CT","DCT")

#Calcuating DCT
Control_Target <- as.character(gene_names[1])
Control_Sample <- as.character(sample_names[1])
#Control_Sample <- c("NAME_1","NAME_2"...)
print(paste("Calculating DCT values... Assuming",Control_Target,"Control"))

for(i in 1:as.numeric(colSums(!is.na(DCT_data)))[1])
{
	if(as.character(DCT_data[[2]][i]) == Control_Target)
	{
		ave_cntl_ct <- DCT_data[[3]][i]; DCT_data[[4]][i] <- ""
	}
	else
	{
		DCT_data[[4]][i] <- DCT_data[[3]][i]-ave_cntl_ct
	}
}

warnings("[!!] Cannot get DCT")

#Calculating DDCT
print(paste("Calculating DDCT values... Assuming",Control_Sample,"control sample.."))

DDCT_data <- DCT_data;DCT <- list()
DDCT_data <- cbind.data.frame(DCT_data,1:len_ave_data,1:len_ave_data,1:len_ave_data)
colnames(DDCT_data) <- c("Sample.Name","Target.Name","AVE.CT","DCT","DDCT","FOLD CHANGE","Up/Down?")

for(samp in 1:length(Control_Sample))
{
  for(i in 1:as.numeric(colSums(!is.na(DDCT_data))[1]))
  {
    if(DCT_data[[4]][i] != "")
    {
      for(j in 1:length(gene_names))
      {
        if(gene_names[j] == as.character(DDCT_data[[2]][i]))
        {
          if(as.character(DDCT_data[[1]][i]) == Control_Sample[samp])
          {
            DCT[j] <- as.numeric(DDCT_data[[4]][i])
            DDCT_data[[5]] <- as.numeric(DDCT_data[[4]][i]) - as.numeric(DCT[j])
          } 
          else 
          {
            DDCT_data[[5]][i] <- as.numeric(DDCT_data[[4]][i]) - as.numeric(DCT[j])
          }
        }
      }
    }
    else
    {
      DDCT_data[[5]][i] <- 0
    }
  }
  #calculating fold change
  print("Calculating Fold Change...")
  for(x in 1:len_ave_data)
  {
    if(DDCT_data[[4]][x] != "")
    {
      DDCT_data[[6]][x] <- 2^(-1*DDCT_data[[5]][x])
      DDCT_data[[6]][x] <- round(DDCT_data[[6]][x],2)
    } 
    else 
    {
      DDCT_data[[6]][x] <- 1
      DDCT_data[[4]][x] <- 0
    }  
  }

  for(x in 1:len_ave_data)
  {
  	if(as.numeric(DDCT_data[[6]][x]) <= 0.5){DDCT_data[[7]][x] <- "Down"}
  	else if (as.numeric(DDCT_data[[6]][x]) >= 1.5){DDCT_data[[7]][x] <- "Up"}
  	else{DDCT_data[[7]][x] <- NA}
  }
 
  #Saving File
  print(DDCT_data[1:6])
  print("Saving calculated data...")
  write.csv(DDCT_data,file=paste(Control_Sample[samp],"Calculated_fold_changes.csv"))
  print("saved!")
}

warnings("Cannot get DDCT/Fold Change")

height = list(); xval = list(); targ <- DDCT_data[[2]][1]

for(j in 1:length(gene_names))
{
	for(i in 1:len_ave_data)
	{
		if(gene_names[j] == DDCT_data[[2]][i] && DDCT_data[[2]][i] != Control_Target)
		{
			height[i] <- DDCT_data[[6]][i]	
			xval[i] <- as.character(DDCT_data[[1]][i])
			targ_b <- targ
			targ <- DDCT_data[[2]][i]
			if(targ != targ_b)
			{
			  
			  if(as.character(typeof(height)) == "double")
			  {
			    barplot(height[!is.na(height)], ylab=paste("Fold Change: ",targ_b))
			    axis(1,a=1:length(height[!is.na(height)]) , labels = xval[!is.na(xval)])
			  }
			  height <- rep(NA,len_ave_data)
			  xval <- rep(NA,len_ave_data)
			}
		}
		
	}
}
if(as.character(typeof(height)) == "double")
{
  barplot(height[!is.na(height)], ylab=paste("Fold Change: ",targ_b))
  axis(1,a=1:length(height[!is.na(height)]) , labels = xval[!is.na(xval)])
}

warnings()
