#Student Name: Pierre Wensel
#Assignment 2
#Course: Bioinformatics GR55 [MUADO-12UV]
#Date: October 27,2023


# I first install dplyr package because I use between() function in problem#4. I used between because
#it was unclear if "between" in problem#4 was to be or not to be interpreted as 
#left X>= LEFT VALUE and X<= RIGHT VALUE:
install.packages("dplyr")    
# Load dplyr package     
library("dplyr") 

#Exercise 1:
#DEFINING VARIOUS FUNCTIONS
#1
function1<-function(p,alpha, x0){
  Q<-x0*((1-p)^(-1/(alpha-1)))
  return(Q)
}

#2
function2<-function(p,alpha, x0, lower.tail=TRUE){
  if(lower.tail==TRUE) {
    Q<-x0*((1-p)^(-1/(alpha-1)))
    return(Q)
  }else{ 
    p <- (1-p)
    Q<-x0*((1-p)^(-1/(alpha-1)))
    return(Q)
    }
}
#3
function3<-function(p,alpha, x0, lower.tail=TRUE){
  if(lower.tail==TRUE) {
    #Calling function1 here
    return(function1(p,alpha,x0))
  }else{ 
    p <- (1-p)
    return(function1(p,alpha,x0))
  }
}

#4
function4<-function(p,alpha, x0, lower.tail=TRUE){
  #Using stopifnot for data validation. If conditions for parameters not met, function will not execute
  stopifnot(between(p,0,1),(alpha > 1),(x0 > 0))
  if(lower.tail==TRUE) {
    return(function1(p,alpha,x0))
  }else{ 
    p <- (1-p)
    return(function1(p,alpha,x0))
  }
}

#5
function5<-function(alpha, x0, lower.tail,n){
  #is.integer(n)is NOT TRUE so it was not used here for data validation of paramter n
  stopifnot(n>=1) 
  x <- c(1:n)
  #Initializing vector
  #Q=c() was not used because I want to allocate only enough memory for vector of size n
  Q <- vector("numeric", n)
  #Iterating with for loop
  for(i in x) {
    #Generaing random number between 0 and 1 from distribution each iteration
    p<-runif(1)
    #Outputting p values to check calculations:
    #print(p)
    #Calling function4 and generating output vector
    Q[i] <- function4(p, alpha, x0, lower.tail)
    #{Q=c(Q,function4(p, alpha, x0, lower.tail))} ALSO WORKS
  }
    return(Q)
}

#6
function6<-function(alpha, x0, lower.tail,n){
  stopifnot(n>=1)
  #Generating a vector of p values with dimensions set by parameter n
  p<-runif(n)
  #Outputting p values to check calculations:
  #print(p)
  #Q<-apply(p,2,function4(p, alpha, x0, lower.tail))DOES NOT WORK DUE TO function4 not being considered a function
  #Q<-apply(p,2,function4) DOES NOT WORK AS dim(x) must have positive length
  #Calling function4 and passing in vector p, and avoiding use of for-loop
  Q<-function4(p, alpha, x0, lower.tail) 
  return(Q)
}

#7
function7<-function(alpha, x0, lower.tail,n){
  stopifnot(n>=1)
  p<-runif(n)
  #Outputting p values to check calculations:
  #print(p)
  Q<-function4(p, alpha, x0, lower.tail)
  #Outputting p values to check calculations:
  #print(Q)
  #return(c(mean(Q), sd(Q))) Returning a vector object would also work if it was permitted by assignment
  #list includes mean and standard deviation
  return(list(mean = mean(Q), standard_deviation = sd(Q))) 
}

#8
function8<-function(alpha, x0, lower.tail,n){
  stopifnot(n>=1)
  p<-runif(n)
  #Outputting p values to check calculations:
  #print(p)
  Q<-function4(p, alpha, x0, lower.tail)
  #print(Q)
  #Quantile() function creates sample quantiles within a data set with probability[0, 1]. 
  #Such as first quantile is at 0.25[25%], second is at 0.50[50%], and third is at 0.75[75%]
  return(quantile(Q, prob=c(.25,.5,.75), type=1)) 
}

#9
function9<-function(alpha, x0, lower.tail,n){
  stopifnot(n>=1)
  p<-runif(n)
  #Outputting p values to check calculations:
  #print(p)
  Q<-function4(p, alpha, x0, lower.tail)
  color<- readline(prompt="Enter your histogram color as a lower-case word with no parenthesis:")
  #Returning histogram using user-entered color
  return(hist(Q, xlab="Q Values", ylab="frequency", main=mtext(~italic("Colored histogram")), col=color))
}

#Exercise 2

#STEP 0: Import the ‘WA2_OmicsDataExercisesII.txt’ into R 
dataset1<-read.table("WA2_OmicsDataExercisesII.txt", header=T, sep="")
#PLEASE NOTE: THIS dataset1 is used for other ASSIGNMENT STEPS below that rely on its existence:

#STEP1: Find the length of dataset1 columns:
#First, I used the following to get brief overview of data
#str(dataset1)
#summary(dataset1)
#ncol(dataset1)
#nrow(dataset1)
#Then, I used he following to get length of dataset1 dataframe columns
sapply(dataset1, function(x) sum(complete.cases(x)))
#colSums(!is.na(dataset1)) ALSO CAN BE USED

#STEP2:
#Writing a function to accept filename as argument
colmean<-function(genefile){
  #Reading in the file to create dataframe
  df<-read.table(genefile, header=T, sep="")
  #subsetting into df_subset to retrieve only numeric columns from original dataframe df:
  df_subset<-df[,unlist(lapply(df,is.numeric))]
  #Returning a vector of column mean values:
  return(apply(df_subset,2,mean))
}

#Calling function colmean() to calculate mean of columns from file  passed in as argument
colmean("WA2_OmicsDataExercisesII.txt")

#STEP3:
#Subsetting only numeric columns for Sample87 rows
subdataset2<-dataset1[which(dataset1$Sample_ID=="Sample87"), c(3,4,5)]
#subdataset3<-dataset1[dataset1$Sample_ID=="Sample87", c(3,4,5)] ALSO WORKS
#Creating vector of minimum values
minvector<-apply(subdataset2,2,function(x) min(x))
#Creating vector of maximum values
maxvector<-apply(subdataset2,2,max)
#Combining the 2 vectors to form a matrix
newmatrix<-cbind(minvector,maxvector)

#STEP4:
#Calculating the mean expression values for each of the 7 unique genes
tapply(dataset1$Expression_Value,dataset1$Gene_Name, mean)

#STEP5:
#Using for loop to generate list of unique gene names
unique_gene_names=list()
for(i in dataset1$Gene_Name){
  x=unique(dataset1$Gene_Name[i])
  unique_gene_names[i]=x
}
#Displaying the list
unique_gene_names

#STEP6:
#Using tapply function to instead generate list of unique gene names using work from step4
unique_genes_list<-as.list(names(tapply(dataset1$Expression_Value,dataset1$Gene_Name, mean)))
#typeof(unique_genes_list)

#STEP7
#Using if/for-loop approach to calculate mean of dataframe dataset1's gene expression values exceeding gt value (Note: generated from assignment text file previously in STEP0)
avg_gt_version1<-function(x,gt){
  sum<-0
  count<-0
  for(i in 1:length(x)){
    if(x[i]>gt){
      sum<-sum+x[i]
      count<-count+1
      }
  }
  avg<-sum/count
  return(avg)
}

#Using non-for-loop approach to calculate mean of dataframe dataset1's gene expression values exceeding gt value 
avg_gt_version2<-function(x,gt){
  x<-x[x>gt]
  avg<-mean(x)
  return(avg)
}


#Calling the if and for-loop-based function:
avg_gt_version1(dataset1$Expression_Value, 0.25)
#Measuring time for execution of version1 function:
system.time(avg_gt_version1(dataset1$Expression_Value, 0.25))

#Calling the non-if and for-loop-based function:
avg_gt_version2(dataset1$Expression_Value, 0.25)
#Measuring time for execution of version2 function:
system.time(avg_gt_version2(dataset1$Expression_Value, 0.25))

print("Version 2 function without if/for-loop was expected to be faster function, but both functions took 0 time ")

