#Student: Pierre Wensel
#Course: Bioinformatics GR55-MUADO-12UV
#Filename: Intro_R_Wensel.R
#Assignment: Exercise 0
#Please Note: All answers Provided in Clear Text. Additional Comments in "Green" Text Provided After Answer For Highlighting Rationale and Alternative Approaches

#16.1 Part 1:

#Code to execute a script called “myscript.R”
source("myscript.R")

#Code to assign the value A to a variable x
x<-"A"


#Code to generate a sequence from 7 to 30 with increment 3
seq(7, 30, by = 3)

#Code to obtain information about function glm
#  ?glm is another commmand to obtain information about a function glm
help(glm)

#Code to list all the objects in the current environment
ls()

#Code to remove all objects
rm(list = ls())

#Code to specify the following path to the working directory: C:
setwd("C:/")

#getwd() to check the current working directory after setting it

#Create a vector x containing the numbers 1, 2, 1, 1, 1, 2
x = c(1,2,1,1,2) 

#Create a vector y containing the words yes, no, no, yes, no
y= c("yes", "no", "no", "yes", "no")

#Compute the number of elements in vector y
length(y)

#Code to obtain the sequence of integer numbers from 10 to 25
10:25

#Use the function rep() to generate the sequence 1, 2, 1, 2, 1, 2
rep(c(1,2), 3)

#Code to generate the sequence 1, 1, 1, 2, 2, 2
rep(c(1,2), each=3)

#Code to generate a sequence containing 7 yes and 5 no
##This outputs "yes", "no" appropriate times
rep(c("yes","no"), times=c(7,5)) 

#paste(toString(rep(c("yes"),7)), toString(rep(c("no"),5)), sep = ", ") ##this outputs "yes, no" appropriate times
#cat(toString(rep(c("yes"),7)), toString(rep(c("no"),5)), sep = ", ") ##This outputs yes, no appropriate times

#Code to obtain the sequence 40, 35, 30, 25, 20, 15, 10
seq(from = 40, to = 10, by = -5)

#16.2 Part 2:

#Input the text file using read.table, assigning the input to a variable pdata.
pdata <- read.table("pData.txt", header=TRUE, sep="")

#Find the class of the variables pheno and sex. Convert them into factors using as.factor.
class(pdata$pheno)
class(pdata$sex)
as.factor(pdata$pheno)
as.factor(pdata$sex)

#Show the 10 first values of the variable “age”.
pdata$age[1:10]


#Repeat the previous values, each 3 times.
rep(pdata$age[1:10], each=3)

#Create a new data.frame “pdata_subset” containing the first 20 rows.
pdata_subset<-pdata[1:20,]


#Add in the previous dataset a new column of random values “Values”, that goes from 0.05 to 0.95.
pdata_subset$Values<-runif(nrow(pdata_subset), min=0.05, max=0.95)

#pdata_subset$values<-sample(0.05:0.95, size = nrow(pdata_subset), replace = TRUE) #This does not work


#Create a new matrix “Data_matrix” containing the information of the variables “X1”, “X3” and “X5”, 
#uniquely for the last 10 observations of the original “pdata” dataset.
Data_matrix<-unique(as.matrix(pdata[(nrow(pdata)-9):nrow(pdata),which((colnames(pdata) %in% c("X1","X3","X5")))]))

#Data_matrix<-as.matrix(pdata[1:3,which((colnames(pdata)=="X1")|(colnames(pdata)=="X3")|(colnames(pdata)=="X5"))]) #One alternative approach I was taking

#Print a sentence indicating the dimensions (rows and columns) of the matrix. 
#Use the function print, paste, nrow and ncol to do so. 
#The sentence should be “The matrix has XXXX rows and XXX columns”.

print(paste("The matrix has ", toString(nrow(Data_matrix)), " rows and ", toString(ncol(Data_matrix)), " columns", sep=" "))

#Note; toString() was not necessary for this example but could be in other examples.
#Thank you










