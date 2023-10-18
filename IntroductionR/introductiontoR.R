# "#" is used to create comments in your code and 
# are not executed by the R interpreter

## where are we?
getwd()

## the helpful help
?getwd()

#we need to go the folder "IntroductionR", remember the "" are needed for the paths
setwd("IntroductionR/")
getwd()

#What files are in the folder?
list.files()

######################
# Basic calculations #
######################

4 + 4
4 / 2
4 * 4
2 ^ 4 

################
# Conditional  #
################

2+2==5# = means another thing
-1.9>-2
4>7

##################
# Variables in R #
##################
#Most of the time R acts on things that are stored in variables
#You can see all your variables in the Environment

#you can assign variables using <- or =

x<- 7

x

x+3
2 ^ x
x * x

# check what type of data is contained within a variable
class(x)

w <- "europa"

w*2

#because
class(w)


############
## Vectors #
############

#are used to store multiple items into a single variable. 
#A one-dimensional object holding multiple items 

#The c() operator is used to concatenate or combine elements!!!

y <- c(5, 6, 7)
y

##############
# Dataframe ##
##############
#Variables can also hold tables
#Dataframes are two-dimensional objects of rows and columns


#let's do another vector and then combine it with our previous vector

z <- c(8, 9, 10)
our_table <- data.frame(y, z)

our_table

class(our_table)


# The wonderful world of indexing
# Subsetting tables or vectors based on the index position

#we specify the vector name, and then put in brackets [ ] 
#the position(s) we are interested in.


y # the whole vector
y[1] # the first item
y[2] # second item
y[3] # third item
y[-1]

#We can also ask for multiple by using the c() function 


y[c(1,3)] # specifying items 1 and 3

#or also with conditional

y[y >= 6] #Give me all the values of y where y is greater than or equal to 6

y[!y >= 6] # ! means opposite

### subset tables
our_table
our_table[2, 2]


our_table[ , 2] # subset all rows, but only the second column

our_table[3, ] # only row 3, but both columns

## subset by names

colnames(our_table) 

#specify a column we want to pull from a dataframe based on the column name, 
#is to enter the table variable name, followed by a $

our_table$z

#or

our_table[c(2,3), "z"]

#Notice that we gave the column name within quotes here. 
#This is because we want R to know to just interpret the text 
#and not to look for an object stored in the variable name "z".

###################################
# Reading in and writing out data #
###################################

?read.table

gene_annotations_tab <- read.table("gene_annotations500.txt", sep = "\t", header = TRUE)

head(gene_annotations_tab) #print first 5 rows


colnames(gene_annotations_tab)

dim(gene_annotations_tab)


### extract all the rows that don't contain missing values in the KO_ID column
#is.na() function is used to identify missing or NA (Not Available) values

head(gene_annotations_tab$KO_ID)
is.na(head(gene_annotations_tab$KO_ID))

### try to do it by yourself

# KEGG_only_tab <- gene_annotations_tab[] .... continue
 








# the first part of the [] means: get all the rows where the KO_ID 
# column value is not NA (it would give us all the ones that are NA 
# if we didn't include the `!` upfront) 
# the second part after the comma,   we are providing nothing specifying 
#which columns, which as mentioned above means to take all of them



KEGG_only_tab <- gene_annotations_tab[!is.na(gene_annotations_tab$KO_ID), ]


dim(gene_annotations_tab) 
dim(KEGG_only_tab)

## write.table()
?write.table()


write.table(KEGG_only_tab, "KEGG_annotated.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
list.files() # checking it is there now





