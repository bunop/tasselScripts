#! /usr/bin/Rscript --slave --vanilla

# Un programma per calcolare la matrice triangolare dato un file (con una matrice)

# definisco un po' di roba utile
library("optparse")

# Una funzione spero più allegra per capire i parametri di ingresso
Parser = optparse::OptionParser(description="A program to make a triangoular matrix from a file")
Parser = add_option(Parser, c("-i","--inputfile"), type="character", help="Input matrix file", default="")
Parser = add_option(Parser, c("-o","--outputfile"), type="character", help="Output triangoular matrix file", default="")

#questo è quello che arriva in ingresso
argv = commandArgs(trailingOnly = TRUE)

#Se non ho passato nulla, stampo la guida ed esco
if (length(argv)==0) {
  argv = c("--help")
}

# Qui leggo tutti i parametri che ho passato
arguments = parse_args(Parser, args=argv, print_help_and_exit=TRUE, positional_arguments=FALSE)

if (arguments$inputfile == "") {
  stop("You MUST specify a matrix inputfile")
}

if (arguments$outputfile == "") {
  stop("You MUST specify a matrix outputfile")
}

#ok leggo il file in ingresso
read.table(file=arguments$inputfile,header=TRUE)->matrix

#ok dimensione della matrice (che è un dataframe!)
matrix_size = dim(matrix)

if (matrix_size[2] > 10) {
  # Questi sono i valori letti. Head delle prime 10 colonne
  head(matrix[,1:10])
} else {
  #ma se non ho 10 colonne, non mi faccio tutti sti strippi
  head(matrix)
}

#così ho tolto la prima colonna
M <- matrix[,-1]

#M è data frame, ora matrice
M <- as.matrix(M)

#Serve davvero fare un head di M?
# head(M[,1:12])

print("Calculating triangoular matrix...")

#Questa è la vera trasformazione in triangolare. Si mette NA nella parte ALTA
M[upper.tri(M,diag=TRUE)]<-NA

#debug
matrix_size = dim(M)

if (matrix_size[2] > 10) {
  # Questi sono i valori letti. Head delle prime 10 colonne
  head(M[,1:10])
} else {
  #ma se non ho 10 colonne, non mi faccio tutti sti strippi
  head(M)
}

#scrivo il file in output
write.table(M,file=arguments$outputfile,quote=FALSE,row.names=FALSE)
