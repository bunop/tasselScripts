#! /usr/bin/Rscript --slave --vanilla

#Definisci con subset quali traits vuoi rappresentare

# definisco un po' di roba utile
library("optparse")
suppressPackageStartupMessages(library("compare"))
#library("compare")

# direi di caricare le librerie dallo script qqman
source("~/R/Scripts/qqman.r")

# Ho bisogno di definire la tipologia delle colonne per il file GLM e MLM
GLM_columns <- c("Trait", "Marker", "Locus", "Locus_pos", "marker_F", "marker_p", "markerR2", "markerDF", "markerMS", "errorDF", "errorMS", "modelDF", "modelMS" )
MLM_columns <- c("Trait", "Marker", "Locus", "Site", "df", "F", "p", "errordf", "markerR2", "Genetic.Var", "Residual.Var", "X.2LnLikelihood")

#Da qui mi tiro tutti gli strippi
# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# trailingOnly=TRUE means that only arguments after --args are returned
# if trailingOnly=FALSE then you got:
# [1] "--no-restore" "--no-save" "--args" "2010-01-28" "example" "100"

# Una funzione spero più allegra per capire i parametri di ingresso
Parser = optparse::OptionParser(description="A program to make tassel manhattan and qqplot")
Parser = add_option(Parser, c("-i","--inputfile"), type="character", help="Tassel results file")
Parser = add_option(Parser, c("-W","--width"), type="integer", help="Image width [default %default]", default=480)
Parser = add_option(Parser, c("-H","--height"), type="integer", help="Image heigth [default %default]", default=480)
Parser = add_option(Parser, c("-o","--outfile"), type="character", help="Image output file basename [required]", default="")
Parser = add_option(Parser, c("-t","--typegraph"), type="character", help="Graph type (manhattan, qqplot) [default %default]", default="manhattan")
Parser = add_option(Parser, c("--suggestiveline"), type="double", help="Suggestive line (log10, ex 1e-5) [default %default]", default=0 )
Parser = add_option(Parser, c("--genomewideline"), type="double", help="Genomewide line (log10, ex 5e-8) [default %default]", default=0 )
Parser = add_option(Parser, c("--titleprefix"), type="character", help="Title prefix", default="" )
Parser = add_option(Parser, c("--snpsfile"), type="character", help="A list of SNPs to highlight", default="")
Parser = add_option(Parser, c("--gff3file"), type="character", help="The gff3 file with gene annotations", default="")

#questo è quello che arriva in ingresso
argv = commandArgs(trailingOnly = TRUE)

#Se non ho passato nulla, stampo la guida ed esco
if (length(argv)==0) {
  argv = c("--help")
}

# Qui leggo tutti i parametri che ho passato
arguments = parse_args(Parser, args=argv, print_help_and_exit=TRUE, positional_arguments=FALSE)

# debug
# print(argv)
# print (arguments)
# print (arguments$outfile)

# Attenzione: devo passare il file in output
if (arguments$outfile == "") {
  stop("You MUST specify a basename for PNG output files")
}

if (arguments$typegraph != "manhattan" && arguments$typegraph != "qqplot" ) {
  stop("Graph type must be 'manhattan' or 'qqplot' ")
}

#Se non passo argomenti per le righe suggerite, allora le imposto come false
if (arguments$suggestiveline == 0) {
  arguments$suggestiveline = FALSE
} else {
  arguments$suggestiveline = -log10(arguments$suggestiveline)
}

if (arguments$genomewideline == 0) {
  arguments$genomewideline = FALSE
} else {
  arguments$genomewideline = -log10(arguments$genomewideline)
}

#Se devo evidenziare degli SNPs, lo faccio
snps_to_highlight = NULL

if (arguments$snpsfile != "") {
  snps_to_highlight = scan(arguments$snpsfile, character())
}

#posso aver passato un file gff3 con le annotazioni geniche
gff3_data = NULL

if (arguments$gff3file != "") {
  gff3_data = read.table(arguments$gff3file, header=F, sep="\t")
  
  #Assegno le colonne al gff3
  gff3_colnames = c("reference_sequence", "source", "type", "start_position", "end_position", "score", "strand", "phase", "attributes")
  colnames(gff3_data) <- gff3_colnames
  
  #Per il momento mi preparo a gestire solamente le tracce di tipo GENE
  gff3_data = subset(gff3_data, gff3_data$type=="gene")
  
  #la reference sequence dovrebbe essere numerica (nel caso del riso con 12 cromosomi numerici, non ci sono problemi)
  gff3_data = transform(gff3_data, reference_sequence=as.numeric(sub("Chr", "", reference_sequence, ignore.case=T)))
  
}

# i parametri che ho passato stanno in un dataframe.
tassel_data <- read.table(arguments$inputfile, header=T, sep="\t")

# ok provo a leggere le colonne che ci sono dentro
columns = names(tassel_data)

if ( compare(columns, GLM_columns, equal=F)$result == TRUE ) {
  print("GLM model detected")
  
  # Ora cambiano i nomi delle colonne per renderli utilizzabili dalle funzioni qqman
  colnames(tassel_data)[2] = "SNP"
  colnames(tassel_data)[3] = "CHR"
  colnames(tassel_data)[4] = "BP"
  colnames(tassel_data)[6] = "P"
  
  #Assegno il prefisso per il titolo del grafico, se mi manca
  if (arguments$titleprefix == "") {
    arguments$titleprefix = "GLM-P"
  }
  
} else if ( compare(columns, MLM_columns, equal=F)$result == TRUE ) {
  print("MLM model detected")
  
  # Ora cambiano i nomi delle colonne per renderli utilizzabili dalle funzioni qqman
  colnames(tassel_data)[2] = "SNP"
  colnames(tassel_data)[3] = "CHR"
  colnames(tassel_data)[4] = "BP"
  colnames(tassel_data)[7] = "P"
  
  #Assegno il prefisso per il titolo del grafico, se mi manca
  if (arguments$titleprefix == "") {
    arguments$titleprefix = "MLM"
  }
  
} else {
  stop("Make sure your data derive from tassel")
}

head(tassel_data)

# Adesso vorrei vedere quali traits sono definiti
traits = unique(tassel_data$Trait)

# Per ogni trait, un grafico apposta
for (i in 1:length(traits)) {
  filename = paste0(arguments$outfile, "_", traits[i], ".png")
  tassel_trait = subset(tassel_data, tassel_data$Trait==traits[i])
  
  #ok, dirigo le mie storie sun un file PNG
  png(filename, width=arguments$width, height=arguments$height)
  
  #costruisco la stringa per il titolo
  graph_title = paste(arguments$titleprefix, traits[i])
  
  # ok, ora posso chiamare le funzioni che mi interessano (poi dovrò decidere quale funzione chiamare e con che parametri)
  if (arguments$typegraph == "manhattan") {
    # debug
    print(paste("Faccio il manhattan plot per", traits[i], "e lo salvo in", filename))
    
    #setto i margini nuovi (bottom, left, top, right)
    old_par = par(oma=c(1,2,1,3)) #par(oma=c(2,3,2,4))
    
    #Questi sono i margini interni. Sono importanti, altrimenti le scritte vengono tagliate
    #http://rankexploits.com/musings/2011/margin-control-in-r-oma-and-mar-and-the-vanishing-axis-label/
    par(mar=c(5,5,5,5))
    
    #faccio il grafico
    #cex.main raddoppia le dimensioni del titolo principale
    #cex.lab sono le label degli assi
    #cex.axis sono le etichette dei valori degli assi
    manhattan(tassel_trait, colors=c("green","blue"), pch=18, annotate=snps_to_highlight, gff3_data=gff3_data, genomewideline=arguments$genomewideline, suggestiveline=arguments$suggestiveline, main=graph_title, cex.main=2, cex.lab=1.5, cex.axis=1.5)
    
    #ripristino i margini originali
    par(old_par)
    
  }
  
  else {
    # debug
    print(paste("Faccio il qqplot per", traits[i], "e lo salvo in", filename))
    
    #setto i margini nuovi (bottom, left, top, right)
    old_par = par(oma=c(1,2,1,3)) #par(oma=c(2,3,2,4))
    par(mar=c(5,5,5,5))
    
    #cex.main raddoppia le dimensioni del titolo principale
    #cex.lab sono le label degli assi
    #cex.axis sono le etichette dei valori degli assi
    qq(tassel_trait$P, main=graph_title, cex.main=2, cex.lab=1.5, cex.axis=1.5 )
    
    #ripristino i margini originali
    par(old_par)
    
  }
  
  # per stampare i messaggi
  warnings()
  
  #fine
  dev.off()
  
}

