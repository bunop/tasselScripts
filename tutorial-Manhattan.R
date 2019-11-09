# Un allegro tutorial per fare un manhattan plot
# http://gettinggeneticsdone.blogspot.it/2011/04/annotated-manhattan-plots-and-qq-plots.html

#carica prima le funzioni che ti servono per lavorare
source('~/R/Scripts/qqman.r')

#Adesso scarica il file di dati esempio del tutorial
temp <- tempfile()
download.file("http://people.virginia.edu/~sdt5z/0STABLE/plink.assoc.txt.gz",temp)
results <- read.table(gzfile(temp), header=T)
unlink(temp)

#Guarda un po come è fatto il file che hai scaricato:
head(results)

#   CHR        SNP        BP A1     F_A    F_U A2  CHISQ       P     OR
# 1   1 rs10495434 235800006  T 0.16670 0.1951  G 0.2428 0.62220 0.8250
# 2   1  rs6689417  46100028  C 0.04255 0.0000  T 3.4840 0.06195     NA
# 3   1  rs3897197 143700035  G 0.03191 0.0000  A 2.5980 0.10700     NA
# 4   1  rs2282450 202300047  T 0.41840 0.3659  A 0.5154 0.47280 1.2470
# 5   1   rs567279  66400050  0 0.00000 0.0000  T     NA      NA     NA
# 6   1 rs11208515  64900051  C 0.23960 0.2805  G 0.3861 0.53430 0.8082

# Mi ricordano un po' i file che ho con Tassel. Quello che serve però sono
# le colonne CHR per il cromosoma, SNP che è l'id dello SNP, BP che è l'allegra
# posizione dove trovo lo SNP e P che rappresenta il P-value, il valore che
# voglio rappresentare.
head(subset(results, select=c("SNP", "CHR", "BP", "P")))

#          SNP CHR        BP       P
# 1 rs10495434   1 235800006 0.62220
# 2  rs6689417   1  46100028 0.06195
# 3  rs3897197   1 143700035 0.10700
# 4  rs2282450   1 202300047 0.47280
# 5   rs567279   1  66400050      NA
# 6 rs11208515   1  64900051 0.53430

# Questa è la chiamata base per il manhattan plot
manhattan(results)

# Ci sono dei parametri che si riferiscono ai grafici base. Questi non sono
# segnati nella funzione manhattan, ma si passano e vengono iterpretati in modo
# opportuno all'interno di essa, come pch. Vedi ad esempio la documentazione
# per points
manhattan(results,pch=18)

# puoi passare una lista di colori per questo grafico, ci sono dei parametri poi
# per settare le soglie di threshold. In questo caso spengo tutto:
manhattan(results, colors=c("black","#666666","#CC6600"), pch=20, genomewideline=F, suggestiveline=F)

# Si può anche specificare una lista di SNPs da colorare in modo differente dagli
# altri. Questo link non funziona, provo a fissarli io a mano
# snps_to_highlight <- scan("http://www.StephenTurner.us/snps.txt", character())
snps_to_highlight <- results$SNP[7550:7750]

# Il grafico con gli SNPs annotati
manhattan(results, annotate=snps_to_highlight, pch=20, main="Manhattan Plot")

# Questa serve per rappresentare solo uno zoom sul cromosoma:
manhattan(subset(results, CHR==1), pch=20, annotate=snps_to_highlight, main="Chromosome 11")

# Questa invce stampa il qqplot
qq(results$P)

# Questo è invece un'altro tutorial dello stesso autore, che puoi trovare a questo indirizzo:
# http://gettinggeneticsdone.blogspot.it/2010/03/create-annotated-gwas-manhattan-plots.html
# ma preferisco comunque la funzione di prima, anche se più lenta

# This requires ggplot2
require(ggplot2)

# First, load these functions from source: Queste funzioni le hai già in locale
#source("http://dl.dropbox.com/u/66281/0_Permanent/qqman.r")

# Next, load your PLINK results file to a data frame (anche questi dati li
# hai già scaricati e letti nella variabile results):
#mydata=read.table("plink.qassoc", header=TRUE)

# Assuming you already have a vector of rs numbers to highlight. Lui li
# aveva già definiti, io invece li definisco qui
ImportantSNPs = c(3821815, 1851665, 1621816, 1403694, 1656922,  166479)

# Call the manhattan function, with annotate=T.
# The SNPlist argument takes the list of SNPs to highlight.
# Save the plot to an object. Devi usare però una funzione che
# si basa sulla libreria ggplot2
myplot=gg.manhattan(results,annotate=T,SNPlist=ImportantSNPs)

# Finally, save the plot in the current directory using ggsave()
ggsave("manhattan.png",myplot,w=12,h=9,dpi=100)

# Provo ora a utilizzare il modulo GAP, per capire se ci guadagno qualcosa come tempistiche 
# oppure no. Carico la libreria che mi interessa:
library(gap)

# provo a scrivere direttamente un immagine
png("gap_mhtplot.png")

# La funzione che voglio usare ha bisogno solo di tre colonne, cromosomi posizioni e p values
gap_data <- subset(na.omit(results), select=c("CHR", "BP", "P"))

# Ora i tizi definiscono le caratteristiche del grafico. Così creo una lista di 
# 12 elementi, ripetendo "gray10", "gray50" per 6 volte
color <- rep(c("gray10", "gray50"),6)

# Questo comando influisce su come viene rappresentato il grafico. In questo caso
# las=2 imposta le scritte sugli assi in verticale
# xpd=TRUE dice di plottare nella figura e non nel plot (?)
# cex.axis The magnification to be used for axis annotation relative to the current setting of cex.
# cex: serve per riscalare le scritte nel grafico

# é bene però registrare il valore vecchio di PAR
old_par = par()

# Ok assegno il valore nuovo
par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)

# Questa è la funzione che imposta i colori e lo spessore delle linee per il manhattan plot
# Ci puà essere anche un'altra funzione che ha effetto in verticale sui dati. Questa funzione
# viene passata come parametro alla funzione per fare il grafico
#ops <- mht.control(colors=color,yline=1.5,xline=3)
ops <- mht.control(colors=color)

# La funzione per graficare. Pch devono essere i simbolini rappresentati
mhtplot(gap_data,ops,pch=19)

# Quando hai finito coi grafici, devi chiudere il file png
dev.off()

# Ripristina i vecchi valori di par
par(old_par)

# Mi sembra che rispetto alla funzione precedente, questa faccia un po' cagare. Provo a realizzare il qqplot

# Un file dove salvare le immagini
png("qqunif.png")

# Un comando pre creare un vettore di una distribuzione uniforme
#u_obs <- runif(1000)
#r <- qqunif(u_obs,pch=21,bg="blue",bty="n")

# Ma io ho già qualcosa che voglio rappresentare
r <- qqunif(gap_data$P,pch=21,bg="blue",bty="n")

# Questo sembra essere un allegro gruppo per creare una lista di elementi con un certo valore,
# cioè crea vettore di TRUE/FALSE in base al valore dell'i-esimo elemento. Stiamo parlando di
# r che era il grafico qqplot
u_exp <- r$y
hits <- u_exp >= 2.30103

# i punti sopra questo valore, diventano verdi
points(r$x[hits],u_exp[hits],pch=21,bg="green")
dev.off()

# Perlomeno questo grafico è un po' più piacevole. Ho provato poi a caricare il pacchetto postgwas
# ma sembra avere dei problemi a livello delle dipendenze.


