
library("phyloseq")
library("ggplot2")
library("scales")
library("grid")

##loading biom
rich_sparse_biom = file.path("/home/paolo/R/Scripts/sorted_all_lib_merged_by_MSC_otu_table_with_metadata.biom")
biom = import_biom(rich_sparse_biom,  parseFunction=parse_taxonomy_default)

##inspection biom
#print(biom)
ntaxa(biom)
nsamples(biom)
taxa_names(biom)
sample_names(biom)
sample_data(biom)
sample_variables(biom)
otu_table(biom)[1:5, 1:5]

##Operate on biom
colnames(tax_table(biom)) <- c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus_Species_VT")

#Commentate perchè se carico le storie dai file CSV non mi serve riassegnare queste storie (verificare)
#~ sample_data(biom)$Management <- factor(sample_data(biom)$Management, levels = c("upland", "lowland"))
#~ sample_data(biom)$MaturationStage <- factor(sample_data(biom)$MaturationStage, levels = c("tillering", "panicle-formation", "milky-stage-maturation"))

#print ordered biom data
pK = plot_bar(biom, "Kingdom", fill = "Kingdom", facet_grid = Management ~ MaturationStage)
warnings()
pK + geom_bar(aes(color = Kingdom, fill = Kingdom), stat = "identity", position = "stack")
warnings()
pKi = plot_bar(biom, "Kingdom", fill = "Kingdom", facet_grid = MaturationStage  ~ Management)
warnings()
pKi + geom_bar(aes(color = Kingdom, fill = Kingdom), stat = "identity", position = "stack")
warnings()

#~ pO = plot_bar(biom, "Order", fill = "Order", facet_grid = Management ~ MaturationStage)
#~ warnings()
#~ pO + geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack")
#~ warnings()
#~ pOi = plot_bar(biom, "Order", fill = "Order", facet_grid = MaturationStage  ~ Management)
#~ warnings()
#~ pOi + geom_bar(aes(color = Order, fill = Order), stat = "identity", position = "stack")
#~ warnings()
#~ 
#~ pF = plot_bar(biom, "Family", fill = "Family", facet_grid = Management  ~ MaturationStage)
#~ warnings()
#~ pF + geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack")
#~ warnings()
#~ pFi = plot_bar(biom, "Family", fill = "Family", facet_grid = MaturationStage  ~ Management)
#~ warnings()
#~ pFi + geom_bar(aes(color = Family, fill = Family), stat = "identity", position = "stack")
#~ warnings()
