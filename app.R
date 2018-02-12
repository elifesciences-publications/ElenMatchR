###ElenMatchR V0.3

library(shiny)
library(randomForest)
library(data.table)
library(plotly)
library(ggtree)
library(ggdendro)
library(DT)
library(tidyverse)
library(Biostrings)

#Preprocess genome data and store as RDS for file size / ease reasons, see: https://www.rdocumentation.org/packages/base/versions/3.4.3/topics/readRDS
#for (scaffold in list.files("resources/scaffolds/")){
#  tmp<-readDNAStringSet(paste0("resources/scaffolds/",scaffold))
#  saveRDS(tmp, paste0("resources/RDS/scaffolds/",scaffold,".RDS"))
#}
#for (gtf in list.files("resources/gtfs//")){
#  tmp<-read_tsv(paste0("resources/gtfs/", gtf), comment="##", col_names = F) %>% filter(!is.na(X2)) #exploit parsing error of seq portion
#  colnames(tmp)<-c("seqname", "source","feature","start","end","score","strand","frame","attribute")
#  saveRDS(tmp, paste0("resources/RDS/gtfs/",gtf,".RDS"))
#}
#for (matrix in list.files("/Volumes/turnbaughlab/qb3share/jbisanz/Eggerthella_Genomes_2018/po5/outs/")){
#  if(!grepl("-graph", matrix)){
#  tmp<-read_tsv(paste0("/Volumes/turnbaughlab/qb3share/jbisanz/Eggerthella_Genomes_2018/po5/outs/", matrix), col_names = T) %>% 
#    select(-1,-2,-3) %>%
#    mutate(OG=paste0("OG",gsub("pid([0-9][0-9])_cov([0-9][0-9])..+","\\1\\2",matrix),"_", 1:nrow(.))) %>%
#    select(OG, everything())
#  colnames(tmp)<-gsub("\\.faa","", colnames(tmp))
#  saveRDS(tmp, paste0("resources/RDS/matrix/",matrix,".RDS"))
#  }
#}
#tree<-ape::read.tree("/Volumes/turnbaughlab/qb3share/jbisanz/Eggerthella_Genomes_2018/nsegata-phylophlan-1d174e34b2ae/output/allstrains_jan24/allstrains_jan24.tree.nwk")
#tree<-ape::root(tree, outgroup="Escherichia_coli_K12MG1655", resolve.root=T)
#saveRDS(tree, "resources/RDS/phylotree.RDS")

#Static Vars on load
phenos<-read_tsv("Phenotypes.txt", skip=1)
tree<-readRDS("resources/RDS/phylotree.RDS")

# User interface ########## ########## ########## ########## ########## ########## ##########
ui <- navbarPage(HTML("ElenMatchR: Comparative Genomics Tool v0.3"),
                 tabPanel("Settings",
                          h1("Instructions"),
                          p("Welcome to ElenMatchR v0.3. This tool is for the linking of phenotypes to genes in Corriobacteriia and specifically Eggerthella lenta. This tool is described in Bisanz et al., `Establishing a toolkit for genetics and ecology of the Coriobacteriia and Eggerthella lenta: prevalent symbionts of the gut microbiome` (check back soon for reference). To use your own data, upload a copy of the provided template after selecting Custom_Phenotype from the Phenotype dropdown box. Pick your desired clustering thresholds (decrease for comparisons between more distantly related Coriobacteriia) and random forest parameters. After clicking execute analysis, the remaining tabs will be populated with results. Please report any problems to Jordan.Bisanz@ucsf.edu."),
                          hr(),
                          h1("Phenotypes"),
                          selectInput("phenotype", "Phenotype:", selected="Digoxin_Reduction",
                                      choices=c("Digoxin_Reduction",
                                                "Tetracycline_Resistance",
                                                "Custom_Phenotype")),
                          fileInput("phenofile","Or select Custom_Phenotype above and upload your own:", multiple = FALSE),
                          helpText("Download template", a("here.", href="https://raw.githubusercontent.com/jbisanz/ElenMatchR/master/Phenotypes.txt")),
                          hr(),
                          h1("Clustering Parameters"),
                          selectInput("pid", "Min Percent Identity", selected="60", 
                                      choices=c("30",
                                                "40",
                                                "50",
                                                "60",
                                                "70",
                                                "80",
                                                "90",
                                                "95",
                                                "99")),
                          selectInput("pcov", "Min Query Coverge", selected="80",
                                      choices=c("50",
                                                "60",
                                                "70",
                                                "80",
                                                "90",
                                                "95")),
                          hr(),
                          h1("Random Forest Parameters"),
                          numericInput("mtry", "mtry", value=0),
                          helpText("Number of variables used at each split, if zero defaults to sqrt(#features)."),
                          numericInput("ntree", "ntree", value=1000),
                          helpText("Number of trees generated."),
                          hr(),
                          h1("Output Options"),
                          numericInput("nfeats","Number of Features to display", value=15),
                          selectInput("cladogram",label="Should phylogenetic tree be shown as cladogram?", choices=c("Cladogram","PhylogeneticTree"), selected="PhylogeneticTree"),
                          helpText("This number of features will be included in output."),
                          hr(),
                          actionButton("goButton", "Execute Analysis"),
                          hr(),
                          p("Bisanz et al., 2018")
                 ),
                 tabPanel("Phenotypes", dataTableOutput("phenos"), p("This is a copy of the phenotype data that was used for prediction. Please ensure it is correct.")),
                 tabPanel("Classifier Accuracy", tableOutput("ConfusionMatrix") , p("This shows RF classifier accuracy. The accuracy is not of much importance for the purpose of feature selection but it is provided for reference.")),
                 tabPanel("Annotation Table", dataTableOutput("ImportanceTable"), p("The following table contains the predictive accuracy of each of the top variables as well as their annotation information. Sequences can be retrived from the Candidate Sequences tab.") ),
                 tabPanel("Importance Scatter Plot",  p("The top features are plotted here with their importance (mean decrease GINI). This will give an idea of relative strength of your top predictors."), downloadButton("DLImpScatter"), plotOutput("ImpScatter")),
                 tabPanel("Heat Map",  p("Heat map shows occurence between groups. See Annotation table tab for more information."), downloadButton("DLHeatMap"), plotOutput("HeatMap")),
                 tabPanel("Pan-genome Tree", p("Pan-genome tree based on gene occurence (UPGMA clustering of binary distance)."), downloadButton("DLpangenometree"), plotOutput("pangenometree")),
                 tabPanel("Phylogenetic Tree", p("Phylophlan based (~400 conserved gene sequences) phylogenetic tree (Segata et al., 2013: 10.1038/ncomms3304)."), downloadButton("DLphylotree") , plotOutput("phylotree")),
                 tabPanel("Candidate Sequences", p("Sequences corresponding to Importance Table are available for download as FASTA files."), hr(), h2("Nucleotide Sequences"),  downloadButton("DLnucleotide"), hr(), h2("Amino Acid Sequences"), downloadButton("DLaminoacid"),
                          p("FASTA headers denote the orthologous cluster from which the gene occured and annotation data.")),
                 hr()
)

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 25*1024^2)
  
  observeEvent(input$goButton, {
    withProgress(message = 'Analyzing Data', value = 0, {
      n=14
    message(
      paste(sep = "", date(), "\n",
            "protocol: ", session$clientData$url_protocol, "\n",
            "hostname: ", session$clientData$url_hostname, "\n",
            "pathname: ", session$clientData$url_pathname, "\n",
            "port: ",     session$clientData$url_port,     "\n",
            "search: ",   session$clientData$url_search,   "\n")
    )
    
    
    
  #Get and display phenotypes
message("Loading Phenotypes")
incProgress(1/n, detail = "Loading phenotypes...")

  if(is.null(input$phenofile)){
    message(paste0("Loading ", input$phenotype, " from file"))
    phenotypes<-phenos %>% dplyr::select(Strain, Genome_ID, ShortName, input$phenotype) %>% filter(!is.na(get(input$phenotype)))
    #phenotypes<-phenos[,which(colnames(phenos) %in% c("Strain", input$phenotype))]
  } else { 
    message(paste0("Loading phenotypes from ", input$phenofile[1]))
    phenofiles<-input$phenofile
    phenotypes<-read_tsv(phenofiles$datapath, skip=1) %>% dplyr::select(Strain, Genome_ID, ShortName, Custom_Phenotype) %>% dplyr::select(Strain, Genome_ID, ShortName, input$phenotype) %>% filter(!is.na(get(input$phenotype)))
  }
  output$phenos<-renderDataTable({phenotypes}, extensions = 'Buttons', options=list("dom" = 'T<"clear">lBfrtip', buttons = list('copy', 'csv', 'excel', 'pdf', 'print')))

  ###Get and collapse genes
  incProgress(2/n, detail = "Loading gene matrix...")
  matrix<-input$matrices
  message(paste0("Loading gene matrix from ", matrix$datapath ," at ", input$pcov, "% coverage and ", input$pid, " % identity."))
  genes<-readRDS(paste0("resources/RDS/matrix/pid",input$pid,"_cov",input$pcov,".proteinortho.RDS"))
  
  message(paste("Genes are of format:", class(genes), "with dimensions:", dim(genes)))
  
  incProgress(3/n, detail = "Converting and subsetting gene matrix...")
  message("Converting gene matrix")
  
  genes<-genes %>% dplyr::select(OG, which(colnames(genes) %in% phenotypes$Genome_ID))   #subset to only the strains in the analysis
  bin_genes<-genes %>% mutate_at(phenotypes$Genome_ID, function(x) if_else(x=="*", 0, 1))
  bin_genes$Observed<-apply(bin_genes[,2:ncol(bin_genes)], 1, sum) 
  
  message("Of ", length(bin_genes$Observed), " input features: ", 
          sum(bin_genes$Observed==0), " are not found in input strains and ",
          sum(bin_genes$Observed==length(phenotypes$Genome_ID)), " are found invariably leaving ",
          sum(bin_genes$Observed!=0 & bin_genes$Observed!=length(phenotypes$Genome_ID)), " possible features." )   # remove features present in all or none of the strains as they are not informative

  bin_genes<-bin_genes %>% filter(Observed!=0 & Observed!=length(phenotypes$Genome_ID)) %>% dplyr::select(-Observed)
  
  lookup<-tibble(OG=bin_genes$OG, cluster=apply(bin_genes[,2:ncol(bin_genes)],1, function(x) paste(x, collapse="|")))
  
  incProgress(5/n, detail = "Collapsing the covarying features...")
  message("Collapsing the covarying features...")
  
  derep<-bin_genes %>% as.data.frame() %>% column_to_rownames("OG") %>% as.matrix
  rownames(derep)<-lookup[match(rownames(derep), lookup$OG),]$cluster # if somethingif fishy double check this!!!!!!!
  derep<-derep[!duplicated(rownames(derep)),]
  
  message("After dereplication  of ", nrow(bin_genes), " features , there are ", nrow(derep), " features")

  ###Random Forest
  message(paste0(date(),": Carrying out random forest"))
  incProgress(6/n, detail = "Carrying out random forest...")
  
  if(input$mtry==0){RF<-randomForest(t(derep), as.factor(phenotypes[match(colnames(derep), phenotypes$Genome_ID),] %>% pull(input$phenotype)),  importance=T, ntree=input$ntree)}
  else{RF<-randomForest(t(derep), as.factor(phenotypes[match(colnames(derep), phenotypes$Genome_ID),] %>% pull(input$phenotype)),  importance=T, ntree=input$ntree, mtry=input$mtry)}
  
  message(paste0(date(),": Completed random forest"))
  incProgress(7/n, detail="Completed random forest")
  
  Importance<-as.data.frame(RF$importance)
  
  
  message("Rendering Confusion Matrix")
  incProgress(8/n, detail = "Rendering Confusion Matrix")
  
  output$ConfusionMatrix<-renderTable(cbind(rownames(RF$confusion), RF$confusion))
  
  
  ###Dereplicate
  message("Re-replicating features and sorting")
  incProgress(9/n, detail = "Re-replicating features and sorting")
  
  Importance<-Importance %>% as.data.frame %>%
    rownames_to_column("cluster") %>%
    full_join(lookup) %>%
    dplyr::select(OG, cluster, everything()) %>%
    left_join(genes) %>%
    arrange(desc(MeanDecreaseGini, decreasing=T)) %>%
    top_n(input$nfeats, MeanDecreaseGini)
  
  
  message("Annotating top features and getting sequences")
  incProgress(10/n, detail = "Annotating top features")
  
  # pull annotations from annotation line of GTF
  
  for(genome in phenotypes$Genome_ID){
    gtf<-readRDS(paste0("resources/RDS/gtfs/",genome,".gff.RDS")) %>% 
      filter(feature=="CDS") %>% 
      mutate(Gene=gsub("^ID=","", attribute) %>% gsub(";..+","", .)) %>%
      mutate(Information=paste(Gene,seqname,start,end,strand,attribute, sep="|"))
    Importance<-Importance %>% mutate_at(genome, function(x) {
              gtf[match(x,gtf$Gene),]$Information
          }
      )
  }
  output$ImportanceTable<-renderDataTable({Importance}, extensions = 'Buttons', options=list("dom" = 'T<"clear">lBfrtip', buttons = list('copy', 'csv', 'excel', 'pdf', 'print')))
  
  
  message("Making Scatter plot")
  incProgress(11/n, detail = "Plotting")

  IMPPLOT<-(ggplot(Importance %>% mutate(Rank=1:nrow(.)), aes(x=Rank,y=MeanDecreaseGini, label=OG)) 
    +  theme_bw() 
    + geom_point() 
    + geom_line() 
    + theme(axis.text.x = element_text(angle = 45, vjust = 0.1)) 
    + labs(x="Ranked Importance", y="Importance (Mean Decrease Gini)")  
    + geom_text(angle=45, size=4, alpha=0.2,hjust=0,vjust=0, color="red")
    + scale_y_continuous(limits=c(0, 1.2*max(Importance$MeanDecreaseGini)))
  )
  
  output$ImpScatter <- renderPlot({IMPPLOT}, height=700)
  
  output$DLImpScatter <- downloadHandler(
    filename = function() { paste0("ScatterPlot_", input$phenotype, "_cov", input$pcov, "_pid", input$pid, '.pdf') },
    content = function(file) {
      device <- function(..., width, height) grDevices::pdf(..., width = width, height = height)
      ggsave(file, plot = IMPPLOT, device = device, width=10, height=7.5)
    }
  )
  
  
  message("Making Heat Map")
  
  Importance %>% 
    dplyr::select(OG,phenotypes$Genome_ID) %>%
    mutate_at(phenotypes$Genome_ID, function(x) if_else(is.na(x), "Absent", "Present")) %>%
    gather(-OG, key="Genome_ID", value="Presence") %>%
    mutate(OG=factor(OG, levels=rev(Importance$OG))) %>%
    left_join(phenotypes) %>%
    ggplot(aes(x=ShortName,y=OG, fill=Presence)) +
         theme_bw() +
         geom_tile() +
         facet_grid(~get(input$phenotype), scales="free", space="free") +
         theme(axis.text.x=element_text(angle=45, hjust=1)) +
         scale_fill_manual(labels=c("Absent","Present"), values=c("aliceblue","steelblue")) +
         ylab("orthologous gene cluster") +
         xlab("strain") -> PLOT
  
  #PLOT
  output$HeatMap <- renderPlot({PLOT}, height=700)
    
  ###
  output$DLHeatMap <- downloadHandler(
    filename = function() { paste0("HeatMap_", input$phenotype, "_cov", input$pcov, "_pid", input$pid, '.pdf') },
    content = function(file) {
      device <- function(..., width, height) grDevices::pdf(..., width = width, height = height)
      ggsave(file, plot = PLOT + theme(axis.text.y=element_text(size=5)), device = device, width=10, height=7.5)
    }
  )
  ####
  
  message("Indicating Presence in Pan-genome Tree")
  
  pangenome<- bin_genes %>% 
    as.data.frame() %>% 
    column_to_rownames("OG") %>%
    t(.) %>%
    dist(., method="binary") %>%
    hclust(., method="average")

  #(ggdendrogram(pangenome) 

  #https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2
  #phenotypes$Strains<-rownames(phenotypes)
  dendr    <- dendro_data(pangenome, type="rectangle") # convert for ggplot
  clust.df <- data.frame(label=phenotypes$Genome_ID, Phenotype=phenotypes %>% pull(input$phenotype))
  # dendr[["labels"]] has the labels, merge with clust.df based on label column
  dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
  # plot the dendrogram; note use of color=cluster in geom_text(...)
  dendr$labels$label<-phenotypes[match(dendr$labels$label, phenotypes$Genome_ID),]$ShortName
  
    ggplot() +
      theme_classic() +
      geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) +
      geom_text(data=label(dendr), aes(x, y, label=label, color=Phenotype), hjust=1, vjust=0, angle=45, nudge_y=-0.01, size=3) +
      scale_y_continuous(limits=c(-0.1, 1), breaks=seq(0,1,0.2)) +
      scale_x_continuous(expand=c(0.1,0.1)) +
      ylab("Binary Distance") +
      xlab("Strain") +
      theme(axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line.y=element_blank()) +
      theme(legend.position=c(0.8,0.8)) -> pantree
      
    
  
   output$pangenometree <- renderPlot({pantree}, height=600)
  
  ###
  output$DLpangenometree <- downloadHandler(
    filename = function() { paste0("PanGenomeTree_", input$phenotype, "_cov", input$pcov, "_pid", input$pid, '.pdf') },
    content = function(file) {
      device <- function(..., width, height) grDevices::pdf(..., width = width, height = height)
      ggsave(file, plot = pantree, device = device, width=10, height=7.5)
    }
  )
  ####
  
  
  message("Indicating Presence in Phylogenetic Tree")
  
  tree$tip.label<-gsub("\\.faa","", tree$tip.label)
  filt.tree<-drop.tip(tree, tree$tip.label[!tree$tip.label %in% phenotypes$Genome_ID])
  filt.tree$tip.label<-phenotypes[match(filt.tree$tip.label, phenotypes$Genome_ID),]$ShortName
  
  if(input$cladogram=="Cladogram"){
    plottree<-ggtree(filt.tree, branch.length = "none") + geom_tiplab() #+ xlim(0, 18)
    }
  else{plottree<-ggtree(filt.tree) + geom_tiplab()}

  
  plottree<-(plottree %<+% (phenotypes %>% select(ShortName, input$phenotype)))
  plottree +
    geom_tippoint(aes_string(color=input$phenotype)) + 
    theme(legend.position=c(0.1,0.7)) +
    scale_x_continuous(expand=c(0.2,0.2)) -> plottree

  output$phylotree <- renderPlot({plottree})
  
  ###
  output$DLphylotree <- downloadHandler(
    filename = function() { paste0("PhyTree_", input$phenotype, "_cov", input$pcov, "_pid", input$pid, '.pdf') },
    content = function(file) {
      device <- function(..., width, height) grDevices::pdf(..., width = width, height = height)
      ggsave(file, plot = plottree, device = device, width=10, height=7.5)
    }
  )
  ####
  
  incProgress(12/n, detail = "Extracting Candidate Sequences")
  message("Pulling sequences")
  
  Importance %>%
    gather(key="Strain", value="Annotation", phenotypes$Genome_ID) %>%
    filter(!is.na(Annotation)) %>%
    mutate(Gene=gsub("\\|..+","",Annotation)) -> seqlookup
  
  grid<-str_split(seqlookup$Annotation, "\\|", n=6, simplify=T) %>%
    as.tibble() %>%
    setNames(c("Gene","Contig","start","stop","strand","scrap"))
    
  seqlookup<-seqlookup %>% left_join(grid) %>% arrange(Strain) %>%
   # mutate(Sstart=if_else(strand=="+", start, stop)) %>% #flip start/stops for rev strand
   # mutate(Send=if_else(strand=="+", stop, start)) %>%
    mutate(NucSeq=NA, AASeq=NA)

for(i in 1:nrow(seqlookup)){
  
  if(i==1){ #avoid repeated opening/closining of files
    contigs<-readRDS(paste0("resources/RDS/scaffolds/",seqlookup[i,]$Strain,".scaffolds.fasta.RDS"))
    names(contigs)<-gsub("..+\\|","", names(contigs))
    loadedstrain=seqlookup[i,]$Strain
  } else if(i==1 | loadedstrain!=seqlookup[i,]$Strain){ #avoid repeated opening/closining of files
    contigs<-readRDS(paste0("resources/RDS/scaffolds/",seqlookup[i,]$Strain,".scaffolds.fasta.RDS"))
    names(contigs)<-gsub("..+\\|","", names(contigs))
    loadedstrain=seqlookup[i,]$Strain
  }
  
  seqlookup[i,]$NucSeq<-
  as.character(subseq(contigs[[seqlookup[i,]$Contig]], start=as.numeric(seqlookup[i,]$start), end=as.numeric(seqlookup[i,]$stop)))
  seqlookup[i,]$NucSeq<-if_else(seqlookup[i,]$strand=="+", seqlookup[i,]$NucSeq, as.character(reverseComplement(DNAString(seqlookup[i,]$NucSeq))))
  seqlookup[i,]$AASeq<-as.character(translate(DNAString(seqlookup[i,]$NucSeq)))
}
  
  seqlookup<-seqlookup %>% arrange(desc(OG)) %>%
    arrange(desc(MeanDecreaseGini))
  
  AAseqs<-AAStringSet(seqlookup$AASeq)
  names(AAseqs)<-paste(seqlookup$OG, seqlookup$Annotation)
  
 NUCseqs<-AAStringSet(seqlookup$NucSeq)
  names(NUCseqs)<-paste(seqlookup$OG, seqlookup$Annotation)
  
  output$DLnucleotide <- downloadHandler(
    filename = function() { paste0(input$phenotype, "_cov", input$pcov, "_pid", input$pid, '.ffn') },
    content = function(file) {
      writeXStringSet(NUCseqs, file)
    }
  )
  
  
  output$DLaminoacid <- downloadHandler(
    filename = function() { paste0(input$phenotype, "_cov", input$pcov, "_pid", input$pid, '.faa') },
    content = function(file) {
      writeXStringSet(AAseqs, file)
      #device <- function(..., width, height) grDevices::pdf(..., width = width, height = height)
      #ggsave(file, plot = pantree, device = device, width=10, height=7.5)
    }
  )
  
  
  
  message("Analysis Complete")
  })
})
  

  }


# Run the application 
shinyApp(ui = ui, server = server)

