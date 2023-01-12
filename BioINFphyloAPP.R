#-------------------------
# Project name: Phylogenetic tree maker 
# App name: bioINFphylo
# Producer: Sal Jahromi
# ------------------------





# install.packages("shiny")
# install.packages("shinydashboard")
# install.packages("shinythemes")
# install.packages("shinybusy")
# 
# library(shiny)
# library(shinythemes)
# library(shinybusy)
# library(shinydashboard)
# 
# #Installing the package dendextend
# install.packages("dendextend")
# install.packages("ggmsa")
# 
# install.packages("BiocManager")
# BiocManager::install("Biostrings")
# library("Biostrings")
# 
# install.packages("phangorn")
# library("phangorn")
# 
# #Installing the package phytools
# install.packages("phytools")
# 
# #Installing the package taxize
# install.packages("taxize")
# 
# install.packages("seqinr")
# install.packages("tidyverse")
# #lapply is part of the apply family of functions, one of those functions that act as a for loop
# BiocManager::install("msa")
# 



package_list <- c("shiny", "shinydashboard","shinythemes","shinybusy","dendextend","ggmsa","BiocManager","Biostrings","phangorn","phytools","taxize","seqinr","tidyverse", "msa")
new.packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


lapply(c("ape","seqinr","Biostrings", "tidyverse","msa", "phangorn",
         "phytools","stringr","taxize","shiny","shinythemes","shinybusy","shinydashboard"),
       library, character.only = TRUE)




ui<- dashboardPage(skin = "green",
                   dashboardHeader(title = "bioINFphylo"),
                   dashboardSidebar(),
                   dashboardBody(
                     
                     box(
                       
                       title = ("Inputs"),
                       
                       radioButtons("seqType", "Sequence Type " , choices = c("Amino Acid" = "AA", "Nucleotide (DNA)" = "DNA")),
                       
                       fileInput(
                         inputId = "file",
                         label = "Drag and drop here",
                         multiple = TRUE,
                         buttonLabel = "Browse...",
                         placeholder = "No file selected"
                       ),
                       
                       
                       tags$h3("Input:"),
                       
                       selectInput("msaMethod", "MSA Method: ", choices = c("ClustalW", "ClustalOmega", "Muscle")),
                       
                       selectInput("treeReconstructionMethod", "Tree Reconstruction Method", choices = c("UPGMA", "Neighbor Joining")),
                       
                       numericInput("bsIterations", "Bootstrap Iterations ", "500"),
                       
                       checkboxInput("showFileSeq", "Show File Sequences ", FALSE),
                       checkboxInput("bootstrapBool", "Bootstrap Tree ", FALSE),
                       
                       actionButton("action", "Generate Tree!"),
                       id = "box1",
                       collapsible = TRUE,
                     ),
                     
                     box(
                       title = "Output Results",
                       verbatimTextOutput("empty"),
                       verbatimTextOutput("txtout"),
                       plotOutput("txtout2"),
                       id = "box2",
                       collapsible = TRUE,
                       collapsed = FALSE,
                     )
                   )
)


#server side
server <- function(input, output, session){ 
  
  #Text in output section
  output$empty <- renderText("Nothing to show...")
  
  #When "Generate tree" button is pressed, code is run.
  observeEvent(input$action,{
    #the output section text changes to this
    output$empty <- renderText("Phylo Result")
    
    #This receives the inputted information and stores into variables
    msaMethod<- input$msaMethod #the msa method
    seqType <- input$seqType #the type of sequence(amino acid/ dna)
    bsIterations <- input$bsIterations # number of interations for bootstrap
    treeReconstructionMethod <- input$treeReconstructionMethod #tree reconstruction method
    
    showFileSeq <- input$showFileSeq #checkbox to whether show file seq in output box or not
    bootstrapBool <- input$bootstrapBool #checkbox to whether bootstrap graph or not
    
    # print(seqType)
    # print(bsIterations)
    
    show_modal_spinner(spin = "circle",text = "Reading file...") #loading circle
    Sys.sleep(1)
    substituition_model = "" #will store the sub model depending on the sequence type (JTT or JC69)
    
    if(seqType == "AA"){ #if file has AA seq, then use readAAstringset
      update_modal_spinner(text = "Reading Amino Acids...")
      Sys.sleep(1)
      fileR <- readAAStringSet(input$file$datapath)
      substituition_model = "JTT"
    }
    else if(seqType == "DNA"){#else if file has DNA seq, then use readDNAstringset
      update_modal_spinner(text = "Reading DNA...")
      Sys.sleep(1)
      fileR <- readDNAStringSet(input$file$datapath)
      substituition_model = "JC69"
    }
    
    
    # print(substituition_model)  
    
    if(showFileSeq == TRUE){ #show the file sequence if the checkbox is true
      output$txtout <- renderPrint({
        print(fileR)
      })
    }
    
    
    # names(hba_aa) <- str_remove_all(str_split_fixed(names(hba_aa),pattern = "\\[", n=2)[,2],"\\]")
    
    
    update_modal_spinner(text = "Generating MSA...")
    Sys.sleep(1)
    #generate msa
    fileR_msa <- msa(fileR, method = msaMethod) #msa made with inputted msa method
    
    
    update_modal_spinner(text = "Seqinr alignment & phydat element...")
    Sys.sleep(1)
    #Convert to a seqinr alignment and then th a phyDat class element
    fileR_phy <- as.phyDat(msaConvert(fileR_msa, "seqinr::alignment"), type = seqType)
    
    
    
    #making distance matrix
    update_modal_spinner(text = "building distance matrix...")
    fileR_dml <- dist.ml(fileR_phy, model = substituition_model)#builds distance matrix with substitution_model based on whether file is AA or DNA
    
    
    
    
    
    
    
    ## If the tree is neighbor joining tree, use the following code made for NJ trees
    update_modal_spinner(text = "Reconstructing tree...")
    Sys.sleep(1)
    
    if(treeReconstructionMethod == "Neighbor Joining"){
      fileR_tree <- midpoint(phangorn::NJ(fileR_dml))    #generating neighbor joining tree and midpointing
      
      
      
      if(bootstrapBool == FALSE){
        output$txtout2 <- renderPlot({
          plot(fileR_tree)
        })
      }
      else if(bootstrapBool == TRUE){ #if bootstrap check box is selected -> bootstrap the tree
        update_modal_spinner(text = "Bootstrapping...")
        Sys.sleep(1)
        
        fileR_tree_bs <- bootstrap.phyDat(fileR_phy,
                                            FUN = function(x)NJ(dist.ml(x,model = substituition_model)), bs = bsIterations )
      
        
        output$txtout2 <- renderPlot({
          plotBS(fileR_tree,fileR_tree_bs,"phylogram")
        })
        
      }
      
    }
    
    ## If the tree is upgma tree, use the following code made for upgma trees
    
    else if(treeReconstructionMethod == "UPGMA"){
      
      
      fileR_tree <- (phangorn::upgma(fileR_dml)) #no midpoint because upgma tree is already midpointed
      
      if(bootstrapBool == FALSE){
        output$txtout2 <- renderPlot({
          plot(fileR_tree)
        })
      }
      
      else if(bootstrapBool == TRUE){ #if bootstrap check box is selected -> bootstrap the tree
        update_modal_spinner(text = "Bootstrapping...")
        Sys.sleep(1)
        
        fileR_tree_bs <- bootstrap.phyDat(fileR_phy,
                                            FUN = function(x)upgma(dist.ml(x,model = substituition_model)), bs = bsIterations )
        
        
        output$txtout2 <- renderPlot({
          plotBS(fileR_tree,fileR_tree_bs,"phylogram")
        })
        
      }
    }
    
    update_modal_spinner(text = "Almost Done...")
    Sys.sleep(1)
    
    
    
    
    remove_modal_spinner()
    
  })
  
}

shinyApp(ui, server)
