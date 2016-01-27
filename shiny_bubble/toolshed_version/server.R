shinyServer(function(input, output, session) {  

################################################################################  
  jsnetwork <- reactive({
    str_x=paste0(input$main.x)
    str_y=paste0(input$main.y)
    str_color=paste0(input$main.color)
    str_size=paste0(input$main.size)
    str_cutoff= paste0(input$main.cutoff)
    str_label= input$main.label
    scl_size = input$main.scale
    main.data2 <- main.data[!(main.data$PreyGene %in% input$main.exclude),]
    main.data2 <- subset(main.data2, main.data2[(colnames(main.data2)=="log2(FoldChange)")] >= input$main.change)
    cytoscape <- subset(main.data2, SaintScore>=as.numeric(str_cutoff), 
                        select = c(str_x,str_y,"Bait","PreyGene",str_size))
    id <- c(as.character(cytoscape$PreyGene),as.character(cytoscape$Bait))
    source <- as.character(cytoscape$Bait)
    name <- c(as.character(cytoscape$PreyGene),as.character(cytoscape$Bait))
    target <- as.character(cytoscape$PreyGene)
    node.data <- data.frame(id,name)
    edge.data <- data.frame(source,target)
    node.data$color <- rep(input$node.color, nrow(node.data))
    node.data$color[node.data$name %in% cytoscape$Bait] <- input$bait.color # do the same but with shape
    node.data$shape <- rep(input$node.shape, nrow(node.data))
    node.data$shape[node.data$name %in% cytoscape$Bait] <- input$bait.shape# do the same but with shape
    edge.data$color <- rep(input$edge.color, nrow(edge.data))
    network <- createCytoscapeJsNetwork(node.data, edge.data)
    rcytoscapejs(network$nodes, network$edges, showPanzoom=TRUE, 
                 layout=input$net.layout, highlightConnectedNodes=FALSE)
  })
################################################################################
  jsnetwork_flat <- reactive({
    str_x=paste0(input$main.x)
    str_y=paste0(input$main.y)
    str_color=paste0(input$main.color)
    str_size=paste0(input$main.size)
    str_cutoff= paste0(input$main.cutoff)
    str_label= input$main.label
    scl_size = input$main.scale
    main.data2 <- main.data[!(main.data$PreyGene %in% input$main.exclude),]
    main.data2 <- subset(main.data2, main.data2[(colnames(main.data2)=="log2(FoldChange)")] >= input$main.change)
    cytoscape <- subset(main.data2, SaintScore>=as.numeric(str_cutoff), 
                        select = c(str_x,str_y,"Bait","PreyGene",str_size))
    id <- c(as.character(cytoscape$PreyGene),as.character(cytoscape$Bait))
    source <- as.character(cytoscape$Bait)
    name <- c(as.character(cytoscape$PreyGene),as.character(cytoscape$Bait))
    target <- as.character(cytoscape$PreyGene)
    node.data <- data.frame(id,name)
    edge.data <- data.frame(source,target)
    node.data$color <- rep("#888888", nrow(node.data))
    node.data$color[node.data$name %in% cytoscape$Bait] <- "#FF0000"
    network <- toJSON(as.list(createCytoscapeJsNetwork(node.data, edge.data)))
  })
################################################################################
  bubblebeam <- reactive({
    str_x=paste0(input$main.x)
    str_y=paste0(input$main.y)
    str_color=paste0(input$main.color)
    str_size=paste0(input$main.size)
    str_cutoff= paste0(input$main.cutoff)
    str_label= input$main.label
    scl_size = input$main.scale
    main.data2 <- main.data[!(main.data$PreyGene %in% input$main.exclude),]
    main.data2 <- subset(main.data2, main.data2[(colnames(main.data2)=="log2(FoldChange)")] >= input$main.change)
    
    
    c <- subset(main.data2, SaintScore>=as.numeric(str_cutoff), select = c(str_x,str_y,"Bait","PreyGene",str_size))
    colnames(c) <- c("x","y","Bait","PreyGene","size")
    p <- ggplot(data=c, x=x, y=y,size=size)+ geom_point(data=c, aes(x=x, y=y,size=size), color=input$bubble.color) + scale_size(range=scl_size)
    p <- p + geom_point(data=c, aes(x=x, y=y, size=size), colour=input$outline.color, shape=21) + labs(x=str_x, y=str_y, size=str_size)
    if(length(levels(c$Bait) > 1)) {p <- p + facet_wrap(~Bait, scales="free_y")}
    if(str_label== TRUE & length(c$x)>=1) {p <- p + geom_text(data=c, aes(x=x,y=y,label=PreyGene, size=max(size)/5, vjust=0, hjust=0),colour=input$label.color, show_guide=FALSE)}
    if(input$theme== "b/w") {p <- p + theme_bw()}
    
    
    if(str_color=="crapome") {
      a <- subset(main.data2, CrapomePCT <80 & SaintScore >=as.numeric(str_cutoff), select = c(str_x,str_y,"Bait","PreyGene",str_size,"CrapomePCT"))
      b <- subset(main.data2, CrapomePCT >=80 & SaintScore >=as.numeric(str_cutoff), select = c(str_x,str_y,"Bait","PreyGene",str_size,"CrapomePCT"))
      colnames(a) <- c("x","y","Bait","PreyGene","size", "CrapomePCT")
      colnames(b) <- c("x","y","Bait","PreyGene","size","CrapomePCT")
      p <- ggplot(data=a, x=x, y=y,size=size) + geom_point(data=a,aes(x=x,y=y,size=size),color=I(input$filt.color)) +
        scale_size(range=scl_size) + geom_point(data=a, aes(x=x, y=y, size=size), colour=input$outline.color, shape=21)
      if(length(levels(a$Bait) > 1)) {p <- p + facet_wrap(~Bait, scales="free_y")}
      if(str_label== TRUE & length(a$x)>=1) {p <- p + geom_text(data=a, aes(x=x,y=y,label=PreyGene, size=max(size)/5, vjust=0, hjust=0),colour=input$label.color)}
      
      p <- p + geom_point(data=b, aes(x=x, y=y, size=size, color=CrapomePCT)) + 
        scale_colour_gradient(limits=c(80, 100), low=input$filt.color, high=input$bubble.color) + 
        labs(colour="CRAPome Probability \nof Specific Interaction (%)", x=str_x, y=str_y,size=str_size) + 
        geom_point(data=b, aes(x=x, y=y, size=size), colour=input$outline.color, shape=21)
      if(str_label== TRUE & length(b$x)>=1) {p <- p + geom_text(data=b, aes(x=x,y=y,label=PreyGene, size=max(size)/5, vjust=0, hjust=0),colour=input$label.color, show_guide=FALSE)}
      if(input$theme== "b/w") {p <- p + theme_bw()}
    }
    p
}) 
################################################################################
table_display <- reactive({
  str_cutoff= paste0(input$main.cutoff)
  main.data2 <- main.data[!(main.data$PreyGene %in% input$main.exclude),]
  main.data2 <- subset(main.data2, main.data2[(colnames(main.data2)=="log2(FoldChange)")] >= input$main.change)
  table <- subset(main.data2, SaintScore>=as.numeric(str_cutoff))
  table
}) 
################################################################################
pathway_graph <- reactive({
  str_val <- input$path_x
  str_cutoff= paste0(input$main.cutoff)
  main.data2 <- main.data[!(main.data$PreyGene %in% input$main.exclude),]
  main.data2 <- subset(main.data2, main.data2[(colnames(main.data2)=="log2(FoldChange)")] >= input$main.change)
  table <- subset(main.data2, SaintScore>=as.numeric(str_cutoff))
  preys <- unique(as.character(table$PreyGene))
  EG_IDs <- list()
  for(i in 1:length(preys)){
    EG_IDs[i] <- mygene::query(preys[i])$hits$entrezgene[1]
  }

  pathways <- enrichKEGG(EG_IDs, 
                         organism = paste0(input$path_org), 
                         pvalueCutoff = as.numeric(input$path_pval), # NOT SURE IF pvaluecutoff WORKS
                         pAdjustMethod = paste0(input$path_adj),
                         readable=TRUE)
  
  pathways <- as.data.frame(summary(pathways))
  pathways$x <- factor(pathways$Description,levels=rev(pathways$Description))

  if(str_val=="pvalue"){
    p <- ggplot(data=pathways, aes(y=-log(pvalue), x=x)) + geom_bar(stat="identity") + coord_flip()}
  if(str_val=="p.adjust"){
    p <- ggplot(data=pathways, aes(y=-log(p.adjust), x=x)) + geom_bar(stat="identity") + coord_flip()}
  p
})
################################################################################
pathway_table <- reactive({
  str_val <- input$path_x
  str_cutoff= paste0(input$main.cutoff)
  main.data2 <- main.data[!(main.data$PreyGene %in% input$main.exclude),]
  main.data2 <- subset(main.data2, main.data2[(colnames(main.data2)=="log2(FoldChange)")] >= input$main.change)
  table <- subset(main.data2, SaintScore>=as.numeric(str_cutoff))
  preys <- unique(as.character(table$PreyGene))
  EG_IDs <- list()
  for(i in 1:length(preys)){
    EG_IDs[i] <- mygene::query(preys[i])$hits$entrezgene[1]
  }

  pathways <- enrichKEGG(EG_IDs, 
                         organism = paste0(input$path_org), 
                         pvalueCutoff = as.numeric(input$path_pval),
                         pAdjustMethod = paste0(input$path_adj),
                         readable=TRUE)
  summary(pathways)
})
################################################################################
ontology_graph <- reactive({
  str_val <- input$path_x
  str_cutoff= paste0(input$main.cutoff)
  main.data2 <- main.data[!(main.data$PreyGene %in% input$main.exclude),]
  main.data2 <- subset(main.data2, main.data2[(colnames(main.data2)=="log2(FoldChange)")] >= input$main.change)
  table <- subset(main.data2, SaintScore>=as.numeric(str_cutoff))
  preys <- unique(as.character(table$PreyGene))
  EG_IDs <- list()
  for(i in 1:length(preys)){
    EG_IDs[i] <- mygene::query(preys[i])$hits$entrezgene[1]
  }

  pathways <- enrichGO(EG_IDs, 
                         organism = paste0(input$path_org),
                       pvalueCutoff = as.numeric(input$path_pval),
                         pAdjustMethod = paste0(input$path_adj),
                         readable=TRUE,
                        ont = input$GO_ont)
  pathways <- as.data.frame(summary(pathways))
  pathways$x <- factor(pathways$Description,levels=rev(pathways$Description))

  if(str_val=="pvalue"){
    p <- ggplot(data=pathways, aes(y=-log(pvalue), x=x)) + geom_bar(stat="identity") + coord_flip()}
  if(str_val=="p.adjust"){
    p <- ggplot(data=pathways, aes(y=-log(p.adjust), x=x)) + geom_bar(stat="identity") + coord_flip()}
  p
})
################################################################################
ontology_table <- reactive({
  str_val <- input$path_x
  str_cutoff= paste0(input$main.cutoff)
  main.data2 <- main.data[!(main.data$PreyGene %in% input$main.exclude),]
  main.data2 <- subset(main.data2, main.data2[(colnames(main.data2)=="log2(FoldChange)")] >= input$main.change)
  table <- subset(main.data2, SaintScore>=as.numeric(str_cutoff))
  preys <- unique(as.character(table$PreyGene))
  EG_IDs <- list()
  for(i in 1:length(preys)){
    EG_IDs[i] <- mygene::query(preys[i])$hits$entrezgene[1]
  }

  pathways <- enrichGO(EG_IDs, 
                       organism = paste0(input$path_org), 
                       pvalueCutoff = as.numeric(input$path_pval),
                       pAdjustMethod = paste0(input$path_adj),
                       readable=TRUE,
                       ont = input$GO_ont)
  summary(pathways)
  })
################################################################################
   # Plot Bubble Graph
   output$bubbles=renderPlot({
     print(bubblebeam())
   })
  # Render Cytoscape Network
  output$network=renderRcytoscapejs({
    jsnetwork()
  })
  # Render SAINT Table
  output$table=renderDataTable({
    table_display()
  })
  # Render Pathway Analysis Graph
  output$pathPlot=renderPlot({
  print(pathway_graph())
  })
  # Render Pathway Analysis Table
  output$pathTable = downloadHandler(
  filename = function() {
    paste(Sys.Date(), "_KEGG.txt",sep='')}, 
  content= function(file){
    x = pathway_table()
    write.table(x, file, quote=FALSE,
                sep = "\t", row.names=FALSE)})
  # Render Gene Ontology Graph
  output$ontPlot=renderPlot({
  print(ontology_graph())
  })
  # Render Gene Ontology Table
  output$ontTable = downloadHandler(
  filename = function() {
    paste(Sys.Date(), "_GO_",as.character(input$GO_ont),".txt",sep='')}, 
  content= function(file){
    x = ontology_table()
    write.table(x, file, quote=FALSE,
                sep = "\t", row.names=FALSE)})
  # Save SAINT table
  output$saveTable = downloadHandler(
    filename = function() {
    paste(Sys.Date(), "_dataTable.txt",sep='')}, 
    content= function(file){
      x = table_display()
      write.table(x, file, quote=FALSE,
        sep = "\t", row.names=FALSE)})
      
  # Download
   output$main.down = downloadHandler(filename = function() {
                                      paste(Sys.Date(), "_bubbleplot",input$main.file,sep='')},
                                      content = function(file){
                                        p=bubblebeam()
                                        ggsave(file, plot = p,scale=2,dpi=300)})

observeEvent(input$saveImage, { #acknowledge Augustin Luna (Sloan Kettering) for help getting this working
  # NOTE: Message cannot be an empty string "", nothing will happen    
  session$sendCustomMessage(type="saveImage", message="NULL")
})

output$saveJSON = downloadHandler(filename="network.json",
  content = function(file){
    write(jsnetwork_flat(),file=file)
    })

 })