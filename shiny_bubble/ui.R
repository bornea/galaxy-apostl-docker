shinyUI(
  pageWithSidebar(

  headerPanel(img(src="APOSTL_icon.png", height = 100, width = 150)),
    
    sidebarPanel(
             sliderInput("main.cutoff", "Saint Score Cutoff", min=0, max=1, value=0.8),
             sliderInput("main.change", "Log2(Fold Change) Cutoff", 
                         min=round(min(main.data[(colnames(main.data)=="log2(FoldChange)")]),1),
                         max=round(max(main.data[(colnames(main.data)=="log2(FoldChange)")]),1), 
                         value=round(min(main.data[(colnames(main.data)=="log2(FoldChange)")]),1)
                         ),
             checkboxInput(inputId = "main.label",label="Bubble Labels",value = FALSE),
             selectInput("main.x","X axis",selected = "ln(NSAF)",choices=c("ln(NSAF)","SpecSum", "log2(FoldChange)", "SaintScore", "logOddsScore", "-log10(BFDR)")),
             selectInput("main.y","Y axis",selected = "log2(FoldChange)",choices=c("ln(NSAF)","SpecSum", "log2(FoldChange)", "SaintScore", "logOddsScore", "-log10(BFDR)")),
             selectInput("main.size", "Size", selected = "SpecSum", choices=c("ln(NSAF)","SpecSum", "log2(FoldChange)", "SaintScore", "logOddsScore")),
             sliderInput("main.scale", "Scale Bubble Size", min=0.1, max=100, value=c(1,10)),
             selectInput("main.exclude", "Click or search to select proteins to exclude", multiple=TRUE, choices=preys),
             selectInput("main.file", "Bubble Plot File Type", choices=c(".pdf",".png",".tif"), selected=".png"),
             downloadButton('main.down', 'Download Bubble Plot'),
             actionButton("saveImage", "Download network as PNG"),
             downloadButton("saveJSON", "Download network as JSON"),
             helpText(
                p("APOSTL support is provided by the Haura and Rix labs:"),
                p("- Adam Borne (Adam.Borne@moffitt.org)\n\n"),
                p("- Brent Kuenzi (Brent.Kuenzi@moffitt.org)\n\n"),
                p("- Paul Stewart (Paul.Stewart@moffitt.org)"))
             ),
  mainPanel(
  tabsetPanel(
    tabPanel("Bubble Graph",
    plotOutput("bubbles",width="100%",height="500px", click = "plot1_click"),
    column(4,
    radioButtons("main.color","Color",choices=c("crapome","fixed"),selected="crapome",inline=TRUE),
    selectInput("bubble.color", "Bubble Color", multiple=FALSE, choices=colors, selected="#FF0000")),
    column(4,
    selectInput("outline.color", "Outline Color", multiple=FALSE, choices=c("white","black"), selected="black"),
    selectInput("filt.color","CRAPome Filtered Bubble Color",choices=colors,selected="#D2B48C")),
    column(4,
    selectInput("label.color", "Label Color", multiple=FALSE, choices=c("white","black"), selected="black"),
    selectInput("theme","Select Theme",choices=c("Default","b/w"),selected="Default")),
    p("Welcome to APOSTL! Please set your cutoffs on the sidebar and your graph will be generated above."), strong("Please see below for a description:"),br(),
    strong("
    1) Saint Score Cutoff:"),p("
      - Define the cutoff for interaction confidence as specified by the SaintScore"),
    strong("
    2) Log2(Fold Change) Cutoff:"),p("
      - Define the cutoff for fold enrichment over your control purifications"),
    strong("
    3) Axes:"),p("
      - You are able to specify both your x and y axes with a number of parameters calculated either by SAINTexpress or APOSTL"),
    strong("
    4) Size:"),p("
      - These same options can be used to specify the size scaling of the bubbles within the graph"),
    strong("
    5) Color:"),p("
      - Define how the bubbles should be colored"),p("
      - 'crapome' requires you to define a 'crapome file' from", a("Workflow 1",href= "http://www.crapome.org/?q=wk_1_1_search"), "on the", a("crapome website.",href=
      "www.crapome.org")),
    strong("
    6)  Exclusion List:"), p("- Select proteins you wish to exclude from the graph space and network. Excluded proteins should be reported when publishing"),
    strong("7) Bubble Plot File Type:"), p("- Choose file type to download 300 dpi bubble graph (.tif,.png, or.pdf)")
    ),
    tabPanel("Cytoscape Network",
    rcytoscapejsOutput("network",width="125%",height="600px"), 
    p("Above is a cytoscape network of nodes/edges passing the filtering criteria. You are able to choose a 
      large amount of colors and some standard shapes."),p("BETA Release: The network refreshes itself and DOES NOT save user 
      layouts when altering parameters. Make sure to finalize filtering criteria and aesthetics before reorganizing 
      the network. There remains some instability in the node shapes with a large number of nodes present. 
      Try multiple shapes to resolve the issue.
      "),
    column(4,
    selectInput("node.color", "Select color for preys", multiple=FALSE, choices=colors, selected="#3B444B"),
    selectInput("bait.color", "Select color for baits", multiple=FALSE, choices=colors, selected="#8B0000")
    ),
    column(4,
    selectInput("node.shape", "Select shape for preys", multiple=FALSE, choices=shapes, selected="ellipse"),
    selectInput("bait.shape", "Select shape for baits", multiple=FALSE, choices=shapes, selected="ellipse")
    ),
    column(4,
    selectInput("edge.color", "Select color for edges", multiple=FALSE, choices=colors, selected="#C08081"),
    selectInput("net.layout", "Select layout algorithm", multiple=FALSE, choices=layouts, selected="cose")
    )
    ),tags$head(tags$script(src="cyjs.js")),
    tabPanel("Data Table",
             dataTableOutput('table'),
             downloadButton("saveTable", "Save Table")),
    tabPanel("Pathway Analysis",
             plotOutput("pathPlot",width="100%",height="500px"),
             column(4,
             selectInput("path_org","Organism", choices=c("mouse","yeast","human"),multiple=FALSE,selected="human"),
             sliderInput("path_pval", "pValue Cutoff",min=0,max=1, value=0.05)),
             column(4,
             selectInput("path_adj","pAdjustMethod", 
                         choices=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                         multiple=FALSE, selected="bonferroni")),
             column(4,
             selectInput("path_x","select x-axis",choices=c("pvalue","p.adjust")),
             downloadButton("pathTable", "Download Data"))),
    tabPanel("Gene Ontology",
             plotOutput("ontPlot",width="100%",height="500px"),
             column(4,
             selectInput("path_org","Organism", choices=c("mouse","yeast","human"),multiple=FALSE,selected="human"),
             sliderInput("path_pval", "pValue Cutoff",min=0,max=1, value=0.05)),
             column(4,
             selectInput("path_adj","pAdjustMethod", 
                         choices=c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                         multiple=FALSE, selected="bonferroni"),
             downloadButton("ontTable", "Download Data")),
             column(4,
             selectInput("path_x","select x-axis",choices=c("pvalue","p.adjust")),
             selectInput("GO_ont","Select Ontology", choices=c("MF","BP","CC"))))
    ))
))