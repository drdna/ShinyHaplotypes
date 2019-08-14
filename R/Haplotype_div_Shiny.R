'%!in%' <- function(x,y)!('%in%'(x,y))
require('ggplot2')
library(data.table)
library(shiny)
library(dplyr)
library(graphics)


# First block of code sets up initial window

# This can be any file in the directory containing the diffs files

# Change "Chr6" & "AR4" to match relevent fields from any filename in the input directory

df <- read.table(paste("~/SHINYSLIDE/", "Chr6", ".", "AR4", ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

strains <- row.names(df)

# Folowing line needs to be changed to include whatever sample is present in the above-specified filename

strains[[length(strains)+1]] <- "AR4"   #add reference strain to the strain list

sstrains <- sort(strains)

chromosomes <- as.character(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7"))

plotLength <- 600           # initial datapoints to plot

includedClades <- as.character(c(unique(df$clade)))

cladeNames <- c("Brachiaria1" = "B1",
                "Brachiaria2" = "B2",
                "Cynodon" = "C",
                "Echinochloa" = "Ec",
                "Eleusine1" = "E1",
                "Eleusine2" = "E2",
                "Eragrostis" = "Er",
                "Hakonechloa" = "Ha",
                "Lolium" = "L",
                "Lolium2" = "L2",
                "Leptochloa" = "Lep",
                "Leersia" = "Le",
                "Oryza" = "O",
                "Paspalum" = "Pa",
                "Panicum" = "P",
                "Setaria" = "S",
                "Stenotaphrum" = "St",
                "Triticum" = "T",
                "Unknown" = "U")

cladeNamesR <- c("B1" = "Brachiaria1",
                 "B2" = "Brachiaria2",
                 "C" = "Cynodon",
                 "Ec" = "Echinochloa",
                 "E1" = "Eleusine1",
                 "E2" = "Eleusine2",
                 "Er" = "Eragrostis",
                 "Ha" =  "Hakonechloa",
                 "L" = "Lolium",
                 "L2" = "Lolium2",
                 "Lep" = "Leptochloa",
                 "Le" = "Leersia",
                 "O" = "Oryza",
                 "Pa" = "Paspalum",
                 "P" = "Panicum",
                 "S" = "Setaria",
                 "St" = "Stenotaphrum",
                 "T" = "Triticum",
                 "U" = "Unknown")


# perform filtering to stop missing clades from showing in display window

filtered <- as.vector(cladeNamesR[includedClades])

cladeNamesFiltered <- sort(cladeNames[filtered])


ui <- fluidPage(

  sidebarLayout(
    sidebarPanel(
      radioButtons(inputId = "inputStrain", "Select strain:", 
                   inline = TRUE, choices = sstrains,
                   selected = sstrains[[1]]
      ),
      
      checkboxGroupInput("inputClades", "Plot clade(s):", 
                         inline = TRUE, choices = cladeNamesFiltered,
                         selected = cladeNamesFiltered[[1]]
      ),
      
      radioButtons(inputId = "inputChromo", "Chromosome to plot:",
                         inline = TRUE, choices = chromosomes,
                         selected = chromosomes[[7]]
      ),

      sliderInput("region", "Select haplotype range:",
                  min = 0, max = 1500, value = c(1,1000)
      )

    ),

    mainPanel(
      plotOutput("Windows")
    )
  )
)

server <- function(input, output) {
  
  output$Windows <- renderPlot({
    
    # Read in haplotype information 
        
    df <- read.table(paste("~/SHINYSLIDE/", input$inputChromo, ".", input$inputStrain, ".diffs", sep = ""), header = TRUE, row.names = 1, check.names = FALSE)
    
    sstrains <- row.names(df)
    
    # Get sliding window information
    
    positionVector <- df[1:ncol(df)-1]
    
    SNPsites <- colnames(positionVector)
    
    revWindowEdges <- rev(SNPsites)
    
    lastWindow <- as.numeric(revWindowEdges[[1]])
    
    haplOnly <- filter(df, clade %in% input$inputClades)
    
    haplPlot <- haplOnly[ ,-length(haplOnly)]
    
    tdf <- t(haplPlot)

    palette(c("darkolivegreen2", "darkorange", "darkorchid1", "darkseagreen", "darkslategray1", "deeppink", "deepskyblue", "dodgerblue", "firebrick", "firebrick1", "gold", "forestgreen", "slateblue", "snow3", "sienna4", "black", "snow4"))

    
    matplot(tdf, xlim = input$region, ylim = c(0,0.2), type = "l", lty = 1, lwd = 2, main = "Haplotype Divergence",
            col = haplOnly$clade, ylab ='Divergence (%)', xlab = 'Haplotype position')
            legend("top", legend = unique(haplOnly$clade), lty=c(1,1), cex=0.8, lwd=10, col = unique(haplOnly$clade), horiz = TRUE)
    
 
  })
  
}

shinyApp(ui = ui, server = server)

