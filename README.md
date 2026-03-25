# ShinyHaplotypes

This project using a sliding window approach to compare haplotype similarity among various fungal isolates on a chromosome-by-chromosome basis. Plot lines are colored according to phylogenetic groupings, which allows one to visualize admixed chromosome regions. Input files are in the format used by Chromopainter (https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromopainter_info.html)
Perl scripts are used to covert the chromopainter files into datasets that can be imported into a R Shiny app which allows one to visualize haplotype similarity between isolates in an interactive manner.

![ShinyHaplotypes](/ShinyHaplotypes.jpg)

## Getting Started

The scripts can be cloned or uploaded to a local directory. 

### Prerequisites

Perl 5.0 or greater; RSTUDIO; R packages: data.table, dplyr, graphics, Shiny, 

### Installation

Download the ZIP archive to your local machine. Extract the archive:

```
unzip ShinyHaplotypes-master.zip
```

Change to the ShinyHaplotypes-master directory:

```
cd ShinyHaplotypes-master
```

Now you're ready to start using the tools.

## Generate the ShinyHaplotypes data files:
Usage: perl ShinyData.pl <Haplotypes_file> <window-size> <step-size> <output-directory> <scale? (default = no)>
```
perl Slide_compare.pl chr1_haplotypes.txt 1000 200 SHINYOUT no
```

### Running the SHINY app

1. Open the ShinyPlot.R script in RStudio
2. Click "Run App"
3. Select the directory containing the datafiles
4. Start exploring haplotype diversity by selecting different toggles and sliders

### Analysis of specific chromosome regions/strains

1. Identify a region with an interesting pattern
2. Use the mouse to click-drag to select any number of plot traces. This will capture information about the chromosome region specified, as well as the strains whose plot lines are interested by the box boundaries
3. Use the pop-up menu to select the type of analysis required (list of capture strains; pairwise divergence analysis; phylogenetic tree construction; PCA analysis)
4. Use buttons at bottom of rendered plot to export graphics, or underlying dataset for analysis using third party programs

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

