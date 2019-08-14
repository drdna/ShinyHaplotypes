# ShinyHaplotypes

This project using a sliding window approach to compare haplotype similarity among various fungal isolates on a chromosome-by-chromosome basis. Plot lines are colored according to phylogenetic groupings, which allows one to visualize admixed chromosome regions. Input files are in the format used by Chromopainter (https://people.maths.bris.ac.uk/~madjl/finestructure-old/chromopainter_info.html)
Perl scripts are used to covert the chromopainter files into datasets that can be imported into a R Shiny app which allows one to visualize haplotype similarity between isolates in an interactive manner.   

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

## Running the programs

First we need to take metadata from the strain.idfile and add it to the haplotype information contained in the main chromopainter datafiles. The following commands show you how to run the program using the provided sample datasets. First we run the Convert_sliding4R.pl script. 

Usage: Convert_sliding4R.pl #idfile# #cpfile# #outfile-name#
```
perl PERL/Convert_sliding4R.pl strain.idfile HAPLOTYPES/chr1.cp chr1_haplotypes.txt
```
This can repeated using the datasets for all chromosomes.

Next, we run the sliding window analysis that compares haplotype similarity along the chromosomes. Window-size (number of variant sites compared) and window step size are user definable. Values of 1000 and 200, respectively produce plots with sufficient definition and navigability. Usage: Slide_compare.pl <outfile-from-Convert_sliding4R> <window-size> <step-size> <output-directory>

```
perl Slide_compare.pl chr1_haplotypes.txt 1000 200 SHINYOUT
```

### Running the SHINY app

Open the Haplotype_div_Shiny.R script in RStudio. If necessary change the path to the directory containing the .diffs datafiles. Save any changes and then click 'Run App.' This will open an interactive window, where you can select which strain to compare, with which host-specialized populations, and on which chromosome. Also available is the ability to zoom in on a specific regions of the haplotype.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

