#' Generate a Manhattann Plot from a raw MatrixEQTL output file
#'
#' @param matrixeqtl.output.file The assumption is that MatrixEQTL was run with a single phenotype, and produced a tab delimited output file with the following columns:
#' SNP - character column, an example of data in this column is "chr8-81050409"
#' p-value - numeric column, an example of data in this column is 0.00005.
#' Please provide the file path. 
#' @param genome.wide.sig.cutoff A horizontal red line will be drawn on the output to indicate the p-value cut-off for genome-wide significance. The default value is 5e-08.
#' @param return.or.save This argument instructs the function to either save the plot to output.file, or to return the plot using `return()`. 
#' @param output.file This argument is a character string to indicate the file in which to save the output file, given that `return.or.save` is set to "save".
#' Please indicate /path/to/file_name without a file extension. The plot will be saved as a 200dpi .png file. 
#'
#' @examples
#' 
#' # return to object 
#' my.plot = plot.manhattan.from.matrixeqtl.output.file(
#' matrixeqtl.output.file="/path/to/matrix_eqtl_outfile.tsv", 
#' genome.wide.sig.cutoff=5e-08, 
#' return.or.save="return"
#' )
#' 
#' # view plot
#' my.plot
#' 
#' # save to file
#' plot.manhattan.from.matrixeqtl.output.file(
#' matrixeqtl.output.file="/path/to/matrix_eqtl_outfile.tsv", 
#' genome.wide.sig.cutoff=5e-08, 
#' return.or.save="save", 
#' output.file = "path/to/output_plot_file"
#' )

plot.manhattan.from.matrixeqtl.output.file = function(matrixeqtl.output.file="", genome.wide.sig.cutoff=5e-08, return.or.save="save", output.file=""){
  
  # test user input
  if(length(matrixeqtl.output.file)==0){
    stop("Please provide a MatrixEQTL input file")
  }
  if( (return.or.save=="save") & (length(output.file)==0) ){
    stop("Please provide the argument outfile.name if you wish to save the output plot. Otherwise, set return.or.save as 'print' to save the plot to an object or print without saving.")
  }
  
  # read in file
  matrixeqtl.output.df = readr::read_delim(matrixeqtl.output.file, delim="\t", show_col_types = FALSE)
  
  # test user input
  if( !("SNP" %in% colnames(matrixeqtl.output.df)) | !("p-value") %in% colnames(matrixeqtl.output.df) ){
    stop("Missing columns in input file, please ensure your dataframe has at least a SNP and p-value column.")
  }
  
  # test user input
  if(dim(matrixeqtl.output.df)[1]==0){
    stop("Please provide non-empty input file")
  }
  
  # pre-process input df
  matrixeqtl.output.df = matrixeqtl.output.df %>% 
    # generate CHR and BP cols from SNP  
    tidyr::separate(col=SNP, into=c("CHR", "BP"), sep="-") %>% 
    dplyr::mutate(CHR = gsub("chr","", CHR)) %>% 
    # ensure X and Y chr labels are encoded as numbers
    dplyr::mutate(CHR = replace(CHR, CHR == "X", 23)) %>% 
    dplyr::mutate(CHR = replace(CHR, CHR == "Y", 24)) %>% 
    # make CHR and BP numeric
    dplyr::mutate(CHR = as.numeric(CHR)) %>% 
    dplyr::mutate(BP = as.numeric(BP)) %>% 
    # rename for compatibility 
    dplyr::rename(
      P=`p-value` 
    ) 
  
  # Add cols that will enable correct spacing of chromosomes along the x-axis 
  matrixeqtl.output.df = matrixeqtl.output.df %>% 
    
    # Calculate size of each chromosome
    dplyr::group_by(CHR) %>% 
    dplyr::summarise(chr_length=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    dplyr::mutate(total=cumsum(chr_length)-chr_length) %>%
    dplyr::select(-chr_length) %>%
    
    # Bind back to results dataframe
    dplyr::left_join(matrixeqtl.output.df, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    tidyr::drop_na(CHR) %>% 
    tidyr::drop_na(BP) %>% 
    dplyr::arrange(CHR, BP) %>% 
    dplyr::mutate(BP_cumulative=BP+total)
  
  # Create x-axis spacing df 
  x_axis_spacing = matrixeqtl.output.df %>% 
    dplyr::group_by(CHR) %>% 
    dplyr::summarize(center=(max(BP_cumulative) + min(BP_cumulative))/2)
  
  # Manhattan plot
  manhattan.plot = ggplot(matrixeqtl.output.df, aes(x=BP_cumulative, y=-log10(P))) +
    
    # Show all points
    geom_point(aes(color=as.factor(CHR)), alpha=0.75, size=0.8) +
    scale_color_manual(values = rep(c("grey", "black"), 22)) +
    
    # custom X axis:
    scale_x_continuous(label = x_axis_spacing$CHR, breaks = x_axis_spacing$center) +
    scale_y_continuous(limits = c(0, 15), expand = c(0, 0)) +     # remove space between plot area and x axis
    
    # Add line to indicate genome wide sig.
    geom_hline(yintercept = -log10(5e-8), colour="red", alpha=0.5)  + 
    xlab("Genomic position") +
    
    # Aesthetics
    theme_bw(base_size = 8) +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  
  if(return.or.save=="return"){return(manhattan.plot)}
  if(return.or.save=="save"){
    
    # define function to save the plot as a png
    save.plot <- function(plot.to.save, output.name) {
      png(paste(output.name), width = 5, height = 3, units="in", res=200)
      print(plot.to.save)
      dev.off() 
    }
    
    # save plot wit user-provided outfile.name
    save.plot(plot.to.save=manhattan.plot, output.name=paste0(output.file, ".png"))
    
  }
}


