# Winnie Huang and Stephanie Morrison - 04/23/2020
# Script to automate reading and analysis of WormLab .csv files
# V2.0: Take user input instead of mucking around in program for each use case
# Does both Thrashing and Crawling analysis
# 4/29/2020 Update: read_summary_data will ignore rows in .csv files past
# "Pirouette Time (s)" - fixes error from "mobility" stat rows just added
# by Wormlab team. 

################ INITIALIZATION ################################################

# clear R memory/variables
#rm(list = ls())

# check if packages are installed; if not, downloads them
if (!require(plyr)){
  install.packages("plyr")
} 
if (!require(ggplot2)){
  install.packages("ggplot2")
} 
if (!require(stringr)){
  install.packages("stringr")
} 

# loads packages
library(plyr)
library(ggplot2)
library(stringr)

# dir = "C:/Users/Hobo Hut/Desktop/R Worms/Crawling/unc-116/Crawling-Renamed"

# List of parameters used for Crawling analysis
parameters = c("Mean.Worm.Length..um.", "Mean.Width..um.", "Mean.Area..um.2.",
               "Center.point.Speed..um.s...Center.point.trajectory.time.",
               "Wavelength..um.", "Mean.Amplitude..um.", "Max.Amplitude..um.")

# Below for "old version" of Worm Lab program?
# parameters = c("Mean.Worm.Length..um.", "Mean.Width..um.", "Mean.Area..um.2.",
#                "Speed..um.s...track.length.time.",
#                "Wavelength..um.", "Mean.Amplitude..um.", "Max.Amplitude..um.")

# parameters = list("Thrash.min")

# Printed on startup - update as necessary
intro <- str_c(
  "************************************************************************\n",
  "* Worm Lab Analysis Program -                                          *\n",
  "* Winnie Huang and Stephanie Morrison, April 2020                      *\n",
  "************************************************************************\n",
  "*       *** Please see README.pdf for detailed instructions. ***       *\n",
  "************************************************************************\n",
  "* All .csv files should be in a single folder with no other documents, *\n",
  "* and should be all thrashing or all crawling data -                   *\n",
  "* no mixes of the two.                                                 *\n",
  "* Files should be in the format: 'STRAIN_TRIAL_DATE.csv'               *\n",
  "************************************************************************\n"
)
#############      FUNCTIONS     ###############################################

# Makes a data frame from the .csv file that matches "filename"
# If further updates to WormLab analysis software break this, look to see if the
# "skip=4" parameter should be changed to more or fewer lines of .csv data
read_summary_data <- function(filename){
  # Note: the original .csv files have a blank line at row 3. If the .csv is
  # saved using a spreadsheet program (e.g. Excel), this row is filled with 
  # empty cells, breaking the previous iteration of this program.
  # Skipping to line 4 treats blank rows and rows of blank cells the same.
  # nrows = 22 ignores any rows past "Pirouette Time (s)"
  file_data <- na.omit(data.frame(t(read.csv(filename, skip = 4, nrows = 22, 
                                             row.names = 1,
                                             header = FALSE))))
  
  # rename rows to fix gaps/errors
  rownames(file_data) <- 1:nrow(file_data)
  
  # extract vectors of thrash count and duration
  number_of_thrashes <- as.numeric(as.vector(file_data$Turn.Count))
  track_duration <- as.numeric(as.vector(file_data$Track.duration..s.))
  
  # Putting this data back in as a numeric prevents the graphs from breaking
  file_data$Duration.s <- track_duration
  
  # calculate thrashes per minute - used in Thrashing analysis
  file_data$Thrash.min <- number_of_thrashes/track_duration*60
  
  # ID every line of the data frame w/ filename, strain, date, trial
  ID_list <- strsplit(filename, "_")
  file_data$Worm.ID <- rep(filename, nrow(file_data))
  date <- strsplit(ID_list[[1]][3], '.csv')
  file_data$Date <- rep(date[[1]], nrow(file_data))
  file_data$Strain <- rep(ID_list[[1]][1], nrow(file_data))
  file_data$Trial <- rep(ID_list[[1]][2], nrow(file_data))
  
  return(file_data)
}


# Multiplot function to print multiple graphs on a single .pdf (or other dev)
# ggplot objects can be passed in ..., or to plotlist 
#                                               (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout matrix is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the layout matrix
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    # byrow = TRUE means that things will fill left to right, top to bottom
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow = TRUE)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    # Makes a viewpoint (basically a grid over a page) for each item in the
    # layout matrix
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Print each plot in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      # matchidx gets row and col for the item in "layout" that == i, thus
      # indexing through layout one-by-one - gets multiple points if one
      # graph is to go over more than one section
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      # prints plot i at coordinates designated by above
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# When given a df, list of strains, and a stat, make a df for each strain and
# run aov separately on each strain by stat ~ Date
# gather P-value and add to list of P-values in desired strain order
anova_p <- function(alldata, strainlevels, stat){
  p_values <- vector()  # Initialize so we can add in the for loop
  
  for (strain in strainlevels){
    # Make a df with only one strain's data
    strain.df <- alldata[which(alldata$strainsort == strain), ]
    
    # Cast the data in question into a double
    # df[[stat]] is the programmatic version of df$stat (which doesn't work)
    strain.df[[stat]] <- as.double(as.character(strain.df[[stat]]))
    
    # Left (stat) has to be a DOUBLE, not a numeric factor!
    # the formula/string stuff allows a variable to be used here.
    date.aov <- aov(as.formula(str_c(stat, " ~ Date")), strain.df)
    
    # Pull P-value from the ANOVA object
    p <- summary(date.aov)[[1]][["Pr(>F)"]][[1]]
    # Add p-value to the list of p_values
    p_values <- c(p_values, p)
  }
  
  # List of p-values; ordered by strain in same sequence as everywhere else
  return(p_values)
}


# GENERIC FUNCTION TO CALCULATE STAT BLOCK
# Calculates n, Mean, SD, SEM, and p-values for all strains
# Pass in the data frame with all data as well as the ordered list
# of *display* strain names AND *the stat to be assessed*
stat_calc <- function(alldata, strainlevels, stat){
  # Inner functions:
  
  # Calculates standard error of the mean;
  # copied wholesale from original program
  se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
  
  # Calculates p_value given two vectors
  p_value <- function(x,y) t.test(x,y)$p.value
  
  
  # list of vectors of all "stat" data in preferred strain order - using
  # display lists of strains in both df and levels
  stat_list <- lapply(strainlevels, function(levels)
                      as.numeric(subset(alldata, 
                                        alldata$strainsort == levels)[[stat]]))
  
  # vectors of calculated values for all strains - lengh, mean, sd, sem
  n_stat <- sapply(stat_list, length)
  mean_stat <- sapply(stat_list, mean)
  sd_stat <- sapply(stat_list, sd)
  SEM_stat <- sapply(stat_list, se)
  param = rep(stat, length(stat_list))
  
  
  # name rows for the p-value matrix - anonymous function call allows string
  # concatenation to be applied in a vectorized manner!
  p_column <- lapply(strainlevels, function(level) str_c("P value vs. ", level
                                                         , " (t-test)"))
  
  # p_matrix will hold all the p values of each strain vs. each other strain
  p_matrix <- matrix(nrow = length(stat_list), ncol = length(stat_list))
  
  # columns in matrix are named - this follows the data into the final result
  colnames(p_matrix) <- p_column
  
  # calculate all p-values by iterating through all strain:strain possibilites
  # Might be able to avoid the for() loop here, but it's fine for now.
  for (i in 1:length(stat_list)){
    for (j in 1:length(stat_list)){
      p_matrix[i,j] <- p_value(stat_list[[i]], stat_list[[j]])
    }
  }

  # Get vector of date-comparison ANOVA p-values for each strain
  anova.p <- anova_p(alldata, strainlevels, stat)
  
  # make final data frame
  results <- data.frame("Strain" = strainlevels, "Parameter" = param,
                        "n" = n_stat,
                        "Mean" = mean_stat, "SD" = sd_stat,
                        "SEM" = SEM_stat)
  
  # adds p-value matrix onto right side of results
  results <- cbind(results, p_matrix)
  
  results$ANOVA <- anova.p
  
  anovacolname <- "Day-to-Day Comparison ANOVA P-Value (Within Strain)"
  
  names(results)[length(names(results))] <- anovacolname
  
  #explicit return for clarity
  return(results)
}


# Writes user-provided strain info (order and preferred display name) into two
# .txt files for use in future script runs to skip user input steps
write_strains <- function(ordered, display){
  capture.output(writeLines(ordered), file = "#OrderedStrains.txt")
  capture.output(writeLines(display), file = "#DisplayStrainNames.txt")
}


# Reads .txt file that was written in previous script run and returns a 
# character vector for use in skipping user input.
# *Only reads one file at a time* in contrast to write_strains()
read_strain <- function(filename){
  capture <- scan(file = filename, what = character(), sep = "\n")
  return(capture)
}


##########    USER INPUT     ###################################################
writeLines(intro)


# sets working directory to address provided by user
setwd(readline(prompt="Copy/Paste Address of Folder Containing .csv Files: "))

# # For debug purposes
# setwd(dir)


## STEP ONE: Get file names, split into date, strain, trial
# makes character vector with names of all .csv files in directory
data_files <- tolower(dir(pattern = ".csv"))

# makes matrix with all filenames split by "_"
# add a check here to throw an error if there are more/fewer than three columns?
filename_split <- str_split(data_files, "_", simplify = TRUE)

if (length(filename_split) %% 3 != 0){
  writeLines(str_c("\n*** Warning: Non-standard file names detected. ***\n",
            "Recommend stopping analysis and re-checking file names.\n"))
}

# gets factor of all the strain names, which we are assuming is in column 1
# of the matrix, i.e. the first substring in all the filenames
strain_names <- as.factor(filename_split[,1])


## STEP TWO: Gather info from User about strain order, names, experiment name,
## type of analysis to be done

# This looks for any files in dir that start with # and end with txt
strain_files <- dir(pattern = c("^#", "txt$"))

userinput <- TRUE

# Check if the vector is populated with the correct number
# If yes, ask the user if they want to use that data
if (length(strain_files) == 2){
  inputloop <- TRUE
  # Gotta sanitize that input.
  while(inputloop){
    writeLines(str_c("Strain Name information found - Would you like to use ",
                       "saved strain names for this analysis?"))
    answer <- tolower(readline(prompt="Y or N: "))
    if (answer %in% c("y", "yes")){
      # I know, hardcoding the filenames. Sue me!
      strains_ordered <- read_strain("#OrderedStrains.txt")
      display_strain_names <- read_strain("#DisplayStrainNames.txt")
      inputloop <- FALSE
      userinput <- FALSE
    } else if (answer %in% c("n", "no")){
      inputloop <- FALSE
    }
  }
}

if (userinput) {  # Only called if no .txt files exist or user not using files
  # Ask the user to assign an order to the strains - important that this is
  # using the original strain names for use in the stats later
  
  # list of all strains to be selected from; depletes as selections continue
  strain_menu <- levels(strain_names)
  
  # the list generated by the following for/menu loop
  strains_ordered <- character(length(strain_menu))
  
  # Asks the user to select the desired order of strains; removes already
  # selected options from the menu
  for(i in 1:nlevels(strain_names)){
    menu_title <- str_c("Select the strain to be displayed in slot ", i, " of ",
                        nlevels(strain_names), "\n('0' will exit program)")
    selection <- menu(strain_menu,title = menu_title)
    strains_ordered[i] <- strain_menu[selection]
    strain_menu <- strain_menu[-selection]
  }
  
  # Creates a character vector and populates it with user-supplied
  # display names in the order designated by the user in previous for loop
  display_strain_names <- character(nlevels(strain_names))
  for(i in 1:nlevels(strain_names)){
    display_strain_names[i] <- readline(prompt =
                          str_c("Enter the desired display name for strain '",
                                                strains_ordered[i], "': " ))
  }
  # Ask user if they want to save the strain name data.
  inputloop = TRUE
  # Gotta sanitize that input
  while(inputloop){
    writeLines(str_c("Would you like to save the strain name data ",
                      "for future analysis?"))
    answer <- tolower(readline(prompt="Y or N: "))
    if (answer %in% c("y", "yes")){
      write_strains(strains_ordered, display_strain_names)
      inputloop = FALSE
    } else if (answer %in% c("n", "no")){
      inputloop = FALSE
    }
  }
}


# Ask user to name the experiment
writeLines(str_c("\nWhat is the name of this experiment?", 
                 "\n(This will be appended to each result",
                 " document generated, so make it short.)"))
experiment <- readline(prompt="Experiment Name: ")


## Ask user if this is thrashing or crawling data
menu_title <- str_c("Are these files Crawling or Thrashing data?")
analysis_choice <- menu(c("Crawling", "Thrashing"),title = menu_title)
if (analysis_choice == 2){
  experiment <- str_c(experiment, "_thrashing")
  parameters <- list("Thrash.min")
} else {
  experiment <- str_c(experiment, "_crawling")
}


## STEP 3: 
############ END USER INPUT - BEGIN PROCESSING DATA + OUTPUT ###################

# Initialize all_data variable
all_data <- data.frame()
# Populate all_data with all information from the .csv files
for(i in 1:length(data_files))
{
  all_data <- rbind(all_data, read_summary_data(data_files[i]))
}


# Create Analysis directory and move wd into it for creating files
dir.create(str_c(experiment, "_analysis"))
setwd(str_c(experiment, "_analysis"))


# The following block allows all_data to be indexed by display_strain_name:
# Add factor to all_data with preferred strain order
all_data$strainsort <- factor(all_data$Strain, levels = strains_ordered)
# Data is sorted by desired order
all_data <- arrange(all_data, strainsort)
# Original strain names are replaced with "Display" names from user
levels(all_data$strainsort) <- display_strain_names
  

# SPLIT HERE FOR CRAWLING VS THRASHING
if (analysis_choice == 1){
  # CRAWLING ANALYSIS
  # Make a data frame with all stats calculated
  # "parameters" is defined at the top of this document and can be changed as 
  # desired to add an arbitrary number of stat calculations.
  all_stats <- data.frame()
  for (param in parameters){
    all_stats <- rbind(all_stats, stat_calc(all_data, 
                                            display_strain_names, param))
    all_stats[nrow(all_stats)+1,] <- NA  # Creates a blank row between stats
  }
  
  write.csv(all_stats, file = str_c(experiment, "_allstats.csv"), 
            row.names = FALSE, na="")
  
  write.csv(all_data, file = str_c(experiment, "_pooleddata.csv"), 
            row.names = FALSE)
  
  # Stats have to be changed into doubles to be used as continuous data
  all_data$Mean.Worm.Length..um. <- as.double(
    as.character(all_data$Mean.Worm.Length..um.))
  
  all_data$Center.point.Speed <- as.double(as.character(
    all_data$Center.point.Speed..um.s...Center.point.trajectory.time.))
  
  # Get ANOVA p values into all_data.
  # Doing this AFTER writing the .csvs because I am rounding the p-values and
  # don't want to present to the user confusing/misleading data.
  all_data$ANOVAlength <- all_data$strainsort
  p_vals <- round(anova_p(all_data, display_strain_names, 
                          "Mean.Worm.Length..um.") , digits = 6)
  levels(all_data$ANOVAlength) <- factor(p_vals)
  
  
  all_data$ANOVAspeed <- all_data$strainsort
  p_vals <- round(anova_p(all_data, display_strain_names, "Center.point.Speed")
                  , digits = 6)
  levels(all_data$ANOVAspeed) <- factor(p_vals)
  
  
  # Create graphs in ggplot2:
  # Main plots - mean worm length and speed / strain
  length_box <- ggplot(all_data, aes(strainsort, Mean.Worm.Length..um.)) +
    geom_boxplot() +
    geom_jitter(color = "orange", size = 1) +
    labs(title = "Mean Worm Length",
         x = "Strain Name", y = "Mean Worm Length (um)") +
    theme (axis.text.x = element_text(angle = 20, hjust = 1)) +
    expand_limits(y=0)
  
  length_hist <- ggplot(all_data, aes(Mean.Worm.Length..um.)) +
    geom_histogram(color = "black", fill = "white", bins = 30) +
    facet_grid(strainsort ~ .) +
    labs(title = "Mean Worm Length (Histogram)",
         x = "Mean Worm Length (um)", y = "Count")
  
  speed_box <- ggplot(all_data, aes(strainsort, Center.point.Speed)) +
    geom_boxplot() +
    geom_jitter(color = "orange", size = 1) +
    labs(title = "Mean Worm Speed",
         x = "Strain Name", y = "Mean Worm Speed (um/s)") +
    theme (axis.text.x = element_text(angle = 20, hjust = 1)) +
    expand_limits(y=0)
  
  speed_hist <- ggplot(all_data, aes(Center.point.Speed)) +
    geom_histogram(color = "black", fill = "white", bins = 30) +
    facet_grid(strainsort ~ .) +
    labs(title = "Mean Worm Speed (Histogram)",
         x = "Mean Worm Speed (um/s)", y = "Count")
  
  plots <- list(length_box, length_hist, speed_box, speed_hist)
  
  pdf(str_c(experiment, "_plots.pdf"))
  multiplot(plotlist = plots, cols = 2)
  dev.off()
  
  # day-to-day comparison plots - first is the same box/scatter by day and
  # strain, second is a histogram of distributions per strain per day
  # day_box won't break when strains have different days than each other,
  # but the histogram could get weird in that case.
  # ANOVA not added for Crawling because there are multiple plots
  lengthday_box <- ggplot(all_data, aes(Date, Mean.Worm.Length..um.)) +
    geom_boxplot(aes(fill = strainsort)) +
    geom_jitter(shape = 1) +
    facet_wrap(vars(strainsort), scales = "free_x") +
    labs(title = "Mean Worm Length by Day of Assay",
         x = "Date of Assay", y = "Mean Worm Length") +
    theme (axis.text.x = element_text(angle = 30, hjust = 1),
           legend.position = "none") +
    geom_text(aes(y = max(Mean.Worm.Length..um.)*1.05, x=0.5, 
                  label = ANOVAlength), 
              data = all_data, hjust = "inward")
  
  lengthday_hist <- ggplot(all_data, aes(Mean.Worm.Length..um.)) +
    geom_histogram(aes(fill = strainsort), color = "black", bins=20) +
    facet_grid(strainsort ~ Date) +
    labs(title = "Mean Worm Length by Day of Assay - Histogram",
         x = "Mean Worm Length", y = "Count") +
    theme(legend.position = "none")
  
  speedday_box <- ggplot(all_data, aes(Date, Center.point.Speed)) +
    geom_boxplot(aes(fill = strainsort)) +
    geom_jitter(shape = 1) +
    facet_wrap(vars(strainsort), scales = "free_x") +
    labs(title = "Mean Worm Speed by Day of Assay",
         x = "Date of Assay", y = "Mean Worm Speed") +
    theme (axis.text.x = element_text(angle = 30, hjust = 1),
           legend.position = "none") +
    geom_text(aes(y = max(Center.point.Speed)*1.05, x=0.5, label = ANOVAspeed), 
              data = all_data, hjust = "inward")
  
  speedday_hist <- ggplot(all_data, aes(Center.point.Speed)) +
    geom_histogram(aes(fill = strainsort), color = "black", bins=20) +
    facet_grid(strainsort ~ Date) +
    labs(title = "Mean Worm Speed by Day of Assay - Histogram",
         x = "Mean Worm Speed", y = "Count") +
    theme(legend.position = "none")
  
  day_list <- list(lengthday_box, lengthday_hist, speedday_box, speedday_hist)
  
  pdf(str_c(experiment, "_daycomparisonplots.pdf"))
  multiplot(day_list)
  dev.off()
  # End of Crawling analysis and output
} else if (analysis_choice == 2){
  
  # DO THRASHING
  # Stats are calculated and written to .csv
  thrash_stats <- stat_calc(all_data, display_strain_names, "Thrash.min")
  write.csv(thrash_stats, file = str_c(experiment, "_stats.csv"), 
            row.names = FALSE, na = "")
  
  write.csv(all_data, file = str_c(experiment, "_pooleddata.csv"))
  
  # Get ANOVA p values into all_data. 
  # Doing this AFTER writing the .csvs because I am rounding the p-values and
  # don't want to present to the user confusing/misleading data.
  all_data$ANOVAp <- all_data$strainsort
  p_vals <- round(anova_p(all_data, display_strain_names, "Thrash.min")
                  , digits = 6)
  levels(all_data$ANOVAp) <- factor(p_vals)

  
  # Create graphs in ggplot2:
  
  # Main plots - thrash.min by strain
  thrash_plot1 <- ggplot(all_data, aes(strainsort, Thrash.min)) +
    geom_boxplot() +
    geom_jitter(color = "orange", size = 1) +
    labs(title = "Thrashes per Minute",
         x = "Strain Name", y = "# Thrashes / min") +
    theme (axis.text.x = element_text(angle = 20, hjust = 1)) +
    expand_limits(y=0)
  
  thrash_plot2 <- ggplot(all_data, aes(Thrash.min)) +
    geom_histogram(color = "black", fill = "white", bins = 30) +
    facet_grid(strainsort ~ .) +
    labs(title = "Thrashes per Minute (Histogram)",
         x = "# Thrashes per minute", y = "# of Worms")
  
  plots <- list(thrash_plot1, thrash_plot2)
  
  pdf(str_c(experiment, "_plots.pdf"))
  multiplot(plotlist = plots, cols = 2)
  dev.off()
  
  # day-to-day comparison plot - first is the same box/scatter by day and
  # strain, second is a histogram of thrashing distributions per strain per day
  # day_box won't break when strains have different days than each other, 
  # but the histogram could get weird in that case.
  day_box <- ggplot(all_data, aes(Date, Thrash.min)) +
    geom_boxplot(aes(fill = strainsort)) +
    geom_jitter(shape = 1) +
    facet_wrap(vars(strainsort), scales = "free_x") +
    labs(title = "Thrashes per Minute by Day of Assay", 
         x = "Date of Assay", y = "# Thrashes per Minute") +
    theme (axis.text.x = element_text(angle = 30, hjust = 1),
           legend.position = "none") +
    # This places the ANOVA p-value on top of each chart at the y-value of the
    # maxiumum value. Might cover some data.
    # Uncover data by commenting out or adding some distance to y = 
    geom_text(aes(y = max(Thrash.min)*1.05, x=0.5, label = ANOVAp), 
              data = all_data, hjust = "inward")
  
  day_plot <- ggplot(all_data, aes(Thrash.min)) +
    geom_histogram(aes(fill = strainsort), color = "black", bins=20) +
    facet_grid(strainsort ~ Date) +
    labs(title = "Thrashes per Minute by Day of Assay - Histogram",
         x = "# Thrashes per Minute", y = "# of Worms") +
    theme(legend.position = "none")
  day_list <- list(day_box, day_plot)
  
  pdf(str_c(experiment, "_daycomparisonplots.pdf"))
  multiplot(day_list)
  dev.off()
  # End of thrashing analysis and output
}
