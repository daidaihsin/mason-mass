library(dplyr)
library(foreach)
library(doParallel)
library(ggplot2)
library(gridExtra)

#[TODO] specify log file path
parallel_log = ""

##### IMPORT DATA #####
rows_csv_colname = c("entry", "nParticle", "pwflag", "charge", "pt", "eta", "phi", "px", "py", "pz")
other_csv_colname = c("value")
csv_col_classes = c(rep("numeric",4), rep("NULL",3), rep("numeric",3))
  
import_data_source <- function(working_dir, raw_csv_filename, compressed_csv_filename, exist_header=F, nrow_to_read=NA){
  cat('IMPORT START: ', strftime(Sys.time(),"%Y-%m-%d %H:%M:%S"), '\n')
  cat('length(working_dir): ', length(working_dir), '-- working_dir: ', working_dir, '\n')
  cat('length(raw_csv_filename): ', length(raw_csv_filename), '-- raw_csv_filename:', raw_csv_filename, '\n')
  
  if(working_dir != ''){
    setwd(working_dir)
  }
  
  if(raw_csv_filename != '') {
    ifelse(is.na(nrow_to_read),
           collision_data <- read.csv(raw_csv_filename, header = exist_header),
           collision_data <- read.csv(raw_csv_filename, header = exist_header, nrows = nrow_to_read))
    
    colnames(collision_data) <- rows_csv_colname
  }else{
    compressed_csv_filename <- gsub("\\.csv$", "", compressed_csv_filename)
    
    filenames <- gsub("\\.csv$", "", list.files(pattern="\\.csv$"))
    for(i in filenames){
      ifelse(is.na(nrow_to_read),
             ifelse(i == compressed_csv_filename, assign(i, read.csv(paste(i, ".csv", sep=""), header=exist_header, col.names = rows_csv_colname, colClasses = csv_col_classes), envir = .GlobalEnv),
                                                  assign(i, read.csv(paste(i, ".csv", sep=""), header=exist_header, col.names = other_csv_colname), envir = .GlobalEnv)),
             ifelse(i == compressed_csv_filename, assign(i, read.csv(paste(i, ".csv", sep=""), header=exist_header, col.names = rows_csv_colname, nrows = nrow_to_read, colClasses = csv_col_classes), envir = .GlobalEnv),
                                                  assign(i, read.csv(paste(i, ".csv", sep=""), header=exist_header, col.names = other_csv_colname), envir = .GlobalEnv))
      )
   
    }
   
    #cat('min_group_px: ', min(ori_rows$px), '--min_group_py: ', min(ori_rows$py),'--min_group_pz: ', min(ori_rows$pz),'\n')
    #cat('max_group_px: ', max(ori_rows$px), '--max_group_py: ', max(ori_rows$py),'--max_group_pz: ', max(ori_rows$pz),'\n')
    
    collision_data <- get(compressed_csv_filename)
  }
  
  # take only charged particals to calculate meson mass
  collision_data <- collision_data[collision_data$pwflag==0 & collision_data$charge !=0, ]
  
  cat('nrow(collision_data): ', nrow(collision_data), '\n')
  cat('IMPORT FINISH: ', strftime(Sys.time(),"%Y-%m-%d %H:%M:%S"), '\n')
  return (collision_data)
}


decode_compressed_data <- function(compressed_data_df){
  cat('DECODING START: ', strftime(Sys.time(),"%Y-%m-%d %H:%M:%S"), '\n')
 
  n_cores <- detectCores() - 5
  cl = makeCluster(n_cores, outfile = parallel_log)
  for(col in c("px", "py", "pz")){
    #for(col in colnames(compressed_data_df)){
    cat('col: ', col, '\n')
    tmp_df = get(col)
    compressed_data_df[,col]<-parSapply(cl, compressed_data_df[,col], function(x){tmp_df[x+1,]})
    #compressed_data_df[,col]<-sapply(compressed_data_df[,col], function(x){tmp_df[x+1,]})
  }
  stopCluster(cl)
  cat('DECODING FINISH: ', strftime(Sys.time(),"%Y-%m-%d %H:%M:%S"), '\n')
  
  #check
  wrong_column <- sapply(compressed_data_df, function(x) sum(is.na(x)))
  if(any(wrong_column>0)) {
    msg <- paste0('DECODING ERROR WITH COLUMN: ', lapply(wrong_column[wrong_column>0], function(x) paste0(x, ", ")), '\n')
    stop(msg)
  }
  
  return (compressed_data_df)
}



##### CALCULATE MESON MASS #####
meson_mass_calculation_dict <- function(same_entry_df){
  # return mass
  mass_list <- c()

  n_charged_partical <- nrow(same_entry_df)
  n_charged_partical_to_pair <- nrow(same_entry_df) - 1
  
  # only one partical
  if(n_charged_partical_to_pair < 1){
    return (as.data.frame(mass_list))
  }
  
  # iterate two particals and calculate mass by using E^2 = P^2 + m^2
  for(i in 1 : n_charged_partical_to_pair){
    first_partical <- same_entry_df[i,]
    
    #dictionary mapping
    first_partical$px <- px[first_partical$px+1, ]
    first_partical$py <- py[first_partical$py+1, ]
    first_partical$pz <- pz[first_partical$pz+1, ]
    first_partical$momentum_square <- (first_partical$px)^2 + (first_partical$py)^2 + (first_partical$pz)^2
    
    for(j in (i+1) : n_charged_partical){
      second_partical <- same_entry_df[j,]
      
      #dictionary mapping
      second_partical$px <- px[second_partical$px+1, ]
      second_partical$py <- py[second_partical$py+1, ]
      second_partical$pz <- pz[second_partical$pz+1, ]
      
      second_partical$momentum_square <- (second_partical$px)^2 + (second_partical$py)^2 + (second_partical$pz)^2
      
      first_add_second <- data.frame(momentum_square = (first_partical$px + second_partical$px)^2 + (first_partical$py + second_partical$py)^2 + (first_partical$pz + second_partical$pz)^2,
                                     pt_square = (first_partical$px + second_partical$px)^2 + (first_partical$py + second_partical$py)^2)
      
      if(first_partical$charge != second_partical$charge){
        #case1: pion_kaon
        first_partical$energy <- sqrt(first_partical$momentum_square + 0.13957018^2)
        second_partical$energy <- sqrt(second_partical$momentum_square + 0.493677^2)
        mass <- sqrt((first_partical$energy + second_partical$energy)^2 - first_add_second$momentum_square)
        pt <- sqrt(first_add_second$pt_square)
        
        if(pt>5){
          mass_list <- append(mass_list, mass)
        }
        
        #case2: kaon_pion
        first_partical$energy <- sqrt(first_partical$momentum_square + 0.493677^2)
        second_partical$energy <- sqrt(second_partical$momentum_square + 0.13957018^2)
        mass <- sqrt((first_partical$energy + second_partical$energy)^2 - first_add_second$momentum_square)
        pt <- sqrt(first_add_second$pt_square)
        
        if(pt>5){
          mass_list <- append(mass_list, mass)
        }
      }
    }
  }
  
  return (as.data.frame(mass_list))
}


meson_mass_calculation <- function(same_entry_df){
  # return mass
  mass_list <- c()
  
  n_charged_partical <- nrow(same_entry_df)
  n_charged_partical_to_pair <- nrow(same_entry_df) - 1
  
  # only one partical
  if(n_charged_partical_to_pair < 1){
    return (as.data.frame(mass_list))
  }
  
  # iterate two particals and calculate mass by using E^2 = P^2 + m^2
  for(i in 1 : n_charged_partical_to_pair){
    first_partical <- same_entry_df[i,]
    first_partical$momentum_square <- (first_partical$px)^2 + (first_partical$py)^2 + (first_partical$pz)^2
    
    for(j in (i+1) : n_charged_partical){
      second_partical <- same_entry_df[j,]
      second_partical$momentum_square <- (second_partical$px)^2 + (second_partical$py)^2 + (second_partical$pz)^2
      
      first_add_second <- data.frame(momentum_square = (first_partical$px + second_partical$px)^2 + (first_partical$py + second_partical$py)^2 + (first_partical$pz + second_partical$pz)^2,
                                     pt_square = (first_partical$px + second_partical$px)^2 + (first_partical$py + second_partical$py)^2)
      
      if(first_partical$charge != second_partical$charge){
        #case1: pion_kaon
        first_partical$energy <- sqrt(first_partical$momentum_square + 0.13957018^2)
        second_partical$energy <- sqrt(second_partical$momentum_square + 0.493677^2)
        mass <- sqrt((first_partical$energy + second_partical$energy)^2 - first_add_second$momentum_square)
        pt <- sqrt(first_add_second$pt_square)
        
        if(pt>5){
          mass_list <- append(mass_list, mass)
        }
        
        #case2: kaon_pion
        first_partical$energy <- sqrt(first_partical$momentum_square + 0.493677^2)
        second_partical$energy <- sqrt(second_partical$momentum_square + 0.13957018^2)
        mass <- sqrt((first_partical$energy + second_partical$energy)^2 - first_add_second$momentum_square)
        pt <- sqrt(first_add_second$pt_square)
        
        if(pt>5){
          mass_list <- append(mass_list, mass)
        }
      }
    }
  }
  
  return (as.data.frame(mass_list))
}


calculate_meson_mass_main <- function(collision_data){
  cat('CALCULATE START: ', Sys.time(), '\n')
  collision_data <- collision_data %>% 
    group_by(entry) %>% 
    do(meson_mass_calculation(.))
  
  cat('CALCULATE FINISH: ', Sys.time(), '\n')
  return (collision_data)
}


calculate_meson_mass_main_parallel <- function(collision_data){
  distinct_entry <- unlist(collision_data %>% distinct(entry))
  cat('distinct_entry:', max(distinct_entry), '\n')
  
  n_cores <- detectCores() - 2  
  cl = makeCluster(n_cores, outfile = parallel_log)
  registerDoParallel(cl)
  
  cat('PARALLEL CALCULATE START: ', strftime(Sys.time(),"%Y-%m-%d %H:%M:%S"), '\n')
  meson_mass <- c()
  tryCatch({
    meson_mass <- foreach(i = distinct_entry, .combine = rbind,
                          .export=c('meson_mass_calculation', 'meson_mass_calculation_dict', 'px', 'py', 'pz'), 
                          .packages=c()) %dopar% {
                            # s <- sample(1:10, 1)
                            # if(s > 6) {
                            #   cat("entry: ", i, '\n')
                            # }
                            cat('entry:', i, '\n')
                            meson_mass_calculation_dict(collision_data[collision_data$entry==i,])
                            #meson_mass_calculation(collision_data[collision_data$entry==i,])
                          }
    
    cat('PARALLEL CALCULATE FINISH: ', strftime(Sys.time(),"%Y-%m-%d %H:%M:%S"), '\n')
    #return (as.data.frame(meson_mass))
    return (data.frame(mass_list = meson_mass))
  }, error = function(e) {
    print(e)
  },
  finally = {
    stopCluster(cl)
    rm('meson_mass', 'cl')
  })
  
}



##### PLOT #####
# plot_type: 'h' for histogram, 's' for scatter
# plot_opts: bins, x_tick,  x_range, y_lim, color, title (, y_range, x_lim)
plot_mass_figs <- function(meson_mass_df, plot_type='h', plot_opts){
  
  if(length(plot_opts$x_range)==0){
    plot_x_range <- c(min(meson_mass_df$mass_list), max(meson_mass_df$mass_list))
  }else{
    plot_x_range <- plot_opts$x_range
    meson_mass_df <- meson_mass_df[meson_mass_df$mass_list > plot_x_range[[1]][1] & 
                                     meson_mass_df$mass_list < plot_x_range[[1]][2],]
  }
 
  plot_title <- ifelse(length(plot_opts$title)==0, 'Plot of Mass Distribution', plot_opts$title)
  plot_bins <- ifelse(length(plot_opts$bins)==0, 100, plot_opts$bins)
  plot_color <- ifelse(length(plot_opts$color)==0, 'red', plot_opts$color)
  
  if(plot_type == 's'){
    #y <- meson_mass_df[meson_mass_df$mass_list>1.8 & meson_mass_df$mass_list<1.9,] %>% 
    #  mutate(bar_cat=trunc(mass_list*10^3)/10^3) %>% 
    #  group_by(bar_cat) %>% 
    #  dplyr::summarise(mean=mean(mass_list), std=sd(mass_list), cnt=n(), max=max(mass_list), min=min(mass_list))
    
    h <- hist(meson_mass_df, breaks = plot_bins, plot=FALSE)
    sort_count <- sort(h$counts, T)
    top_1 <- h$mids[which(h$counts==sort_count[1])]
    top_2 <- h$mids[which(h$counts==sort_count[2])]
    top_3 <- h$mids[which(h$counts==sort_count[3])]
    top_4 <- h$mids[which(h$counts==sort_count[4])]
    top_5 <- h$mids[which(h$counts==sort_count[5])]
    
    top5_df <- data.frame(x=c(top_1, top_2, top_3, top_4, top_5), y=sort_count[1:5])

    h_df <- data.frame(mids=h$mids, counts=h$counts)
  
    p <- ggplot(data =h_df, aes(x = mids, y = counts)) + 
      geom_point(shape = 3,  color = plot_color) + 
      ggtitle(plot_title) +
      labs(x = "mass", y = "count") +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_classic() 
    
    p <- p + 
      geom_point(data =top5_df, aes(x = x, y = y), shape = 1, color = 'blue', size = 5) + 
      geom_text(data = top5_df, aes(x = x, y = y, label = x, vjust = -1.5, hjust = 0.5))
  }else{
    meson_mass_df <- data.frame(mass_list = meson_mass_df)
    #p <- ggplot(data = meson_mass_df, aes(meson_mass_df$mass_list)) + 
    #  geom_histogram(bins = plot_bins) +
    #  theme_classic()
    
    p <- ggplot(data = meson_mass_df, aes(meson_mass_df$mass_list)) + 
      #stat_bin(geom = "step", bins = plot_bins, color = plot_color) +
      stat_bin(geom = "bar", bins = plot_bins, color = plot_color, alpha=0.3) +
      ggtitle(plot_title) + 
      labs(x = "mass", y = "count") +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme_classic()
  }

  ifelse(is.na(plot_opts$x_tick), 
         p <- p, 
         p <- p + scale_x_continuous(breaks = seq(plot_x_range[[1]][1], plot_x_range[[1]][2], plot_opts$x_tick)))
  
  ifelse(is.na(plot_opts$ylim),
         p <- p, 
         p <- p + coord_cartesian(ylim = unlist(plot_opts$ylim)))
  
  return (p)
}


plot_column_distribution_figs <- function(data_df, plot_bins=50, plot_grid_ncols=2){
  p_t <- ggplot(data_df, aes(x=pt)) + 
    geom_histogram(bins=plot_bins, colour="black", fill="#FF6666", alpha=.2) +
    #geom_histogram(aes(y=..density..), colour="black", fill="white")+
    #geom_density(alpha=.2, fill="#FF6666") +
    theme_classic()
  p_x <- ggplot(data_df, aes(x=px)) + 
    geom_histogram(bins=plot_bins, colour="black", fill="#E69F00", alpha=.2) +
    #geom_histogram(aes(y=..density..), colour="black", fill="white")+
    #geom_density(alpha=.2, fill="#FF6666") +
    theme_classic()
  p_y <- ggplot(data_df, aes(x=py)) + 
    geom_histogram(bins=plot_bins, colour="black", fill="#56B4E9", alpha=.2) +
    #geom_histogram(aes(y=..density..), colour="black", fill="white")+
    #geom_density(alpha=.2, fill="#FF6666") +
    theme_classic()
  p_z <- ggplot(data_df, aes(x=pz)) + 
    geom_histogram(bins=plot_bins, colour="black", fill="lightblue", alpha=.2) +
    #geom_histogram(aes(y=..density..), colour="black", fill="white")+
    #geom_density(alpha=.2, fill="#FF6666") +
    theme_classic()
  
  grid.arrange(p_t, p_x, p_y, p_z, ncol=plot_grid_ncols)
}



##### MAIN #####
#1: working_dir
#2: raw_csv_filename
#3: compressed_csv_filename
#4: exist_header=F
#5: nrows=NA
args = commandArgs(trailingOnly=TRUE)

#[TODO] case 1: original raw data
args = c('', 'xxx.csv', '', TRUE, NA)

#[TODO] case 2: directory contains compressed data and dictionary
args = c('', '', 'Rows', TRUE, NA)


#import data
collision_data <- import_data_source(args[1], args[2], args[3], as.logical(args[4]), as.numeric(args[5]))

#calculate mass
mass_result <- calculate_meson_mass_main_parallel(collision_data)

#plot mass distribution
plot_type <- 's'
plot_opts <- data.frame(x_range = I(list(c(1.8, 1.9))), bins=100, x_tick=0.02)
plot_mass_figs(mass_result, plot_type, plot_opts)

#plot column distribution
#distribution_cols <- c( 'px', 'py', 'pz')
#plot_column_distribution_figs(collision_data, distribution_cols, plot_bins=50, plot_grid_ncols=2)
plot_column_distribution_figs(collision_data, plot_bins=50, plot_grid_ncols=2)
