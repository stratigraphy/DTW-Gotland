#written by Michiel Arts michiel.arts@stratigraphy.eu
#Ensure that the right working directories are set to allow one to load the right files!!!
#load the R packages ####
library(ggplot2)
library(astrochron)
library(dtw)
library(plotly)
library(rlist)
library(data.table)
library(WaverideR)
library(dtwclust)
library(ggrepel)
library(dplyr)
library(matrixStats)
library(sf)
library(colorspace)
library(ggspatial)
library(rnaturalearth)
library(DescTools)
library(matrixStats)
library(jsonlite)

# Set the parent directory containing all sub folders
# change out 

# Set the parent directory
# parent_dir <- "https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/tree/main/SI_4_logs_used"
# subdirs <- list.dirs(parent_dir, recursive = FALSE, full.names = TRUE)

# Initialize list
all_txt_data <- list()
# 
# for (i in seq_along(subdirs)) {
#   folder_path <- subdirs[i]
#   folder_name <- basename(folder_path)
#   
#   txt_files <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE, ignore.case = TRUE)
#   if (length(txt_files) == 0) next
#   
#   for (file_path in txt_files) {
#     # Read first few lines to detect if header exists
#     first_line <- readLines(file_path, n = 1)
#     # Simple heuristic: check if first line contains letters (not only numbers)
#     has_header <- grepl("[A-Za-z]", first_line)
#     
#     if (has_header) {
#       # Read with header
#       df <- read.table(file_path, header = TRUE, sep = "", stringsAsFactors = FALSE)
#     } else {
#       # Read without header, then assign colnames manually
#       df <- read.table(file_path, header = FALSE, sep = "", stringsAsFactors = FALSE)
#       # Assign colnames (assuming exactly two columns)
#       colnames(df)[1:2] <- c("DEPTH", "GR")
#     }
#     
#     # Assign a name for the list element
#     file_name <- tools::file_path_sans_ext(basename(file_path))
#     list_name <- paste0(file_name)
#     
#     all_txt_data[[list_name]] <- df
#   }
#}

####

library(gh)
library(gert)

# Clone the repo to a temp directory

Sys.setenv(GITHUB_PAT = "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")
pat <- Sys.getenv("GITHUB_PAT")




repo_path <- tempfile("DTW-Gotland_")
git_clone("git@github.com:stratigraphy/DTW-Gotland.git", repo_path)

all_txt_data <- list()


# 3. Modify the recursive function to include the Authorization header
list_github_txt_files <- function(api_url, token) {
  # Use an Authorization header for authentication
  response <- httr::GET(
    url = api_url,
    httr::add_headers(
      Authorization = paste("token", token),
      Accept = "application/vnd.github.v3+json"
    )
  )
  
  # Check for HTTP errors (like 403 Forbidden)
  if (httr::http_error(response)) {
    stop(paste("GitHub API Error:", httr::status_code(response), httr::content(response, "text", encoding = "UTF-8")))
  }
  
  # Parse the JSON content
  data <- jsonlite::fromJSON(httr::content(response, "text", encoding = "UTF-8"), flatten = TRUE)
  
  txt_files <- character(0)
  
  for (i in seq_along(data$type)) {
    if (data$type[i] == "file" && grepl("\\.txt$", data$name[i], ignore.case = TRUE)) {
      # The 'download_url' is directly provided by the API response
      txt_files <- c(txt_files, data$download_url[i])
    } else if (data$type[i] == "dir") {
      # The API response already contains the 'url' for the subdirectory
      subdir_url <- data$url[i]
      txt_files <- c(txt_files, list_github_txt_files(subdir_url, token))
    }
  }
  return(txt_files)
}


# Top-level directory API URL
root_api <- "https://api.github.com/repos/stratigraphy/DTW-Gotland/contents/SI_4_logs_used"

# You'll need to load the required libraries
if (!requireNamespace("httr", quietly = TRUE)) install.packages("httr")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
library(httr)
library(jsonlite)

# Call the function with the root URL and the PAT
txt_urls <- list_github_txt_files(root_api, pat)

txt_urls[1:5]  # preview first few

# Read each text file into a named list
well_logs_list <- list()

for (file_url in txt_urls) {
  first_line <- tryCatch(readLines(file_url, n = 1), error = function(e) return(NA))
  if (is.na(first_line)) next
  
  has_header <- grepl("[A-Za-z]", first_line)
  
  df <- tryCatch({
    if (has_header) {
      read.table(file_url, header = TRUE, sep = "", stringsAsFactors = FALSE)
    } else {
      df <- read.table(file_url, header = FALSE, sep = "", stringsAsFactors = FALSE)
      colnames(df)[1:2] <- c("DEPTH", "GR")
      df
    }
  }, error = function(e) NULL)
  
  if (!is.null(df)) {
    file_name <- tools::file_path_sans_ext(basename(file_url))
    well_logs_list[[file_name]] <- df
    list_name <- paste0(file_name)
    all_txt_data[[list_name]] <- df
  }
}

# Check results
str(all_txt_data)
sort(names(all_txt_data))
length(all_txt_data)
# rename to well-log list
well_logs_list <- all_txt_data

#remove the _GR. extension from the file names in the list of well-logs
well_logs_list <- all_txt_data
for (col in 1:length(well_logs_list)) {
  names(well_logs_list)[[col]] <-
    gsub("_GR.*", "", names(well_logs_list)[[col]])
}

well_logs_list <- well_logs_list[order(names(well_logs_list))]

#load the d13C record of the Altajme core ####

Altamje_d13C <- read.csv("https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/main/Altamje_d13C.csv")

Altamje_d13C <- demean(Altamje_d13C, genplot = F, verbose = T)
Altamje_d13C[, 2]  <- Altamje_d13C[, 2] - min(Altamje_d13C[, 2])
Altamje_d13C[, 2]  <- Altamje_d13C[, 2] / max(Altamje_d13C[, 2])
Altamje_d13C <- iso(Altamje_d13C,
                    xmax = 324,
                    xmin = 0,
                    genplot = FALSE)



names(well_logs_list)
graphics.off()
plot(well_logs_list$STSUTARVE_2018[,c(2,1)],type="l",ylim=c(500,0),xlim=c(0,250))
abline(h=c(451.23,	411.17,	313.35,	292,	239.34,	115.36),col="red")



# Load required package
library(plotly)

# Extract the data
data <- well_logs_list$SKALS_1
x <- data[, 2]
y <- data[, 1]

# Create a plotly interactive plot
fig <- plot_ly(
  x = x,
  y = y,
  type = 'scatter',
  mode = 'lines',
  hoverinfo = 'x+y',
  line = list(color = 'blue')
)

# Reverse y-axis if needed (to match ylim = c(500, 0))
fig <- fig %>% layout(
  yaxis = list(autorange = "reversed", title = "Depth"),
  xaxis = list(title = "Value"),
  title = "Interactive Well Log Plot"
)

fig

# load the well-tops/tie-points
well_tops_rotated <- read.csv("https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/main/well_tops_rotated.csv")

exclude_wells <- c("Skaggs_1", "Sanda_1", "Rings_1", "Grotlingbo_2")
wells_incl <- well_tops_rotated[!well_tops_rotated$Well %in% exclude_wells, ]
wells_incl

readjusted_matrix <- wells_incl

readjusted_matrix_means <-
  matrix(
    data = NA,
    nrow = nrow(readjusted_matrix),
    ncol = (ncol(readjusted_matrix) - 4)
  )

# fill the matrix with thickness of the beds
for (i in 4:(ncol(readjusted_matrix) - 1)) {
  a <- (i - 3)
  top <- as.matrix(unlist(readjusted_matrix[i + 1]))
  top[is.na(top)] <- 0
  bottom <- as.matrix(unlist(readjusted_matrix[i]))
  readjusted_matrix_means[, a] <- bottom - top
  
  
}
readjusted_matrix_means

readjusted_matrix_means <- as.data.frame(readjusted_matrix_means)
readjusted_matrix_means$top <-
  as.matrix(unlist(readjusted_matrix[, (ncol(readjusted_matrix))]))
colnames(readjusted_matrix_means) <-
  colnames(readjusted_matrix[, 4:ncol(readjusted_matrix)])


#calculate the strechfactor fot the beds
strech_fact <-
  1 / ((readjusted_matrix_means / matrix(
    rep(
      colMaxs(as.matrix(readjusted_matrix_means[c(2:37),]), na.rm = TRUE),
      each = nrow(readjusted_matrix_means)
    ),
    nrow = nrow(readjusted_matrix_means),
    ncol = ncol(readjusted_matrix_means)
  )))

strech_fact <- as.data.frame(strech_fact)
strech_fact <- cbind(readjusted_matrix[, 3], strech_fact)
colnames(strech_fact) <-
  names(readjusted_matrix[, 3:ncol(readjusted_matrix)])

well_padded <- readjusted_matrix
well_padded <- as.data.frame(as.matrix(well_padded))
well_padded <- well_padded[, 4:ncol(well_padded)]
well_padded <-
  lapply(well_padded, function(well_padded)
    as.numeric(as.character(well_padded)))
well_padded <- as.data.frame(well_padded)
colmaxs <- colMeans(as.matrix(well_padded), na.rm = TRUE)


well_padded$top_surface <-
  well_padded[, ncol(well_padded)] - colmaxs[length(colmaxs)]
base_data <- well_padded[, 1] + 20
well_padded <- cbind(base_data, well_padded)
well_padded <- cbind(strech_fact[, 1], well_padded)
colnames(well_padded)[1] <- c("well")

maxvals <-
  colMeans(as.matrix(readjusted_matrix_means), na.rm = TRUE)
maxvals <- c(20, maxvals)
maxvals <- round((maxvals / 0.1), 0)
maxvals_depth <- maxvals

for (i in ((length(maxvals) - 1):1)) {
  maxvals_depth[i] <- maxvals_depth[i] + maxvals_depth[i + 1]
}
maxvals_depth[length(maxvals_depth) + 1] <- 0



# stretch the data to a common depth scale using the mean thickness
prestreched_data <- list()

for (hln in 1:nrow(well_padded)) {
  norm_log <- matrix(ncol = 3)
  norm_log <- as.data.frame(norm_log[-1, ])
  colnames(norm_log) <- c("norm_depth", "proxy", "depth_old")
  index <- grep(well_padded[hln, 1], names(well_logs_list), ignore.case = TRUE)
  well_1 <- well_logs_list[[index[1]]]
  well_1 <- sortNave(well_1, genplot = F, verbose = F)
  well_1 <- linterp(well_1, 0.1, genplot = F, verbose = F)
  well_1 <- demean(well_1, genplot = F, verbose = T)
  well_1[, 2]  <- well_1[, 2] - min(well_1[, 2])
  well_1[, 2]  <- well_1[, 2] / max(well_1[, 2])
  well_1 <- demean(well_1, genplot = F, verbose = T)
  
  #head(well_1)
  #plot(well_1,type="l")
  
  for (i in 1:(length(maxvals_depth) - 1)) {
    well_1_sel <- well_1[well_1[, 1] < well_padded[hln, i + 1], ]
    
    well_1_sel <-
      well_1_sel[well_1_sel[, 1] > well_padded[hln, i + 2], ]
    new_depths <-
      seq(from = maxvals_depth[i + 1] + 1,
          to = maxvals_depth[i],
          by = 1)
    
    if (nrow(well_1_sel) == 0  | any(is.na(well_1_sel))) {
      pad_NA <- rep(NA, times = length(new_depths))
      seq_1 <- cbind(new_depths, pad_NA, pad_NA)
      
    } else{
      app <-
        approx(
          x = well_1_sel[, 1],
          y = well_1_sel[, 2],
          method = "linear",
          n = length(new_depths),
          rule = 1,
          f = 0,
          ties = mean,
          na.rm = TRUE
        )
      seq_1 <- cbind(new_depths, app$y, app$x)
    }
    
    colnames(seq_1) <- c("norm_depth", "proxy", "depth_old")
    norm_log <- rbind(norm_log, seq_1)
    
  }
  
  norm_log <- norm_log[order(norm_log[, 1], decreasing = FALSE), ]
  prestreched_data <- list.append(prestreched_data, norm_log)
  
}

#if you want to save or load the data
#saveRDS(prestreched_data, file = "prestreched_data.rds")
#prestreched_data <- readRDS("prestreched_data.rds")

#order the data according to the base of the Silurian
prestreched_data <-
  prestreched_data[order(well_padded[, 3], decreasing = T)]
well_padded_ord <- well_padded[order(well_padded[, 3], decreasing = T),]

prestreched_data <- prestreched_data[c(1:34,37,38,35,36)]
well_padded_ord <- well_padded_ord[c(1:34,37,38,35,36),]

names(prestreched_data) <- well_padded_ord[,1]

#create a new matrix to fill with the pre-stretched data
data <- matrix(data = NA, nrow = nrow(prestreched_data[[1]]))
#prestreched_data
#nrow(prestreched_data[[1]])


#fill the matrix with the pre-stretched data which was also normalized (-1 to 1)
for (i in 1:length(prestreched_data)) {
  data_1 <- prestreched_data[[i]]
  data_1 <- data_1[, c(1, 2)]
  data_1[, 2] <- data_1[, 2] - mean(na.omit(data_1[, 2]))
  data_1 <- data_1[, 2]
  data_1 <- data_1 / (max(na.omit(data_1)) - min(na.omit(data_1)))
  data_1 <- data_1 / sd(na.omit(data_1))
  data_1 <- data_1 - mean(na.omit(data_1))
  data_1 <- data_1 / (max(na.omit(data_1)) - min(na.omit(data_1)))
  data_1 <- data_1 - mean(na.omit(data_1))
  data <- cbind(data, data_1)
}
#remove the first empty column
data <- data[, -c(1)]

#create a median curve###
cent <- rowMedians(as.matrix(data), na.rm = TRUE)

#convert the stretched data to a list
data_ls <- as.list(as.data.frame(data))

#pad NA's with the  values of the centroid
for (i in 1:length(data_ls)) {
  data_ls[[i]][is.na(data_ls[[i]])] <- cent[is.na(data_ls[[i]])]
}


exclude_wells <- c("Skaggs_1",
                   "Sanda_1",
                   "Rings_1",
                   "Grotlingbo_2",
                   "Grotlingbo_1",
                   "Grotlingbo_2_Res")

well_logs_list_2 <- well_logs_list[!names(well_logs_list) %in% exclude_wells]
well_tops_rotated_2 <- well_tops_rotated[!well_tops_rotated$Well %in% exclude_wells,]

names_prestretched <- tolower(names(prestreched_data))
names_well_logs <- tolower(names(well_logs_list_2))

# Match names (lowercase) and get index
matched_idx <- match(names_prestretched, names_well_logs)

# Keep only matched and existing entries
valid_idx <- which(!is.na(matched_idx))
matched_order <- matched_idx[valid_idx]

# Reorder well_logs_list_2 to match prestretched_data
well_logs_list_3 <- well_logs_list_2[matched_order]
#well_tops_rotated_2 <- well_tops_rotated_2[matched_order,]


target_order <- names(well_logs_list_3)
well_tops_wells <- well_tops_rotated_2$Well
match_index <- match(toupper(well_tops_wells), toupper(target_order))
sorting_sequence <- order(match_index)
well_tops_rotated_2 <- well_tops_rotated_2[sorting_sequence, ]




# Optional: assign matching names for clarity
names(well_logs_list_3) <- names(prestreched_data)[valid_idx]

for (i in 1:length(well_logs_list_3)) {
  well_1 <- well_logs_list_3[[i]]
  well_1 <- well_1[well_1[, 1] <= max(prestreched_data[[i]][, 3], na.rm = TRUE), ]
  well_1 <- sortNave(well_1, genplot = FALSE)
  well_1 <- linterp(well_1, genplot = FALSE)
  well_1 <- demean(well_1, genplot = F, verbose = T)
  well_1[, 2]  <- well_1[, 2] - min(well_1[, 2])
  well_1[, 2]  <- well_1[, 2] / max(well_1[, 2])
  well_1 <- demean(well_1, genplot = F, verbose = T)
  well_logs_list_3[[i]] <- well_1
}


graphics.off()




###plot figure 2####
{
  
  
  layout.matrix <-
    matrix(c(1:2), nrow = 2, ncol = 1)
  graphics::layout(mat = layout.matrix,
                   heights = c(1),
                   # Heights of the two rows
                   widths = c(1)) # Widths of the two columns
  par(mar = c(2, 4, 2, 2))
  
  plot(
    well_logs_list_3[[1]][, c(2)],
    y = well_logs_list_3[[1]][, c(1)],
    type = "l",
    xlim = c(0, length(well_logs_list_3) * 0.75),
    ylim = c(520, 0),xaxt="n",ylab="depth (m)",xlab=""
  )
  
  for (i in 2:length(well_logs_list_3)) {
    lines(x = well_logs_list_3[[i]][, c(2)] + (i * (0.75) - 1),
          y = well_logs_list_3[[i]][, c(1)],
          col = "black")
  }
  
  #add the d13C record of the Altajme core to the plot
  lines(x = c(Altamje_d13C[, 2] + ((
    length(prestreched_data) * 0.75
  ) - 0.5)),
  y = Altamje_d13C[, 1],
  col = "green")
  
  y_vals <- well_tops_rotated_2[,4:ncol(well_tops_rotated_2)]
  x_vals <- y_vals
  y_vals
  
  for (i in 1:ncol(y_vals)) {
    for (j in 1:length(well_logs_list_3)) {
      
      x_val_well <- (well_logs_list_3[[j]][, c(2)] + (j * (0.75) - 1))
      y_val_well <- (well_logs_list_3[[j]][, c(1)])
      val_well <- cbind(x_val_well, y_val_well)
      val_well <- val_well[complete.cases(val_well), ]
      row_nr <- Closest(val_well[, 2], y_vals[j,i], which = TRUE)
      if (!is.na(row_nr[1])){x_vals[j,i] <- val_well[row_nr[1]]}
    }
    lines(x = x_vals[, i], y = y_vals[, i], col = "red")
    
    
  }
  
  
  
  
  
  
  
  # plot the prestreched data
  plot(
    prestreched_data[[1]][, c(2)],
    prestreched_data[[1]][, c(1)]/10,
    type = "l",
    xlim = c(0, length(data_ls) * 0.75),
    ylim = c(4500/10, 0),xaxt="n",xlab="",ylab="normalized (depth scale)"
  )
  for (i in 2:length(data_ls)) {
    lines(x = prestreched_data[[i]][, c(2)] + (i * (0.75) - 1),
          y = prestreched_data[[i]][, c(1)]/10,
          col = "black")
  }
  # add the formation/tie-points boundaries
  abline(h = maxvals_depth/10, col = "red")
}




# run the DTW with Barycenter averaging ####

dtw_avg <- DBA(
  data_ls,
  centroid = cent,
  max.iter = 2000L,
  delta = 0.0025,
  trace = TRUE,
  norm = "L2",
  normalize = TRUE,
  step.pattern = dtw::symmetric2,
  error.check = TRUE,
  window.size = 200,
  backtrack = TRUE,
  sqrt.dist = TRUE,
  mv.ver = "by-variable"
)

#Denoise and post-process the DBA results and select the 150 tie-points ####
a <- cbind(prestreched_data[[1]][1], dtw_avg)

a <- taner(
  a,
  flow = 1 / 30,
  fhigh = 1 / (max(a[1]) * 2),
  roll = 10 ^ 20,
  detrend = FALSE,
  demean = FALSE
)


a_nl <- noLow(a, smooth = 0.075)
a_mx <- max_detect(a_nl, pts = 5)
a_mn <- min_detect(a_nl, pts = 5)


a_pks <- rbind(a_mx, a_mn)
a_pks[, 2] <- abs(a_pks[, 2])

a_pks <- a_pks[order(-a_pks[, 2]), ]

n_tiepoints <- 150

tie_points <- matrix(data = NA, ncol = 3, nrow = n_tiepoints)
tie_points <- as.data.frame(tie_points)
min_offset <- 0.25

for (i in 1:n_tiepoints) {
  tie_points[i, ] <- a_pks[1, ]
  a_pks <- a_pks[-c(1), ]
  vals <- abs(a_pks[, 1] - tie_points[i, 1])
  a_pks <- a_pks[vals > min_offset, ]
}


layout.matrix <-
  matrix(c(1:3), nrow = 3, ncol = 1)
graphics::layout(mat = layout.matrix,
                 heights = c(1),
                 # Heights of the two rows
                 widths = c(1)) # Widths of the two columns
par(mar = c(4, 4, 0.5, 0))

plot(
  cbind(prestreched_data[[1]][1], dtw_avg),
  type = "l",
  xlab = "normalized depth scale",
  ylab = "normalized gamma=ray"
)
lines(a, col = "red", lwd = 2)

plot(a_nl,
     type = "l",
     xlab = "normalized depth scale",
     ylab = "normalized gamma=ray")

points(a_nl[tie_points[, 1], ],
       col = "blue",
       pch = 19,
       cex = 2)

plot(
  cbind(prestreched_data[[1]][1], dtw_avg),
  type = "l",
  lwd = 2,
  xlab = "normalized depth scale",
  ylab = "normalized gamma=ray"
)
lines(a, col = "red", lwd = 2)
points(a[tie_points[, 1], ],
       col = "blue",
       pch = 19,
       cex = 2)



# re-correlate the wells to the referencecurve ####

a_pks <- tie_points
tie_points <- matrix(data = NA, ncol = 1, nrow = n_tiepoints)

for (i  in 1:nrow(tie_points)) {
  row_nr_1 <-
    DescTools::Closest(as.matrix(prestreched_data[[1]][1]), a_pks[i, 1], which = TRUE)
  tie_points[i, 1] <- row_nr_1[1]
}
nrow(tie_points)

corr_points <- list()
data_tuned_cor <- list()


for (i in 1:length(data_ls)) {
  
  x <- taner(
    cbind(seq(1:length(data_ls[[i]])), data_ls[[i]]),
    flow = 1 / 30,
    fhigh = 1 / (max(a[1]) * 2),
    roll = 10 ^ 20,
    detrend = FALSE,
    genplot = FALSE,
    verbose = FALSE
  )
  
  dtw_basic_1 <- dtw_basic(
    x = x[, 2],
    y = a[, 2],
    norm = "L2",
    normalize = TRUE,
    step.pattern = dtw::symmetric2,
    error.check = TRUE,
    window.size = 200,
    backtrack = TRUE,
    sqrt.dist = TRUE,
    mv.ver = "by-variable"#,
    #open.end = TRUE,
    #open.begin = TRUE
  )
  
  
  tie_points_2 <- tie_points
  for (j  in 1:nrow(tie_points)) {
    row_nr_1 <-
      DescTools::Closest(dtw_basic_1[[2]], tie_points[j, 1], which = TRUE)
    tie_points_2[j, 1] <- row_nr_1[1]
  }
  
  vals <- prestreched_data[[i]][2][dtw_basic_1[[2]][tie_points_2], ]
  depths <-
    prestreched_data[[i]][3][dtw_basic_1[[2]][tie_points_2], ]
  corr_pts <- cbind(depths, vals)
  
  time_cor <-
    cbind(prestreched_data[[1]][1][dtw_basic_1[[2]][tie_points_2], ], prestreched_data[[1]][1][dtw_basic_1[[3]][tie_points_2], ])
  
  time_cor <- sortNave(time_cor[,c(1,2)],genplot = FALSE)
  time_cor <- sortNave(time_cor[,c(2,1)],genplot = FALSE)
  time_cor <- sortNave(time_cor[,c(2,1)],genplot = FALSE)
  
  cbind(prestreched_data[[1]][1], prestreched_data[[i]][2])
  
  data_ls_tuned <-
    tune(cbind(prestreched_data[[1]][1], prestreched_data[[i]][2]),
         time_cor,
         genplot = FALSE)
  
  data_ls_tuned[, 1] <- data_ls_tuned[, 1] / 10
  graphics.off()
  plot(data_ls_tuned,type="l")
  
  
  plot(cbind(prestreched_data[[1]][1], prestreched_data[[i]][2]),type="l")
  
  data_tuned_cor <- list.append(data_tuned_cor, data_ls_tuned)
  
  corr_points <- list.append(corr_points, corr_pts)
}

# plot the tie-points for all the wells ####
graphics.off()
colors <- rainbow(n_tiepoints)

layout.matrix <- matrix(1,
                        nrow = 1,
                        ncol = 1,
                        byrow = TRUE)
layout(layout.matrix, heights = c(1), widths = c(1))

par(mar = c(4, 4, 1, 0.5))


plot(
  x = c(-2.2, length(data_tuned_cor) * 1.2),
  y = c(525, 0),
  col = "white",
  xlab = "",
  ylab = "depth (m)",
  xaxt = "n",
  xaxs = "i",
  yaxs = "i",
  ylim = c(520, 0)
)            # Draw empty plot

points(
  y = (a[tie_points, 1] / 10) + prestreched_data[[1]][1, 3],
  x = a[tie_points, 2] - 1 ,
  col = colors,
  pch = 19,
  cex = 0.75
)


lines(a[, 2] - 1, (a[, 1] / 10) + prestreched_data[[1]][1, 3])


for (i in 1:length(prestreched_data)) {
  well_line <- cbind(well_logs_list_3[[i]][1], well_logs_list_3[[i]][2])
  pts <- corr_points[[i]]
  for (j in 1:nrow(pts)) {
    rownr <- DescTools::Closest(well_line[, 1], pts[j, 1], which = TRUE)
    if (!is.na(rownr[1])) {
      pts[j, 2] <- well_line[rownr[1], 2]
    }
  }
  
  corr_points[[i]] <- pts
  
  points(
    y = pts[, 1],
    x = pts[, 2] + (i * 1.2) - 1,
    col = colors,
    pch = 19,
    cex = 0.75
  )
  
}



depth_mat <- matrix(data = NA,
                    nrow = n_tiepoints,
                    ncol = length(data_tuned_cor))
val_mat <- matrix(data = NA,
                  nrow = n_tiepoints,
                  ncol = length(data_tuned_cor))

for (i in 1:length(prestreched_data)) {
  pts <- corr_points[[i]]
  depth_mat[, i] <- pts[, 1]
  val_mat[, i] <- pts[, 2] + (i * 1.2) - 1
}


depth_mat <- cbind((a[tie_points, 1] / 10) + prestreched_data[[1]][1, 3], depth_mat)
val_mat <-  cbind(a[tie_points, 2] - 1 , val_mat)



for (i in 1:nrow(val_mat)) {
  vals <- as.data.frame(cbind(val_mat[i, ], depth_mat[i, ]))
  lines(vals, col = colors[i])
  
}

for (i in 1:length(prestreched_data)) {
  lines(cbind(well_logs_list_3[[i]][2] + (i * 1.2) - 1, well_logs_list_3[[i]][1]))
}


#pick the formational  boundaries ####
{
names(well_logs_list_3)  
  
  
  
  
# base hogklint
# tracked manually !!!!

#x_vals_base_hog <-  ((seq(1:length(well_logs_list_3)) - 1) * 1.2) - 1

names(well_logs_list_3)

graphics.off()
iop <- 2
well_sel <- well_logs_list_3[[iop]]
plot(well_sel[,2],well_sel[,1],type="l",ylim=c(200,0))

pts <- corr_points[[iop]]

points(
  pts[,2],pts[,1],
  col = colors,
  pch = 19,
  cex = 0.75
)
abline(h=120)






base_hogklint <- c(
  438, # ref
  418,#1
  418,#2
  401,#3
  402,#4
  390,#5
  397,#6
  400.5,#7
  386,#8
  372.5,#9
  381,#10
  372.5,#11
  367.5,#12
  366.75,#13
  352.5,#14
  350,#15
  350,#16
  342,#17
  336,#18
  332.5,#19
  332.5,#20
  337,#21
  327,#22
  328,#23
  330.5,#24
  321,#25
  309.5,#26
  316.25, #27
  304.5,#28
  321,#29
  315,#30
  316.5,#31
  304,#32
  311.5,#33
  313.5,#34
  258,#35
  251,#36
  293.75,#37
  271.6#38
)

vals_ord <- val_mat
depth_mat_ord <- depth_mat

depth_mat_ord <- depth_mat_ord[order(depth_mat[, 1]), ]
vals_ord <- vals_ord[order(depth_mat[, 1]), ]

vals_ord <- vals_ord[, 1:length(base_hogklint)]
depth_mat_ord <- depth_mat_ord[, 1:length(base_hogklint)]


plot(a[, 2] - 1, (a[, 1] / 10) + prestreched_data[[1]][1, 3],
     type="l",ylim=c(425,69))

points(
  y = (a[tie_points, 1] / 10) + prestreched_data[[1]][1, 3],
  x = a[tie_points, 2] - 1 ,
  col = colors,
  pch = 19,
  cex = 0.75
)
#abline(h=73)


str(depth_mat_ord)

plot(a[, 2] - 1, (a[, 1] / 10) + prestreched_data[[1]][1, 3],type="l",
     ylim = c(350, 0))

points(
  y = (a[tie_points, 1] / 10) + prestreched_data[[1]][1, 3],
  x = a[tie_points, 2] - 1 ,
  col = colors,
  pch = 19,
  cex = 0.75
)

abline(h=260)

#base Silurian
base_silurian <- DescTools::Closest(depth_mat_ord[, 1], 493, which = TRUE)

base_L_Visby <- DescTools::Closest(depth_mat_ord[, 1], 482.6  , which = TRUE)

#base upper Visby
base_U_visby <- DescTools::Closest(depth_mat_ord[, 1], 447, which = TRUE)

# base hogklint
# tracked manually see above

#base Tofta
base_Tofta <- DescTools::Closest(depth_mat_ord[, 1], 432, which = TRUE)
#base Hangvar
base_Hangvar <- DescTools::Closest(depth_mat_ord[, 1], 425, which = TRUE)

#base Slite
base_Slite <- DescTools::Closest(depth_mat_ord[, 1], 415, which = TRUE)
#base Frojel
base_Frojel <- DescTools::Closest(depth_mat_ord[, 1], 353  , which = TRUE)
#base Klinteberg
base_Klinteberg <- DescTools::Closest(depth_mat_ord[, 1], 339.7 , which = TRUE)

depth_mat_ord[,1]
depth_mat_ord[131,1]

#base Hemse
base_Hemse <- DescTools::Closest(depth_mat_ord[, 1], 314, which = TRUE)
#base Etelheim  #start Linde ?
base_Etelhem <- DescTools::Closest(depth_mat_ord[, 1], 260, which = TRUE)
#base Nar 
base_Nar <- DescTools::Closest(depth_mat_ord[, 1], 217, which = TRUE)
#base EKE (option 2)
base_EKE <- DescTools::Closest(depth_mat_ord[, 1], 135, which = TRUE)
#base Burgsvik
base_Burgsvik <- DescTools::Closest(depth_mat_ord[, 1], 118, which = TRUE)
#base Hamra Sundre
base_Hamra_Sundre <- DescTools::Closest(depth_mat_ord[, 1], 71, which = TRUE)

# graphics.off()
# plot(well_logs_list_3[[1]][,2],well_logs_list_3[[1]][,1],type="l",ylim=c(400,300))
# abline(h=307.5)
#depth_mat_ord[,1]

readjusted_matrix_ord <- readjusted_matrix
Category <- factor(as.character(readjusted_matrix_ord$Well), levels = c(names(prestreched_data)))
readjusted_matrix_ord <- readjusted_matrix_ord[order(Category), ]

readjusted_matrix_ord <-rbind(readjusted_matrix_ord[1,],readjusted_matrix_ord)
readjusted_matrix_ord[1,] <- NA
readjusted_matrix_ord[1,3] <- c("Ref_curve")


readjusted_matrix_ord$base_silurian <- depth_mat_ord[base_silurian, ]
readjusted_matrix_ord$base_L_visby <- depth_mat_ord[base_L_Visby, ]
readjusted_matrix_ord$base_U_visby <- depth_mat_ord[base_U_visby, ]
readjusted_matrix_ord$base_hogklint <- base_hogklint
readjusted_matrix_ord$base_Tofta <- depth_mat_ord[base_Tofta, ]
readjusted_matrix_ord$base_Hangvar <- depth_mat_ord[base_Hangvar, ]
readjusted_matrix_ord$base_Slite <- depth_mat_ord[base_Slite, ]
readjusted_matrix_ord$base_Frojel <- depth_mat_ord[base_Frojel, ]
readjusted_matrix_ord$base_Klinteberg <- depth_mat_ord[base_Klinteberg, ]
readjusted_matrix_ord$base_Hemse <- depth_mat_ord[base_Hemse, ]
readjusted_matrix_ord$base_Etelhem <- depth_mat_ord[base_Etelhem, ]
readjusted_matrix_ord$base_Nar <- depth_mat_ord[base_Nar, ]
readjusted_matrix_ord$base_EKE <- depth_mat_ord[base_EKE, ]
readjusted_matrix_ord$base_Burgsvik <- depth_mat_ord[base_Burgsvik, ]
readjusted_matrix_ord$base_Hamra_Sundre <- depth_mat_ord[base_Hamra_Sundre, ]

#View(readjusted_matrix_ord)

#str(readjusted_matrix_ord)

#names(well_logs_list_3)[sel_wells-1]

sel_wells <- c(2, 3, 4,12,17,30,33,39)
prestreched_data_2 <- well_logs_list_3[sel_wells-1]
corr_points_2 <- corr_points[sel_wells-1]
readjusted_matrix_ord_2 <- readjusted_matrix_ord[sel_wells, ]


readjusted_matrix_ord_2$distance <- NA
readjusted_matrix_ord_2$distance[1] <- 0

for (i in 2:nrow(readjusted_matrix_ord_2)) {
  readjusted_matrix_ord_2$distance[i] <- sqrt(((readjusted_matrix_ord_2$X[i -
                                                                            1] - readjusted_matrix_ord_2$X[i]) ^ 2
  ) +
    ((readjusted_matrix_ord_2$Y[i -
                                  1] - readjusted_matrix_ord_2$Y[i]) ^ 2
    ))
}

readjusted_matrix_ord_2$distance_sum <- cumsum(readjusted_matrix_ord_2$distance)


}




#plot all the correlations lines in the transect ####
{
layout.matrix <- matrix(1,
                        nrow = 1,
                        ncol = 1,
                        byrow = TRUE)
layout(layout.matrix, heights = c(1), # Heights of the two rows
       widths = c(1))

par(mar = c(4, 4, 1, 0.5))


plot(
  x = c(-1000, max(readjusted_matrix_ord_2$distance_sum) * 1.1),
  y = c(525, 0),
  col = "white",
  xlab = "m along transect",
  ylab = "depth (m)",
  yaxs = "i",
  ylim = c(520, 0)
)            # Draw empty plot



comp_fact <- 4000

for (i in 1:length(prestreched_data_2)) {
  pts <- corr_points_2[[i]]
  points(
    y = pts[, 1],
    x = pts[, 2] * comp_fact + readjusted_matrix_ord_2$distance_sum[i],
    col = colors,
    pch = 19
  )
  
}

depth_mat_2 <- matrix(data = NA,
                      nrow = n_tiepoints,
                      ncol = length(data_tuned_cor))
val_mat_2 <- matrix(data = NA,
                    nrow = n_tiepoints,
                    ncol = length(data_tuned_cor))

for (i in 1:length(prestreched_data_2)) {
  pts <- corr_points_2[[i]]
  depth_mat_2[, i] <- pts[, 1]
  val_mat_2[, i] <- pts[, 2] * comp_fact + readjusted_matrix_ord_2$distance_sum[i]
}


for (i in 1:nrow(val_mat_2)) {
  #i <-35
  vals <- as.data.frame(cbind(val_mat_2[i, ], depth_mat_2[i, ]))
  lines(vals, col = colors[i])
  
}

for (i in 1:length(prestreched_data_2)) {
  lines(
    cbind(
      prestreched_data_2[[i]][2] * comp_fact + readjusted_matrix_ord_2$distance_sum[i],
      prestreched_data_2[[i]][1]
    )
  )
}
}

names(prestreched_data_2)

#plot the formations as polygons #####
{
layout.matrix <- matrix(1,
                        nrow = 1,
                        ncol = 1,
                        byrow = TRUE)
layout(layout.matrix, heights = c(1), # Heights of the two rows
       widths = c(1))

par(mar = c(4, 4, 1, 0.5))


plot(
  x = c(-10000, max(readjusted_matrix_ord_2$distance_sum) * 1.1),
  y = c(525, 0),
  col = "white",
  xlab = "m along transect",
  ylab = "depth (m)",
  ylim = c(520, 0)
)            # Draw empty plot


comp_fact <- 4000

for (i in 1:length(prestreched_data_2)) {
  lines(
    cbind(
      prestreched_data_2[[i]][2] * comp_fact + readjusted_matrix_ord_2$distance_sum[i],
      prestreched_data_2[[i]][1]
    )
  )
  
  
}

depth_mat_sel <- matrix(data = NA, nrow = 8, ncol = 21)
val_mat_sel <- matrix(data = NA, nrow = 8, ncol = 21)

#readjusted_matrix_ord_2
names(readjusted_matrix_ord_2)

for (i in 1:8) {
  dt <-  as.data.frame(prestreched_data_2[[i]])
  dt <- na.omit(dt)
  x_vals <- as.data.frame(dt[, 2])
  y_vals <- as.data.frame(dt[, 1])
  for (j in 4:24) {
    if (!is.na(readjusted_matrix_ord_2[i, j])) {
      row_nr <- DescTools::Closest(readjusted_matrix_ord_2[i, j], y_vals, which = TRUE)
      val_mat_sel[i, j - 3] <- x_vals[row_nr[1], ]
      depth_mat_sel[i, j - 3] <- y_vals[row_nr[1], ]
      y_vals[row_nr[1], ]
      
    }
  }
}


depth_mat_sel <- depth_mat_sel[, c(7:ncol(depth_mat_sel))]
val_mat_sel <- val_mat_sel[, c(7:ncol(val_mat_sel))]
val_mat_sel <- val_mat_sel * comp_fact + readjusted_matrix_ord_2$distance_sum



Altajme_streched <- na.omit(
  cbind(
    well_logs_list_3[[38]][2] * comp_fact + readjusted_matrix_ord_2$distance_sum[8],
    well_logs_list_3[[38]][1]
  )
)


val_mat_sel <- val_mat_sel[, order(depth_mat_sel[1, ], decreasing = TRUE)]
depth_mat_sel <- depth_mat_sel[, order(depth_mat_sel[1, ], decreasing =
                                         TRUE)]

depth_mat_sel

depth_mat_sel[3, 15] <-47
depth_mat_sel[8, 9] <- 130
depth_mat_sel[5, 14] <- 0 
depth_mat_sel[8, 10] <- 0

depth_mat_sel



Altajme_streched_sel_L_SIL <- Altajme_streched[Altajme_streched[, 2] < depth_mat_sel[8, 1], ]
Altajme_streched_sel_L_SIL <- Altajme_streched_sel_L_SIL[Altajme_streched_sel_L_SIL[, 2] >
                                                           depth_mat_sel[8, 2], ]


Altajme_streched_sel_l_visby <- Altajme_streched[Altajme_streched[, 2] <
                                                   depth_mat_sel[8, 2], ]
Altajme_streched_sel_l_visby <- Altajme_streched_sel_l_visby[Altajme_streched_sel_l_visby[, 2] >
                                                               depth_mat_sel[8, 3], ]


Altajme_streched_sel_u_visby <- Altajme_streched[Altajme_streched[, 2] <
                                                   depth_mat_sel[8, 3], ]
Altajme_streched_sel_u_visby <- Altajme_streched_sel_u_visby[Altajme_streched_sel_u_visby[, 2] >
                                                               depth_mat_sel[8, 4], ]

Altajme_streched_sel_Hogklint <- Altajme_streched[Altajme_streched[, 2] <
                                                    depth_mat_sel[8, 4], ]
Altajme_streched_sel_Hogklint <- Altajme_streched_sel_Hogklint[Altajme_streched_sel_Hogklint[, 2] >
                                                                 depth_mat_sel[8, 5], ]




Altajme_streched_sel_Tofta <- Altajme_streched[Altajme_streched[, 2] < depth_mat_sel[8, 5], ]
Altajme_streched_sel_Tofta <- Altajme_streched_sel_Tofta[Altajme_streched_sel_Tofta[, 2] >
                                                           depth_mat_sel[8, 6], ]

Altajme_streched_sel_Hangvar <- Altajme_streched[Altajme_streched[, 2] <
                                                   depth_mat_sel[8, 6], ]
Altajme_streched_sel_Hangvar <- Altajme_streched_sel_Hangvar[Altajme_streched_sel_Hangvar[, 2] >
                                                               depth_mat_sel[8, 7], ]

Altajme_streched_sel_Slite <- Altajme_streched[Altajme_streched[, 2] < depth_mat_sel[8, 7], ]
Altajme_streched_sel_Slite <- Altajme_streched_sel_Slite[Altajme_streched_sel_Slite[, 2] >
                                                           depth_mat_sel[8, 8], ]

Altajme_streched_sel_Mulde_1 <- Altajme_streched[Altajme_streched[, 2] <
                                                   depth_mat_sel[8, 8], ]
Altajme_streched_sel_Mulde_1 <- Altajme_streched_sel_Mulde_1[Altajme_streched_sel_Mulde_1[, 2] >
                                                               depth_mat_sel[8, 9], ]

Altajme_streched_sel_Mulde_2 <- Altajme_streched[Altajme_streched[, 2] <
                                                   depth_mat_sel[8, 9], ]
Altajme_streched_sel_Mulde_2 <- Altajme_streched_sel_Mulde_2[Altajme_streched_sel_Mulde_2[, 2] >
                                                               depth_mat_sel[8, 10], ]

colnames(Altajme_streched_sel_L_SIL) <- c("A", "B")
colnames(Altajme_streched_sel_l_visby) <- c("A", "B")
colnames(Altajme_streched_sel_u_visby) <- c("A", "B")
colnames(Altajme_streched_sel_Hogklint) <- c("A", "B")
colnames(Altajme_streched_sel_Tofta) <- c("A", "B")
colnames(Altajme_streched_sel_Hangvar) <- c("A", "B")
colnames(Altajme_streched_sel_Slite) <- c("A", "B")
colnames(Altajme_streched_sel_Mulde_1) <- c("A", "B")
colnames(Altajme_streched_sel_Mulde_2) <- c("A", "B")

bottom_Sil <- cbind(c(-5000, val_mat_sel[, 1]), c(depth_mat_sel[1, 1], depth_mat_sel[, 1]))
colnames(bottom_Sil) <- c("A", "B")
bottom_L_visby <- cbind(c(-5000, val_mat_sel[, 2]), c(depth_mat_sel[1, 2], depth_mat_sel[, 2]))
colnames(bottom_L_visby) <- c("A", "B")
bottom_u_visby <- cbind(c(-5000, val_mat_sel[, 3]), c(depth_mat_sel[1, 3], depth_mat_sel[, 3]))
colnames(bottom_u_visby) <- c("A", "B")
top_u_visby <- cbind(c(-5000, val_mat_sel[, 4]), c(depth_mat_sel[1, 4], depth_mat_sel[, 4]))
colnames(top_u_visby) <- c("A", "B")
top_Hogklint <- cbind(c(-5000, val_mat_sel[, 5]), c(depth_mat_sel[1, 5], depth_mat_sel[, 5]))
colnames(top_Hogklint) <- c("A", "B")
top_Tofta <- cbind(c(-5000, val_mat_sel[, 6]), c(depth_mat_sel[1, 6], depth_mat_sel[, 6]))
colnames(top_Tofta) <- c("A", "B")
top_Hangvar <- cbind(c(-5000, val_mat_sel[, 7]), c(depth_mat_sel[1, 7], depth_mat_sel[, 7]))
colnames(top_Hangvar) <- c("A", "B")
top_Slite <- cbind(c(-5000, val_mat_sel[, 8]), c(depth_mat_sel[1, 8], depth_mat_sel[, 8]))
colnames(top_Slite) <- c("A", "B")
top_Mulde_1 <- cbind(c(-5000, val_mat_sel[, 9]), c(depth_mat_sel[1, 9], depth_mat_sel[, 9]))
colnames(top_Mulde_1) <- c("A", "B")
top_Mulde_2 <- cbind(c(-5000, val_mat_sel[, 10]), c(depth_mat_sel[1, 10], depth_mat_sel[, 10]))
colnames(top_Mulde_2) <- c("A", "B")

top_Hemse <- cbind(c(-5000, val_mat_sel[, 11]), c(depth_mat_sel[1, 11], depth_mat_sel[, 11]))
colnames(top_Hemse) <- c("A", "B")
base_Etelheim_2 <- cbind(c(-5000, val_mat_sel[, 12]), c(depth_mat_sel[1, 12], depth_mat_sel[, 12]))
colnames(base_Etelheim_2) <- c("A", "B")
base_Eke <- cbind(c(-5000, val_mat_sel[, 13]), c(depth_mat_sel[1, 13], depth_mat_sel[, 13]))
colnames(base_Eke) <- c("A", "B")
base_Burgsvik <- cbind(c(-5000, val_mat_sel[, 14]), c(depth_mat_sel[1, 14], depth_mat_sel[, 14]))
colnames(base_Burgsvik) <- c("A", "B")
base_Hamra_Sundre <- cbind(c(-5000, val_mat_sel[, 15]), c(depth_mat_sel[1, 15], depth_mat_sel[, 15]))
colnames(base_Eke) <- c("A", "B")


L_Sil <- rbind(bottom_Sil, 
               Altajme_streched_sel_L_SIL[rev(order(Altajme_streched_sel_L_SIL[, 2])), ], 
               bottom_L_visby[rev(order(bottom_L_visby[, 1])), ])

polygon(x = L_Sil[, 1],
        y = L_Sil[, 2],
        col = rgb(0, 0, 1, 0.5))

L_Visby <- rbind(bottom_L_visby,
                 Altajme_streched_sel_l_visby[rev(order(Altajme_streched_sel_l_visby[, 2])), ],
                 bottom_u_visby[rev(order(bottom_u_visby[, 1])), ])
polygon(x = L_Visby[, 1],
        y = L_Visby[, 2],
        col = rgb(0, 0, 1, 0.5))

U_Visby <- rbind(bottom_u_visby, Altajme_streched_sel_u_visby[rev(order(Altajme_streched_sel_u_visby[, 2])), ], top_u_visby[rev(order(top_u_visby[, 1])), ])

polygon(x = U_Visby[, 1],
        y = U_Visby[, 2],
        col = rgb(0, 0, 1, 0.5))

Hogklint <- rbind(top_u_visby, Altajme_streched_sel_Hogklint[rev(order(Altajme_streched_sel_Hogklint[, 2])), ], top_Hogklint[rev(order(top_Hogklint[, 1])), ])
polygon(x = Hogklint[, 1],
        y = Hogklint[, 2],
        col = rgb(0, 0, 1, 0.5))

Tofta <- rbind(top_Hogklint, Altajme_streched_sel_Tofta[rev(order(Altajme_streched_sel_Tofta[, 2])), ], top_Tofta[rev(order(top_Tofta[, 1])), ])
polygon(x = Tofta[, 1],
        y = Tofta[, 2],
        col = rgb(0, 0, 1, 0.5))

Hangvar <- rbind(top_Tofta, Altajme_streched_sel_Hangvar[rev(order(Altajme_streched_sel_Hangvar[, 2])), ], top_Hangvar[rev(order(top_Hangvar[, 1])), ])
polygon(x = Hangvar[, 1],
        y = Hangvar[, 2],
        col = rgb(0, 0, 1, 0.5))


top_Slite <- top_Slite[1:7, ]

Slite <- rbind(top_Hangvar,
               Altajme_streched_sel_Slite[rev(order(Altajme_streched_sel_Slite[, 2])), ],
               top_Slite[rev(order(top_Slite[, 1])), ])

polygon(x = Slite[, 1],
        y = Slite[, 2],
        col = rgb(0, 0, 1, 0.5))


Mulde_1 <- rbind(top_Slite, 
                 Altajme_streched_sel_Mulde_1[rev(order(Altajme_streched_sel_Mulde_1[, 2])), ], 
                 top_Mulde_1[rev(order(top_Mulde_1[, 1])), ])
polygon(x = Mulde_1[, 1],
        y = Mulde_1[, 2],
        col = rgb(0, 0, 1, 0.5))

top_Mulde_2_sel <- top_Mulde_2[rev(order(top_Mulde_2[, 1])), ]
top_Mulde_2_sel <- top_Mulde_2_sel[2:8, ]

Mulde_2 <- rbind(
  top_Mulde_1,
  Altajme_streched_sel_Mulde_2[rev(order(Altajme_streched_sel_Mulde_2[, 2])), ],
  c(Altajme_streched_sel_Mulde_2[1, 1], 0),
  c(Altajme_streched_sel_Mulde_2[1, 1] - 5000, 0),
  top_Mulde_2_sel,
  c(-5000,top_Mulde_2_sel[nrow(top_Mulde_2_sel),2])
)
polygon(x = Mulde_2[, 1],
        y = Mulde_2[, 2],
        col = rgb(0, 0, 1, 0.5))







L_Hemse <- rbind(
  top_Mulde_2[1:8, ],
  c(Altajme_streched_sel_Mulde_2[1, 1] - 5000, 0),
  c(Altajme_streched_sel_Mulde_2[1, 1] - 7250, 0),
  na.omit(top_Hemse[rev(order(top_Hemse[, 1])), ])
)
polygon(x = L_Hemse[, 1],
        y = L_Hemse[, 2],
        col = rgb(0, 1, 0, 0.5))



top_Hemse_sel <- na.omit(top_Hemse[rev(order(top_Hemse[, 1])), ])

Etelheim <- rbind(
  na.omit(top_Hemse),
  c(Altajme_streched_sel_Mulde_2[1, 1] - 7250, 0),
  c(top_Hemse_sel[1, 1], 0),
  na.omit(base_Etelheim_2[rev(order(base_Etelheim_2[, 1])), ])
)
polygon(x = Etelheim[, 1],
        y = Etelheim[, 2],
        col = rgb(0, 0, 1, 0.5))

base_Eke_sel <- na.omit(base_Eke[rev(order(base_Eke[, 1])), ])

Nar <- rbind(
  na.omit(base_Etelheim_2),
  c(top_Hemse_sel[1, 1], 0),
  c(base_Eke_sel[1, 1] + 500, 0),
  na.omit(base_Eke[rev(order(base_Eke[1:6, 1])), ])
)


polygon(x = Nar[, 1],
        y = Nar[, 2],
        col = rgb(0, 0, 1, 0.5))

base_Burgsvik[5, 2] <- base_Burgsvik[5, 2]-6
base_Burgsvik[6,2] <- 40

Eke <- rbind(na.omit(base_Eke[c(1:6), ]),
             c(base_Eke_sel[1, 1] + 500, 0),
             c((22 * 1000) + 500, 0),
             na.omit(base_Burgsvik[rev(order(base_Burgsvik[1:6, 1])), ]))



polygon(x = Eke[, 1],
        y = Eke[, 2],
        col = rgb(0, 0, 1, 0.5))

base_base_Burgsvik_sel <- na.omit(base_Burgsvik[rev(order(base_Burgsvik[, 1])), ])

Burgsvik <- rbind(
  na.omit(base_Burgsvik[1:6, ]),
  c((22 * 1000) + 500, 0),
  c(base_Hamra_Sundre[5, 1], 0),
  na.omit(base_Hamra_Sundre[rev(order(base_Hamra_Sundre[1:4, 1])), ])
)

polygon(x = Burgsvik[, 1],
        y = Burgsvik[, 2],
        col = rgb(0, 0, 1, 0.5))


Hamra_Sundre <- rbind(na.omit(base_Hamra_Sundre[1:4, ]),
                      c(base_Hamra_Sundre[5, 1], 0),
                      c(-5000, 0))

polygon(x = Hamra_Sundre[, 1],
        y = Hamra_Sundre[, 2],
        col = rgb(0, 0, 1, 0.5))


# add  the d13C record to the plot
Altamje_d13C = read.csv("D:/PhD/R/R_results/Altajme_d13C.csv",
                        header = T,
                        sep = ";")

Altamje_d13C <- demean(Altamje_d13C, genplot = F, verbose = T)
Altamje_d13C[, 2]  <- Altamje_d13C[, 2] - min(Altamje_d13C[, 2])
Altamje_d13C[, 2]  <- Altamje_d13C[, 2] / max(Altamje_d13C[, 2])
Altamje_d13C <- iso(Altamje_d13C,
                    xmax = 324,
                    xmin = 0,
                    genplot = FALSE)
lines(
  x = (Altamje_d13C[, 2] * 4000) + 62000,
  y = Altamje_d13C[, 1],
  col = "green"
)
}

#plot ages  litholog reference and st.Sutarve-2018 and Altajme ######
{
Llandovery_col <- geo_col(name = "Llandovery")
Wenlock_col <- geo_col(name = "Wenlock")
Telychian_col <- geo_col(name = "Telychian")
Sheinwoodian_col <- geo_col(name = "Sheinwoodian")
Homerian_col <- geo_col(name = "Homerian")
Rhuddanian_col <- geo_col(name = "Rhuddanian")
Aeronian_col <- geo_col(name = "Aeronian")
Ludfordian_col <- geo_col(name = "Ludfordian")
Ludlow_col <- geo_col(name = "Ludlow")
Gorstian_col <- geo_col(name = "Gorstian")
Ord_col <- geo_col(name = "Ordovician")


L_Sil_col <- "#E2A9CD"
L_Visby_col <- "#C279A9"
U_Visby_col <- "#915E81"
Hogklint_col <- "#4065A9"
Tofta_col <- "#49A6CA"
Hangvar_col  <-  "#0C6867"
Slite_col   <- "#1D9A38"
Frojel_col  <- "#999999ff"
Halla_Klinteberg_col  <- "#966919"
L_Hemse_col  <- "#8E8B52"
Etelheim_col  <- "#CDDE94"
Nar_col  <- "#CA970B"
Eke_col  <- "#F49820"
Brugsvik_col <- "#EFE41B"
Hamra_Sundre_col <- "#FCCB68"



base_silurian


tie_points
tie_points_2


shift_down <- 60
graphics.off()
plot(
  y = (a[tie_points, 1] / 10) + shift_down,
  x = a[tie_points, 2] - 1 ,
  col = colors,
  pch = 19,
  cex = 0.75
)

a_sort <- a[tie_points[order(a[tie_points, 1]),],1]/10
a_sort


graphics.off()
plot(
  y = c(a[, 1] / 10)+shift_down,
  x = a[, 2] ,type="l",
  ylim=c(450,300)
)

depth_mat_sel

readjusted_matrix_ord

ref_base_hogklint <-base_hogklint[1]
#base Silurian


ref_base_silurian <- readjusted_matrix_ord$base_silurian[1]

ref_base_l_Visby <-readjusted_matrix_ord$base_L_visby[1]

#base upper Visby
ref_base_u_Visby <- readjusted_matrix_ord$base_U_visby[1]

#base Tofta
ref_base_Tofta <- readjusted_matrix_ord$base_Tofta[1]

#base Hangvar
ref_base_Hangvar <- readjusted_matrix_ord$base_Hangvar[1] 
#base Slite
ref_base_Slite <- readjusted_matrix_ord$base_Slite[1]
#base Halla # start Mulde
ref_base_Halla <- readjusted_matrix_ord$base_Klinteberg[1]
#base hemse near end Mulde
ref_base_Hemse <- readjusted_matrix_ord$base_Hemse[1]
#base Etelhem 
ref_base_Etelhem <- readjusted_matrix_ord$base_Etelhem[1]
#base Nar 
ref_base_Nar <-readjusted_matrix_ord$base_Nar[1]
#base EKE 
ref_base_EKE <- readjusted_matrix_ord$base_EKE[1]
#base Burgsvik
ref_base_Burgsvik <- readjusted_matrix_ord$base_Burgsvik[1]
#base Hamra Sundre
ref_base_Ham_Sund <- readjusted_matrix_ord$base_Hamra_Sundre[1]
#base Frojel
ref_base_Frojel <- readjusted_matrix_ord$base_Frojel[1]
#base klinteberg
ref_base_Klinteberg <- readjusted_matrix_ord$base_Klinteberg[1]


graphics.off()
{
layout.matrix <-
  matrix(c(1:4), nrow = 1, ncol = 4)

layout(layout.matrix,
       heights = c(1),
       widths = c(1,1,1,2))
par(mar = c(4, 4, 4, 0))

ylims <- c(max((a[, 1] / 10))+2.5+shift_down,0)


plot(
  x = c(0, 1),
  y = c(0, max(ylims)),
  col = "white",
  xlab = "",
  ylab = "normalized depth (m)",
  xaxt = "n",
  xaxs = "i",
  yaxs = "i",
  ylim = ylims,
) # Draw empty plot

polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_silurian,
    ref_base_silurian,
    max(ylims),
    max(ylims)
  ),
  col = Ord_col
)

polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_silurian,
    ref_base_silurian,
    ref_base_l_Visby,
    ref_base_l_Visby
  ),
  col = Aeronian_col
)


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_l_Visby,
    ref_base_l_Visby,
    ref_base_u_Visby,
    ref_base_u_Visby
    
  ),
  col = Llandovery_col
)




polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_u_Visby,
    ref_base_u_Visby,
    ref_base_Hemse,
    ref_base_Hemse
  ),
  col = Wenlock_col
)



polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Hemse,
    ref_base_Hemse,
    0,0
  ),
  col = Ludfordian_col
)



par(mar = c(4, 0, 4, 0))

plot(
  x = c(0, 1),
  y = c(0, max(ylims)),
  col = "white",
  xlab = "",
  ylab = "depth (m)",
  xaxt = "n",
  xaxs = "i",
  yaxs = "i",
  yaxt="n",
  ylim = ylims,
) # Draw empty plot



polygon(
  x = c(0, 1, 1, 0),
  y = c(
    527.5,
    527.5,
    ref_base_silurian,
    ref_base_silurian
  ),
  col = Ord_col
)

polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_silurian,
    ref_base_silurian,
    ref_base_l_Visby,
    ref_base_l_Visby
  ),
  col = Aeronian_col
)


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_l_Visby,
    ref_base_l_Visby,
    ref_base_u_Visby+1.5,
    ref_base_u_Visby+1.5
    
  ),
  col = Llandovery_col
)


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_u_Visby+1.5,
    ref_base_u_Visby+1.5,
    (ref_base_Slite+ref_base_Halla)/2,
    (ref_base_Slite+ref_base_Halla)/2
  ),
  col = Sheinwoodian_col
)


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    (ref_base_Slite+ref_base_Halla)/2,
    (ref_base_Slite+ref_base_Halla)/2,
    ref_base_Hemse,
    ref_base_Hemse
  ),
  col = Homerian_col
)
polygon(
  x = c(0, 1, 1, 0),
  y = c(
    
    ref_base_Hemse,
    ref_base_Hemse,
    ref_base_Etelhem,
    ref_base_Etelhem
  ),
  col = Gorstian_col
)

polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Etelhem,
    ref_base_Etelhem,
    0,
    0
  ),
  col = Ludfordian_col
)


par(mar = c(4, 0, 4, 0))

plot(
  x = c(0, 1),
  y = c(0, max(ylims)),
  col = "white",
  xlab = "",
  ylab = "depth (m)",
  xaxt = "n",
  xaxs = "i",
  yaxs = "i",
  yaxt="n",
  ylim = ylims,
) # Draw empty plot


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    527.5,
    527.5,
    ref_base_silurian,
    ref_base_silurian
  )
)

text(x=0.5,y=prestreched_data[[1]][1, 3]/3+(525+(ref_base_silurian))/2,label="Ord")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_silurian,
    ref_base_silurian,
    ref_base_l_Visby,
    ref_base_l_Visby
  ))


text(x=0.5,y=(ref_base_silurian+(ref_base_l_Visby))/2,
     label="L-Sil")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_l_Visby,
    ref_base_l_Visby,
    ref_base_u_Visby,
    ref_base_u_Visby
    
  )
)

text(x=0.5,y=(ref_base_u_Visby+(ref_base_l_Visby))/2,
     label="L-Visby")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_u_Visby,
    ref_base_u_Visby,
    ref_base_hogklint,
    ref_base_hogklint
  ),
  col = 
)

text(x=0.5,y=(ref_base_u_Visby+(ref_base_hogklint))/2,
     label="U-Visby")



polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_hogklint,
    ref_base_hogklint,
    ref_base_Tofta,
    ref_base_Tofta
  ),
  col = 
)

text(x=0.5,y=(ref_base_Tofta+(ref_base_hogklint))/2,
     label="Hogklint")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Tofta,
    ref_base_Tofta,
    ref_base_Hangvar,
    ref_base_Hangvar
  ),
  col = 
)

text(x=0.5,y=(ref_base_Tofta+(ref_base_Hangvar))/2,
     label="Tofta")



polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Hangvar,
    ref_base_Hangvar,
    ref_base_Slite,
    ref_base_Slite
  ),
  col = 
)

text(x=0.5,y=(ref_base_Slite+(ref_base_Hangvar))/2,
     label="Hangvar")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Slite,
    ref_base_Slite,
    ref_base_Halla,
    ref_base_Halla
  ),
  col = 
)

text(x=0.5,y=(ref_base_Slite+(ref_base_Halla))/2,
     label="Slite")



polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Frojel,
    ref_base_Frojel,
    ref_base_Klinteberg,
    ref_base_Klinteberg
  ),
  col = 
)

text(x=0.5,y=(ref_base_Frojel+(ref_base_Klinteberg))/2,
     label="Frojel")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Klinteberg,
    ref_base_Klinteberg,
    ref_base_Hemse,
    ref_base_Hemse
  ),
  col = 
)

text(x=0.5,y=(ref_base_Klinteberg+(ref_base_Hemse))/2,
     label="Klinteberg Fm /Halla Fm")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Hemse,
    ref_base_Hemse,
    ref_base_Etelhem,
    ref_base_Etelhem
  ),
  col = 
)

text(x=0.5,y=(ref_base_Hemse+(ref_base_Etelhem))/2,
     label="L-Hemse Marl Fm")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Etelhem,
    ref_base_Etelhem,
    ref_base_Nar,
    ref_base_Nar
  ),
  col = 
)

text(x=0.5,y=(ref_base_Etelhem+(ref_base_Nar))/2,
     label="Etelhem Fm")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Nar,
    ref_base_Nar,
    ref_base_EKE,
    ref_base_EKE
  ),
  col = 
)

text(x=0.5,y=(ref_base_Nar+(ref_base_EKE))/2,
     label="Nar Fm")

polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_EKE,
    ref_base_EKE,
    ref_base_Burgsvik,
    ref_base_Burgsvik
  ),
  col = 
)

text(x=0.5,y=(ref_base_Burgsvik+(ref_base_EKE))/2,
     label="Eke Fm")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Burgsvik,
    ref_base_Burgsvik,
    ref_base_Ham_Sund,
    ref_base_Ham_Sund
  ),
  col = 
)

text(x=0.5,y=(ref_base_Burgsvik+(ref_base_Ham_Sund))/2,
     label="Burgsvik Fm")


polygon(
  x = c(0, 1, 1, 0),
  y = c(
    ref_base_Ham_Sund,
    ref_base_Ham_Sund,
    0,
    0
  ),
  col = 
)

text(x=0.5,y=(50+(ref_base_Ham_Sund))/2,
     label="Hamra/Sundre Fm")

par(mar = c(4, 0, 4, 4))

plot(a[, 2], (a[, 1] / 10)+shift_down ,xlab="",yaxs = "i",xaxt="n",
     col = "black", lwd = 1,type="l",xlim=c(-1,3.75),ylim=ylims,yaxt="n",xaxs="i")


polygon(
  x = c(-1,3.75,3.75,2.5,2,0.75,0.25, -1),
  y = c(
    527.5,
    527.5,
    L_Sil[9,2],
    L_Sil[9,2],
    L_Sil[4,2],
    L_Sil[4,2],
    ref_base_silurian,
    ref_base_silurian
  ),col=Ord_col
)

polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),
  y = c(
    ref_base_silurian,
    ref_base_silurian,
    L_Sil[4,2],
    L_Sil[4,2],
    L_Sil[9,2],
    L_Sil[9,2],
    L_Visby[9,2],
    L_Visby[9,2],
    L_Visby[4,2],
    L_Visby[4,2],
    ref_base_l_Visby,
    ref_base_l_Visby
    
  ),col=L_Sil_col)


polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),
  y = c(
    ref_base_l_Visby,
    ref_base_l_Visby,
    L_Visby[4,2],
    L_Visby[4,2],
    L_Visby[9,2],
    L_Visby[9,2],
    U_Visby[9,2],
    U_Visby[9,2],
    U_Visby[4,2],
    U_Visby[4,2],
    ref_base_u_Visby,
    ref_base_u_Visby
    
  ),col=L_Visby_col)


polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),  y = c(
    ref_base_u_Visby,
    ref_base_u_Visby,
    U_Visby[4,2],
    U_Visby[4,2],
    U_Visby[9,2],
    U_Visby[9,2],
    Hogklint[9,2],
    Hogklint[9,2],
    Hogklint[4,2],
    Hogklint[4,2],
    ref_base_hogklint,
    ref_base_hogklint
    
  ),col=U_Visby_col)




polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),  y = c(
    ref_base_hogklint,
    ref_base_hogklint,
    Hogklint[4,2],
    Hogklint[4,2],
    Hogklint[9,2],
    Hogklint[9,2],
    Tofta[9,2],
    Tofta[9,2],
    Tofta[4,2],
    Tofta[4,2],
    ref_base_Tofta,
    ref_base_Tofta
    
  ),col=Hogklint_col)

polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),  y = c(
    ref_base_Tofta,
    ref_base_Tofta,
    Tofta[4,2],
    Tofta[4,2],
    Tofta[9,2],
    Tofta[9,2],
    Hangvar[9,2],
    Hangvar[9,2],
    Hangvar[4,2],
    Hangvar[4,2],
    ref_base_Hangvar,
    ref_base_Hangvar
    
  ),col=Tofta_col)


polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),  y = c(
    ref_base_Hangvar,
    ref_base_Hangvar,
    Hangvar[4,2],
    Hangvar[4,2],
    Hangvar[9,2],
    Hangvar[9,2],
    Slite[9,2],
    Slite[9,2],
    Slite[4,2],
    Slite[4,2],
    ref_base_Slite,
    ref_base_Slite
    
  ),col=Hangvar_col)


polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),  y = c(
    ref_base_Slite,
    ref_base_Slite,
    Slite[4,2],
    Slite[4,2],
    Slite[9,2],
    Slite[9,2],
    Mulde_1[9,2],
    Mulde_1[9,2],
    Mulde_1[4,2],
    Mulde_1[4,2],
    ref_base_Frojel,
    ref_base_Frojel
    
  ),col=Slite_col)

polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),  y = c(
    ref_base_Frojel,
    ref_base_Frojel,
    Mulde_1[4,2],
    Mulde_1[4,2],
    Mulde_1[9,2],
    Mulde_1[9,2],
    Mulde_2[9,2],
    Mulde_2[9,2],
    Mulde_2[4,2],
    Mulde_2[4,2],
    ref_base_Klinteberg,
    ref_base_Klinteberg
    
  ),col=Frojel_col)

polygon(
  x = c(-1,0.25,0.75,2,2.5,3.75,3.75,2.5,
        2,0.75,0.25, -1),  y = c(
    ref_base_Klinteberg,
    ref_base_Klinteberg,
    Mulde_2[4,2],
    Mulde_2[4,2],
    Mulde_2[9,2],
    Mulde_2[9,2],0,
    0,
    L_Hemse[4,2],
    L_Hemse[4,2],
    ref_base_Hemse,
    ref_base_Hemse
    
  ),col=Halla_Klinteberg_col)


polygon(
  x = c(-1,0.25,0.75,2,2.5,2,0.75,0.25, -1),
  y = c(
    ref_base_Hemse,
    ref_base_Hemse,
    L_Hemse[4,2],
    L_Hemse[4,2],
    0,
    Etelheim[4,2],
    Etelheim[4,2],
    ref_base_Etelhem,
    ref_base_Etelhem
    
  ),col=L_Hemse_col)

polygon(
  x = c(-1,0.25,0.75,2,2.5,2,0.75,0.25, -1),
  y = c(
    ref_base_Etelhem,
    ref_base_Etelhem,
    Etelheim[4,2],
    Etelheim[4,2],
    0,
    Nar[4,2]-4,
    Nar[4,2]-4,
    ref_base_Nar,
    ref_base_Nar
    
  ),col=Etelheim_col)


polygon(
  x = c(-1,0.25,0.75,2,2.5,2,0.75,0.25, -1),
  y = c(
    ref_base_Nar,
    ref_base_Nar,
    Nar[4,2]-4,
    Nar[4,2]-4,
    0,
    Eke[4,2],
    Eke[4,2],
    ref_base_EKE,
    ref_base_EKE
    
  ),col=Nar_col)


polygon(
  x = c(-1,0.25,0.75,2,2.5,2,0.75,0.25, -1),
  y = c(
    ref_base_EKE,
    ref_base_EKE,
    Eke[4,2],
    Eke[4,2],
    0,
    Burgsvik[4,2],
    Burgsvik[4,2],
    ref_base_Burgsvik,
    ref_base_Burgsvik
    
  ),col=Eke_col)

polygon(
  x = c(-1,0.25,0.75,2,2.5,2,0.75,0.25, -1),
  y = c(
    ref_base_Burgsvik,
    ref_base_Burgsvik,
    Burgsvik[4,2],
    Burgsvik[4,2],
    0,
    Hamra_Sundre[4,2],
    Hamra_Sundre[4,2],
    ref_base_Ham_Sund,
    ref_base_Ham_Sund
    
  ),col=Brugsvik_col)

polygon(
  x = c(-1,0.25,0.75,2,2.5,2,0.75,0.25, -1),
  y = c(
    ref_base_Ham_Sund,
    ref_base_Ham_Sund,
    Hamra_Sundre[4,2],
    Hamra_Sundre[4,2],
    0,
    0,
    0,
    0,
    0    
  ),col=Hamra_Sundre_col)

lines(a[, 2], (a[, 1] / 10)+shift_down,lwd=1)
lines((prestreched_data_2$StSutarve_2018[,2])+1.5,prestreched_data_2$StSutarve_2018[,1],lwd=1)
lines((prestreched_data_2$Altajme[,2])+3,prestreched_data_2$Altajme[,1],lwd=1)

}
}





 #Plot some maps ####
well_tops_rotated_3  <- well_tops_rotated
GROTLLINGBO_2 <- c(1660138,6335746,"Grotlingbo_2",rep(1,7))
GROTLLINGBO_1 <- c(1660000,6335550,"Grotlingbo_1",rep(1,7))
GROTLLINGBO <- rbind(GROTLLINGBO_1,GROTLLINGBO_2)
GROTLLINGBO <- as.data.frame(GROTLLINGBO)
colnames(GROTLLINGBO) <- colnames(well_tops_rotated_3)
well_tops_rotated_3 <- rbind(well_tops_rotated_3,GROTLLINGBO)

well_tops_rotated_3[,1] <- as.numeric(well_tops_rotated_3[,1])
well_tops_rotated_3[,2] <- as.numeric(well_tops_rotated_3[,2])

#load a shapefile
shapename <- st_read('D:/Phd/documents/Altajme/DTW_paper/well_logs/SI_4_logs_used/maps')
local_path <- "DTW_Gotland_Local/SI_4_logs_used/maps/gotland_wells.shp" # Replace 'gotland_wells.shp' with the actual shapefile name

# 3. READ THE SHAPEFILE
shapename <- st_read(local_path)


library(sf)

# 1. Define the necessary files and their base URL
base_url <- "https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/main/SI_4_logs_used/maps/gotland_geology_2"
extensions <- c(".shp", ".shx", ".dbf", ".prj", ".cpg")

# 2. Create a temporary directory to save the files
temp_dir <- tempdir()
local_shp_path <- file.path(temp_dir, "gotland_geology_2.shp")

# 3. Download each component file
for (ext in extensions) {
  file_url <- paste0(base_url, ext)
  local_file_path <- file.path(temp_dir, paste0("gotland_geology_2", ext))
  
  # Use download.file to fetch the raw file content
  # mode="wb" (write binary) is crucial for non-text files like .shp and .dbf
  download.file(url = file_url, destfile = local_file_path, mode = "wb")
}

# 4. Read the shapefile from the local temporary directory
# st_read now sees all related files locally, which avoids the GDAL/HTTPS issues.
shapename <- st_read(local_shp_path)

Visby <- shapename[shapename$id == 1, ]
Hogklint <- shapename[shapename$id == 2, ]
Tofta <- shapename[shapename$id == 3, ]
Hangvar <- shapename[shapename$id == 4, ]
Slite <- shapename[shapename$id == 5, ]
Frojel <- shapename[shapename$id == 6, ]
Halla <- shapename[shapename$id == 7, ]
L_Hemse <- shapename[shapename$id == 8, ]
Etelheim <- shapename[shapename$id == 9, ]
Klinteberg <- shapename[shapename$id == 10, ]
Nar <- shapename[shapename$id == 11, ]
Eke <- shapename[shapename$id == 12, ]
Brugsvik <- shapename[shapename$id == 13, ]
Hamra <- shapename[shapename$id == 14, ]
Visby$id <- c("Visby Fm")

#plot the geological map of Gotland
View(well_tops_rotated_3)
well_tops_rotated_3[,3]

pl <- ggplot() +
  geom_sf(data = shapename) +
  coord_sf(default_crs = sf::st_crs(x = 2400), clip = "off")
pl
pl + geom_sf(data = Visby,
             fill  = rgb(195, 121, 170, 255, maxColorValue = 255)) +
  geom_sf(data = Hogklint,
          fill  = rgb(63, 100, 168, 255, maxColorValue = 255)) +
  geom_sf(data = Tofta,
          fill  = rgb(71, 165, 201, 255, maxColorValue = 255)) +
  geom_sf(data = Hangvar, fill  = "#076766ff") +
  geom_sf(data = Slite, fill  = "#1d9b38ff") +
  geom_sf(data = Frojel, fill  = "#999999ff") +
  geom_sf(data = Halla, fill  = "#976818ff") +
  geom_sf(data = L_Hemse, fill  = "#976818ff") +
  geom_sf(data = Etelheim, fill  = "#8f8d54ff") +
  geom_sf(data = Klinteberg, fill  = "#ccde92ff") +
  geom_sf(data = Nar, fill  = "#cb9807ff") +
  geom_sf(data = Eke, fill  = "#f4981dff") +
  geom_sf(data = Brugsvik, fill  = "#f1e513ff") +
  geom_sf(data = Hamra, fill  = "#fdcc68ff") +
  geom_point(aes(x = well_tops_rotated_3$X, y = well_tops_rotated_3$Y))+  
  geom_label_repel(
    aes(x = well_tops_rotated_3[, 1], y = well_tops_rotated_3[, 2], label = well_tops_rotated_3[, 3]),
    box.padding   = 0.2,
    point.padding = 0.2,
    size = 1,
    max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
    segment.color = 'red',
    segment.size = 0.1,   # ensure segment line is visible
    segment.alpha = 1     # ensure it's not transparent
  ) + annotation_scale(location = "tl") +
  annotation_north_arrow(location = "br",
                         which_north = "true",
                         style = north_arrow_minimal())


View(well_tops_rotated_3)

#plot regional overview map

europe <- ne_countries(
  continent = "Europe",
  scale = 'medium',
  type = 'map_units',
  returnclass = 'sf'
)

your_map <- ggplot(europe) +
  geom_sf() +
  xlim(c(15, 30)) +
  ylim(c(50, 62)) +
  theme_linedraw() + geom_rect(
    aes(
      xmin = 17.6,
      xmax = 19.5,
      ymin = 56.8,
      ymax = 58.1
    ),
    fill = "NA",
    col = "black"
  ) + annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl",
                         which_north = "true",
                         style = north_arrow_minimal())

your_map


your_map <- ggplot(europe) +
  geom_sf() +
  xlim(c(10, 30)) +
  ylim(c(52, 65)) +
  theme_linedraw() + geom_rect(
    aes(
      xmin = 17.6,
      xmax = 19.5,
      ymin = 56.8,
      ymax = 58.1
    ),
    fill = "NA",
    col = "black"
  ) + annotation_scale(location = "bl") +
  annotation_north_arrow(location = "tl",
                         which_north = "true",
                         style = north_arrow_minimal())

your_map



# SI 2. find the base of the Mulde Event and correlate it ####
library(lasr)

# Define the exact file path and construct the raw URL
file_path_in_repo <- "SI_4_logs_used/logs_SI_2/GROTLINGBO-1/GROTLINGBO-1_FDC_Sc23427-2851.las"
las_url <- paste0("https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/main/", file_path_in_repo)
temp_file <- tempfile(fileext = ".las")
download.file(url = las_url, destfile = temp_file, mode = "wb")
las <- read.las(temp_file)
GB_1_log <- as.data.frame(las$log$log.1$data)  



file_path_in_repo <- "SI_4_logs_used/logs_SI_2/GRÖTLINGBO-2/Sc23427-2219.las"
las_url <- paste0("https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/main/", file_path_in_repo)
temp_file <- tempfile(fileext = ".las")
download.file(url = las_url, destfile = temp_file, mode = "wb")
las <- read.las(temp_file)
GB_2_log <- as.data.frame(las$log$data)  


Kauparve_1_RES <- read.csv("https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/main/SI_4_logs_used/logs_SI_2/Kauparve_1_RES.csv")
Kauparve_1_RES <- Kauparve_1_RES[, c(2, 3)]
Kauparve_1_RES <- sortNave(Kauparve_1_RES)
Kauparve_1_RES <- Kauparve_1_RES[Kauparve_1_RES[, 1] < 425, ]

Grotlingbo_2_Res <- read.csv("https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/main/SI_4_logs_used/logs_SI_2/Grotlingbo_2_Res_2.csv")
Grotlingbo_2_Res <- Grotlingbo_2_Res[, c(2, 3)]
Grotlingbo_2_Res <- sortNave(Grotlingbo_2_Res)
Grotlingbo_2_Res <- Grotlingbo_2_Res[Grotlingbo_2_Res[, 1] < 425, ]
Grotlingbo_2_Res <- linterp(Grotlingbo_2_Res, dt = 0.1)

Kauparve_1_GR <- read.csv("https://raw.githubusercontent.com/stratigraphy/DTW-Gotland/main/SI_4_logs_used/OPAB_Digitized/KAUPARVE_1_GR.TXT",sep="\t")
Kauparve_1_GR <- sortNave(Kauparve_1_GR)
Kauparve_1_GR <- Kauparve_1_GR[Kauparve_1_GR[, 1] < 425, ]
Kauparve_1_GR <- linterp(Kauparve_1_GR, dt = 0.1)

graphics.off()

layout.matrix <-
  matrix(c(1:4), nrow = 1, ncol = 4)

layout(layout.matrix,
       heights = c(1, 1, 1,1),
       widths = c(1, 1, 1,1))

par(mar = c(4, 4, 4, 0.5))

plot(
  Kauparve_1_GR[, 2],
  Kauparve_1_GR[, 1],
  type = "l",
  ylim = c(425, 0),
  main = "Kauparve-1 Gamma-ray",
  ylab = "depth (m)",
  xlab = "Gamma-ray"
)
abline(h = 267, lty = 3, col = "red")
abline(h = 239, lty = 3, col = "red")

plot(
  -log(Kauparve_1_RES[, 2]),
  Kauparve_1_RES[, 1],
  type = "l"
  ,
  ylim = c(425, 0),
  main = "Kauparve-1 Resistivity",
  ylab = "depth (m)",
  xlab = "-log(resistibity)"
)
abline(h = 267, lty = 3, col = "red")
abline(h = 239, lty = 3, col = "red")


plot(
  -log(Grotlingbo_2_Res[, 2]),
  Grotlingbo_2_Res[, 1],
  type = "l"
  ,
  ylim = c(425, 0),
  main = "Grotlingbo-2 Resistivity",
  ylab = "depth (m)",
  xlab = "-log(resistibity)"
)
lines(-log(GB_2_log[,2]),GB_2_log[,1],col="red")

abline(h = 255, lty = 3, col = "red")
abline(h = 227, lty = 3, col = "red")

plot(
  GB_1_log[, 3],
  GB_1_log[, 1],
  type = "l"
  ,log="x",
  ylim = c(425, 0),
  xlim=c(2.65,2.35),
  main = "Grotlingbo-1 Density",
  ylab = "depth (m)",
  xlab = "density (gr/cm^3)")

dens_mean <- mwStats(cbind(GB_1_log[, 1],GB_1_log[, 3]),genplot=F, win=2)

lines(dens_mean[,2],dens_mean[,1],col="green",lwd=2)

abline(h = 255, lty = 3, col = "red")
abline(h = 227, lty = 3, col = "red")


