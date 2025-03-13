# Load data

NTEMS <- read.csv("Output/Tabular Data/NTEMS_with_bird_data.csv")
LULC <- read.csv("Output/Tabular Data/LULC_with_bird_data.csv")

# NTEMS covariate histograms

# Set up the layout
par(mfrow=c(1,3))
par(mar=c(4, 4, 2, 1))

# Calculate ranges
prop_con_range <- range(NTEMS$prop_con_1, NTEMS$prop_con_2, NTEMS$prop_con_3, na.rm = TRUE)
clumpy_range <- range(NTEMS$clumpy_1, NTEMS$clumpy_2, NTEMS$clumpy_3, na.rm = TRUE)
age_mn_range <- range(NTEMS$age_mn_1, NTEMS$age_mn_2, NTEMS$age_mn_3, na.rm = TRUE)
age_md_range <- range(NTEMS$age_md_1, NTEMS$age_md_2, NTEMS$age_md_3, na.rm = TRUE)

# Function to create histogram with adjusted y-axis, consistent x-axis, larger title, and color
create_hist <- function(data, title, xlabel, xlim_range, color) {
  # Calculate breaks based on the data range, not age_mn_range
  num_breaks <- 15 # You can adjust this if needed
  breaks <- seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE), length.out = num_breaks + 1)  
  
  h <- hist(data, breaks = breaks, plot = FALSE)  # Use the calculated breaks
  ylim <- c(0, max(h$counts) * 1.1)  # Set y-axis limits with padding
  hist(data, main = title, xlab = xlabel, breaks = breaks,  # Use breaks here as well
       ylim = ylim, xlim = xlim_range, col = color,  
       cex.main = 1.5)  # Plot with adjustments and color
}

# Create histograms for prop_con variables (green)
create_hist(NTEMS$prop_con_1, "prop_con_1", "prop_con_1", prop_con_range, "forestgreen")
create_hist(NTEMS$prop_con_2, "prop_con_2", "prop_con_2", prop_con_range, "forestgreen")
create_hist(NTEMS$prop_con_3, "prop_con_3", "prop_con_3", prop_con_range, "forestgreen")

# Create histograms for clumpy variables (light red)
create_hist(NTEMS$clumpy_1, "clumpy_1", "clumpy_1", clumpy_range, "lightcoral")  
create_hist(NTEMS$clumpy_2, "clumpy_2", "clumpy_2", clumpy_range, "lightcoral")
create_hist(NTEMS$clumpy_3, "clumpy_3", "clumpy_3", clumpy_range, "lightcoral")

# Create histograms for age_mn variables (blue)
create_hist(NTEMS$age_mn_1, "age_mn_1", "age_mn_1", age_mn_range, "lightblue")
create_hist(NTEMS$age_mn_2, "age_mn_2", "age_mn_2", age_mn_range, "lightblue")
create_hist(NTEMS$age_mn_3, "age_mn_3", "age_mn_3", age_mn_range, "lightblue")

# Reset layout
par(mfrow=c(1,1))




####### LULC covariate histograms

# Set up the layout
par(mfrow=c(1,3))
par(mar=c(4, 4, 2, 1))

prop_con_range <- range(LULC$prop_con_1, LULC$prop_con_2, LULC$prop_con_3, na.rm = TRUE)
clumpy_range <- range(LULC$clumpy_1, LULC$clumpy_2, LULC$clumpy_3, na.rm = TRUE)
age_mn_range <- range(LULC$age_mn_1, LULC$age_mn_2, LULC$age_mn_3, na.rm = TRUE)
age_md_range <- range(LULC$age_md_1, LULC$age_md_2, LULC$age_md_3, na.rm = TRUE)

# Function to create histogram with adjusted y-axis, consistent x-axis, larger title, and color
create_hist <- function(data, title, xlabel, xlim_range, color) {
  # Calculate breaks based on the data range, not age_mn_range
  num_breaks <- 15 # You can adjust this if needed
  breaks <- seq(min(data, na.rm = TRUE), max(data, na.rm = TRUE), length.out = num_breaks + 1) 
  
  h <- hist(data, breaks = breaks, plot = FALSE)  # Use the calculated breaks
  ylim <- c(0, max(h$counts) * 1.1)  # Set y-axis limits with padding
  hist(data, main = title, xlab = xlabel, breaks = breaks,  # Use breaks here as well
       ylim = ylim, xlim = xlim_range, col = color, 
       cex.main = 1.5)  # Plot with adjustments and color
}

# Create histograms for prop_con variables (green)
create_hist(LULC$prop_con_1, "prop_con_1", "prop_con_1", prop_con_range, "forestgreen")
create_hist(LULC$prop_con_2, "prop_con_2", "prop_con_2", prop_con_range, "forestgreen")
create_hist(LULC$prop_con_3, "prop_con_3", "prop_con_3", prop_con_range, "forestgreen")

# Create histograms for clumpy variables (light red)
create_hist(LULC$clumpy_1, "clumpy_1", "clumpy_1", clumpy_range, "lightcoral") 
create_hist(LULC$clumpy_2, "clumpy_2", "clumpy_2", clumpy_range, "lightcoral")
create_hist(LULC$clumpy_3, "clumpy_3", "clumpy_3", clumpy_range, "lightcoral")

# Create histograms for age_mn variables (blue)
create_hist(LULC$age_mn_1, "age_mn_1", "age_mn_1", age_mn_range, "lightblue")
create_hist(LULC$age_mn_2, "age_mn_2", "age_mn_2", age_mn_range, "lightblue")
create_hist(LULC$age_mn_3, "age_mn_3", "age_mn_3", age_mn_range, "lightblue")

# Reset layout
par(mfrow=c(1,1))

#### Create histogram of response 
hist(NTEMS$BTNW)
hist(NTEMS$BBWA)
hist(NTEMS$TEWA)
hist(LULC$BTNW)
hist(LULC$BBWA)
hist(LULC$TEWA)



############## Plot response as a function of NTEMS covariates

### Prop_con_1 vs response
x <- NTEMS$prop_con_1
y <- NTEMS$BTNW
plot(x, y, main="prop_con_1 vs BTNW", 
     xlab="prop conifer 150", ylab="BTNW count", 
     pch=19, col="blue")

y <- NTEMS$BBWA
plot(x, y, main="prop_con_1 vs BBWA", 
     xlab="prop conifer 150", ylab="BBWA count", 
     pch=19, col="red")  

y <- NTEMS$TEWA
plot(x, y, main="prop_con_1 vs TEWA", 
     xlab="prop conifer 150", ylab="TEWA count", 
     pch=19, col="forestgreen")  



### Clumpy_1 vs response
x <- NTEMS$clumpy_1
y <- NTEMS$BTNW
plot(x, y, main="clumpy_1 vs BTNW", 
     xlab="clumpy_1", ylab="BTNW count", 
     pch=19, col="blue")  

y <- NTEMS$BBWA
plot(x, y, main="clumpy_1 vs BBWA", 
     xlab="clumpy_1", ylab="BBWA count", 
     pch=19, col="red")

y <- NTEMS$TEWA
plot(x, y, main="clumpy_1 vs TEWA", 
     xlab="clumpy_1", ylab="TEWA count", 
     pch=19, col="forestgreen")  



### Age_mn_1 vs response
x <- NTEMS$age_mn_1
y <- NTEMS$BTNW
plot(x, y, main="age_mn_1 vs BTNW", 
     xlab="age_mn_1", ylab="BTNW count", 
     pch=19, col="blue")

y <- NTEMS$BBWA
plot(x, y, main="age_mn_1 vs BBWA", 
     xlab="age_mn_1", ylab="BBWA count", 
     pch=19, col="red")  

y <- NTEMS$TEWA
plot(x, y, main="age_mn_1 vs TEWA", 
     xlab="age_mn_1", ylab="TEWA count", 
     pch=19, col="forestgreen")  















############## Plot response as a function of LULC covariates


### Prop_con_1 vs response
x <- LULC$prop_con_1
y <- LULC$BTNW
plot(x, y, main="prop_con_1 vs BTNW", 
     xlab="prop conifer 150", ylab="BTNW count", 
     pch=19, col="blue")

y <- LULC$BBWA
plot(x, y, main="prop_con_1 vs BBWA", 
     xlab="prop conifer 150", ylab="BBWA count", 
     pch=19, col="red")  

y <- LULC$TEWA
plot(x, y, main="prop_con_1 vs TEWA", 
     xlab="prop conifer 150", ylab="TEWA count", 
     pch=19, col="forestgreen")  



### Clumpy_1 vs response
x <- LULC$clumpy_1
y <- LULC$BTNW
plot(x, y, main="clumpy_1 vs BTNW", 
     xlab="clumpy_1", ylab="BTNW count", 
     pch=19, col="blue")  

y <- LULC$BBWA
plot(x, y, main="clumpy_1 vs BBWA", 
     xlab="clumpy_1", ylab="BBWA count", 
     pch=19, col="red")

y <- LULC$TEWA
plot(x, y, main="clumpy_1 vs TEWA", 
     xlab="clumpy_1", ylab="TEWA count", 
     pch=19, col="forestgreen")  



### Age_mn_1 vs response
x <- LULC$age_mn_1
y <- LULC$BTNW
plot(x, y, main="age_mn_1 vs BTNW", 
     xlab="age_mn_1", ylab="BTNW count", 
     pch=19, col="blue")

y <- LULC$BBWA
plot(x, y, main="age_mn_1 vs BBWA", 
     xlab="age_mn_1", ylab="BBWA count", 
     pch=19, col="red")  

y <- LULC$TEWA
plot(x, y, main="age_mn_1 vs TEWA", 
     xlab="age_mn_1", ylab="TEWA count", 
     pch=19, col="forestgreen")  


model <- lm(NTEMS$prop_con_1 ~ NTEMS$prop_dec_1)
summary(model)














### Create scatter plot of prop conifer and broadleaf deciduous


# Convert the shapefile into a data frame
AB_point_counts_df <- AB_point_counts_sf %>% 
  st_drop_geometry()  # Remove the spatial geometry to keep only attribute data

# Scatterplot of forest proportions
plot(
  AB_point_counts_df$prop.conifer, 
  AB_point_counts_df$prop.broadleaf,
  xlab = "Proportion of Coniferous Trees",
  ylab = "Proportion of Broadleaf (Deciduous) Trees",
  main = "Forest Type Distribution",
  xlim = c(0, 1),
  ylim = c(0, 1),
  pch = 16,  # Filled circles
  col = "forestgreen"
)

# Adding grid lines for better interpretation
abline(h = 0.5, col = "grey", lty = 2)  # Horizontal line at 0.5
abline(v = 0.5, col = "grey", lty = 2)  # Vertical line at 0.5
grid()




################## Create bar chart showing distribution of types

# Bin the proportions into 0.1 increments
AB_point_counts_df <- AB_point_counts_df %>%
  mutate(
    conifer_bins = cut(
      prop.conifer, 
      breaks = seq(0, 1, by = 0.1), 
      include.lowest = TRUE, 
      labels = seq(0, 0.9, by = 0.1)
    ),
    broadleaf_bins = cut(
      prop.broadleaf, 
      breaks = seq(0, 1, by = 0.1), 
      include.lowest = TRUE, 
      labels = seq(0, 0.9, by = 0.1)
    )
  )

# Summarize the counts for each bin
conifer_counts <- AB_point_counts_df %>%
  count(conifer_bins) %>%
  rename(bin = conifer_bins, count = n)

broadleaf_counts <- AB_point_counts_df %>%
  count(broadleaf_bins) %>%
  rename(bin = broadleaf_bins, count = n)

# Plot bar chart for proportion of coniferous trees
barplot(
  conifer_counts$count,
  names.arg = conifer_counts$bin,
  xlab = "Proportion of Coniferous Trees",
  ylab = "Number of Sites",
  main = "Distribution of Proportion of Coniferous Trees",
  col = "forestgreen",
  border = "white"
)

# Plot bar chart for proportion of broadleaf (deciduous) trees
barplot(
  broadleaf_counts$count,
  names.arg = broadleaf_counts$bin,
  xlab = "Proportion of Broadleaf Trees",
  ylab = "Number of Sites",
  main = "Distribution of Proportion of Broadleaf Trees",
  col = "goldenrod",
  border = "white"
)
