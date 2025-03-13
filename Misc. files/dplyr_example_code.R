# Load necessary library
library(dplyr)

# Create a bird count dataset
set.seed(123)  # For reproducibility
bird_data <- data.frame(
  Location = sample(c("Wetland A", "Forest B", "Grassland C", "Urban D"), 30, replace = TRUE),
  Date = sample(seq(as.Date("2023-05-01"), as.Date("2023-05-30"), by = "day"), 30, replace = TRUE),
  Bird_Species = sample(c("Robin", "Sparrow", "Blue Jay", "Woodpecker", "Cardinal"), 30, replace = TRUE),
  Observer = sample(c("Observer 1", "Observer 2", "Observer 3"), 30, replace = TRUE),
  Count = sample(0:50, 30, replace = TRUE),
  Temperature = round(runif(30, 10, 25), 1),  # Temperature in °C
  Wind_Speed = round(runif(30, 0, 15), 1)     # Wind speed in km/h
)

library(dplyr)

filtered <- bird_data %>%
  select(Location, Date, Bird_Species)

grasslandC <- bird_data %>%
  filter(Location == "Grassland C")

temp_new_column <- bird_data %>%
  mutate(temp_ordinal = case_when(
    Temperature < 5 ~ "Low", 
    Temperature <= 20 ~ "Moderate", 
    Temperature > 20 ~ "High" 
    ))

# Find the average Count of birds observed by each Observer on days with Temperature > 15°C
filtered_grouped_summary <- bird_data %>%
  filter(Temperature > 15) %>%
  group_by(Observer) %>%
  summarize(Average_Count = mean(Count), .groups = "drop")
print(filtered_grouped_summary)











# Log reg example

# Load necessary library
library(dplyr)

# Step 1: Simulate the Dataset with Statistically Significant Data
set.seed(123)  # For reproducibility
bird_data <- data.frame(
  Temperature = runif(100, 5, 30),  # Temperature in °C
  Habitat_Type = sample(c("Forest", "Wetland", "Grassland"), 100, replace = TRUE)
)

# Create a relationship between Temperature and Bird_Detected
# Higher temperature increases the likelihood of detection
bird_data$Bird_Detected <- rbinom(100, 1, plogis(-5 + 0.3 * bird_data$Temperature))

# Step 2: Fit a Logistic Regression Model
# Predicting bird detection based on temperature
model <- glm(Bird_Detected ~ Temperature, data = bird_data, family = "binomial")

# Summarize the model
summary(model)

# Step 3: Make Predictions
# Create a new data frame for prediction
new_data <- data.frame(Temperature = c(10, 15, 20, 25, 30))

# Predict probabilities
new_data$Predicted_Probability <- predict(model, newdata = new_data, type = "response")

# View the predictions
print(new_data)

# Step 4: Visualize the Logistic Regression Curve
# Plot the observed data
plot(bird_data$Temperature, bird_data$Bird_Detected, 
     pch = 19, col = "blue", xlab = "Temperature (°C)", ylab = "Bird Detected (1 = Yes, 0 = No)")

# Add the logistic regression curve
curve(predict(model, newdata = data.frame(Temperature = x), type = "response"), 
      col = "red", lwd = 2, add = TRUE)

