### FIND DIFFERENCES IN TAXA BETWEEN OVERLAPS BEFORE AND AFTER AND UNIQUES BEFORE AND AFTER ###

###All overlaps###
# Load the CSV files
file_before_all <- read.csv("taxonomy_results_with_info/Before_Left_Overlaps_All_taxonomy.csv", row.names = 1)
file_after_all <- read.csv("taxonomy_results_with_info/After_Left_Overlaps_All_taxonomy.csv", row.names = 1)


# Function to find unique and shared organisms
compare_organisms <- function(before, after) {
  # Get row names (organisms) from each dataset
  organisms_before <- rownames(before)
  organisms_after <- rownames(after)
  
  # Find unique and shared organisms
  unique_before <- setdiff(organisms_before, organisms_after)  # Only in "Before"
  unique_after <- setdiff(organisms_after, organisms_before)  # Only in "After"
  shared <- intersect(organisms_before, organisms_after)      # In both
  
  # Return results as a list
  list(
    unique_before = unique_before,
    unique_after = unique_after,
    shared = shared
  )
}

# Compare the organisms
results_all <- compare_organisms(file_before_all, file_after_all)


##############################################################


###Nonresponsive and Healthy Overlaps###
# Load the CSV files
file_before_nonresp_healthy <- read.csv("taxonomy_results_with_info/Before_Left_Overlaps_Nonresponsive_Healthy_taxonomy.csv", row.names = 1)
file_after_nonresp_healthy <- read.csv("taxonomy_results_with_info/After_Left_Overlaps_Nonresponsive_Healthy_taxonomy.csv", row.names = 1)


# Function to find unique and shared organisms
compare_organisms <- function(before, after) {
  # Get row names (organisms) from each dataset
  organisms_before <- rownames(before)
  organisms_after <- rownames(after)
  
  # Find unique and shared organisms
  unique_before <- setdiff(organisms_before, organisms_after)  # Only in "Before"
  unique_after <- setdiff(organisms_after, organisms_before)  # Only in "After"
  shared <- intersect(organisms_before, organisms_after)      # In both
  
  # Return results as a list
  list(
    unique_before = unique_before,
    unique_after = unique_after,
    shared = shared
  )
}

# Compare the organisms
results_nonresp_healthy <- compare_organisms(file_before_nonresp_healthy, file_after_nonresp_healthy)

###########################################################################

###Responsive and Healthy Overlaps###
# Load the CSV files
file_before_resp_healthy <- read.csv("taxonomy_results_with_info/Before_Left_Overlaps_Responsive_Healthy_taxonomy.csv", row.names = 1)
file_after_resp_healthy <- read.csv("taxonomy_results_with_info/After_Left_Overlaps_Responsive_Healthy_taxonomy.csv", row.names = 1)


# Function to find unique and shared organisms
compare_organisms <- function(before, after) {
  # Get row names (organisms) from each dataset
  organisms_before <- rownames(before)
  organisms_after <- rownames(after)
  
  # Find unique and shared organisms
  unique_before <- setdiff(organisms_before, organisms_after)  # Only in "Before"
  unique_after <- setdiff(organisms_after, organisms_before)  # Only in "After"
  shared <- intersect(organisms_before, organisms_after)      # In both
  
  # Return results as a list
  list(
    unique_before = unique_before,
    unique_after = unique_after,
    shared = shared
  )
}

# Compare the organisms
results_resp_healthy <- compare_organisms(file_before_resp_healthy, file_after_resp_healthy)


########################################################################

###Nonresponsive and Responsive Overlaps###
# Load the CSV files
file_before_resp_nonresp <- read.csv("taxonomy_results_with_info/Before_Left_Overlaps_Responsive_Nonresponsive_taxonomy.csv", row.names = 1)
file_after_resp_nonresp <- read.csv("taxonomy_results_with_info/After_Left_Overlaps_Responsive_Nonresponsive_taxonomy.csv", row.names = 1)


# Function to find unique and shared organisms
compare_organisms <- function(before, after) {
  # Get row names (organisms) from each dataset
  organisms_before <- rownames(before)
  organisms_after <- rownames(after)
  
  # Find unique and shared organisms
  unique_before <- setdiff(organisms_before, organisms_after)  # Only in "Before"
  unique_after <- setdiff(organisms_after, organisms_before)  # Only in "After"
  shared <- intersect(organisms_before, organisms_after)      # In both
  
  # Return results as a list
  list(
    unique_before = unique_before,
    unique_after = unique_after,
    shared = shared
  )
}

# Compare the organisms
results_resp_nonresp <- compare_organisms(file_before_resp_nonresp, file_after_resp_nonresp)

##############################################################################

###Nonresponsive Uniques###
# Load the CSV files
file_before_nonresp <- read.csv("taxonomy_results_with_info/Before_Left_Uniques_Unique_Nonresponsive_taxonomy.csv", row.names = 1)
file_after_nonresp <- read.csv("taxonomy_results_with_info/After_Left_Uniques_Unique_Nonresponsive_taxonomy.csv", row.names = 1)


# Function to find unique and shared organisms
compare_organisms <- function(before, after) {
  # Get row names (organisms) from each dataset
  organisms_before <- rownames(before)
  organisms_after <- rownames(after)
  
  # Find unique and shared organisms
  unique_before <- setdiff(organisms_before, organisms_after)  # Only in "Before"
  unique_after <- setdiff(organisms_after, organisms_before)  # Only in "After"
  shared <- intersect(organisms_before, organisms_after)      # In both
  
  # Return results as a list
  list(
    unique_before = unique_before,
    unique_after = unique_after,
    shared = shared
  )
}

# Compare the organisms
results_nonresp <- compare_organisms(file_before_nonresp, file_after_nonresp)

#############################################################################

###Responsive Uniques###
# Load the CSV files
file_before_resp <- read.csv("taxonomy_results_with_info/Before_Left_Uniques_Unique_Responsive_taxonomy.csv", row.names = 1)
file_after_resp <- read.csv("taxonomy_results_with_info/After_Left_Uniques_Unique_Responsive_taxonomy.csv", row.names = 1)


# Function to find unique and shared organisms
compare_organisms <- function(before, after) {
  # Get row names (organisms) from each dataset
  organisms_before <- rownames(before)
  organisms_after <- rownames(after)
  
  # Find unique and shared organisms
  unique_before <- setdiff(organisms_before, organisms_after)  # Only in "Before"
  unique_after <- setdiff(organisms_after, organisms_before)  # Only in "After"
  shared <- intersect(organisms_before, organisms_after)      # In both
  
  # Return results as a list
  list(
    unique_before = unique_before,
    unique_after = unique_after,
    shared = shared
  )
}

# Compare the organisms
results_resp <- compare_organisms(file_before_resp, file_after_resp)

