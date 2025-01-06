output_data_folder_path<-paste0("/Users/capelastegui.f/git/bovine_tree_figure/Outputs/",tdate,"/")
if (!dir.exists(output_data_folder_path)) {
  dir.create(output_data_folder_path, recursive = TRUE)
}


alignment <- readAAStringSet(paste0("/Users/capelastegui.f/git/bovine_tree_figure/Data/04_trees/",tdate,"/",tdate,"_raw_concat_aligned_trimmed_clean.fasta"), format = "fasta")

# Get mutations easily:
reference_name <- "EPI_ISL_19014396_A/Wild_Bird/Wyoming/24_003692_001/2024_2024_01_25"  
reference_name_short <- "EPI_ISL_19014396_A/Wild_Bird/Wyoming"  

texas37_name <- "EPI_ISL_19027114_A/Texas/37/2024_2024_03_28"
texas37_name_short <- "EPI_ISL_19027114_A/Texas/37"

# Extract the reference sequence by its name
reference_sequence <- alignment[[reference_name]]

# Initialize an empty dataframe to store all mismatches
mismatch_df <- data.frame(ref = character(),
                          query = character(),
                          mutation = character(),
                          position = integer(),
                          stringsAsFactors = FALSE)

for (i in seq_along(alignment)) {
  
  # Extract query sequence
  query_sequence <- alignment[[i]]
  query_name <- names(alignment)[i]
  
  # Perform pairwise alignment
  alignment_result <- pairwiseAlignment(reference_sequence, query_sequence)
  
  # Extract mismatches
  mismatches <- mismatchTable(alignment_result)
  
  # Check if there are any mismatches to process
  if (nrow(mismatches) > 0) {
    # Create a mutation column based on mismatches
    mismatches <- mismatches %>%
      mutate(mutation = paste0(PatternSubstring, SubjectStart, SubjectSubstring))  # Adjust mutation format as needed
    
    # Add to the main mismatch dataframe
    temp_df <- data.frame(
      ref = reference_name,
      query = query_name,
      mutation = mismatches$mutation,
      position = mismatches$SubjectStart,
      stringsAsFactors = FALSE
    )
    mismatch_df <- bind_rows(mismatch_df, temp_df)
  } else {
    # If no mismatches, add a row indicating no mutation
    mismatch_df <- bind_rows(mismatch_df, data.frame(
      ref = reference_name,
      query = query_name,
      mutation = "No Mutation",
      position = NA,
      stringsAsFactors = FALSE
    ))
  }
}

mutation_data <-  mismatch_df %>% select(-ref) 

#Re number the mutations so that they align to the postion on each segment, not the number of the superalignment

mutation_data$PA <- ifelse(mutation_data$position <= 716, mutation_data$mutation, NA)
mutation_data$PB2 <- ifelse(mutation_data$position > 716 & mutation_data$position < 1472, mutation_data$mutation, NA)
mutation_data$PB1 <- ifelse(mutation_data$position > 1472 & mutation_data$position < 2233, mutation_data$mutation, NA)
mutation_data$NP <- ifelse(mutation_data$position > 2233 & mutation_data$position < 2731, mutation_data$mutation, NA)
mutation_data$'NA' <- ifelse(mutation_data$position > 2731 & mutation_data$position < 3200, mutation_data$mutation, NA)
mutation_data$HA <- ifelse(mutation_data$position > 3200 & mutation_data$position < 3766, mutation_data$mutation, NA)

#Check the below
mutation_data$NS1 <- ifelse(mutation_data$position > 3766 & mutation_data$position < 3997, mutation_data$mutation, NA)
mutation_data$NEP <- ifelse(mutation_data$position > 3997 & mutation_data$position < 4118, mutation_data$mutation, NA)
mutation_data$M1 <- ifelse(mutation_data$position > 4118 & mutation_data$position < 4370, mutation_data$mutation, NA)
mutation_data$M2 <- ifelse(mutation_data$position > 4370 & mutation_data$position < 4467, mutation_data$mutation, NA)

mutation_data$position_new <- ifelse(
  mutation_data$position <= 716, mutation_data$position,
  ifelse(mutation_data$position > 716 & mutation_data$position < 1472, mutation_data$position - 716,
  ifelse(mutation_data$position > 1472 & mutation_data$position < 2233, mutation_data$position - 1475,
  ifelse(mutation_data$position > 2233 & mutation_data$position < 2731, mutation_data$position - 2232,
  ifelse(mutation_data$position > 2731 & mutation_data$position < 3200, mutation_data$position - 2730,
  ifelse(mutation_data$position > 3200 & mutation_data$position < 3766, mutation_data$position - 3199,
  ifelse(mutation_data$position > 3766 & mutation_data$position < 3997, mutation_data$position - 3765,
  ifelse(mutation_data$position > 3997 & mutation_data$position < 4118, mutation_data$position - 3996,
  ifelse(mutation_data$position > 4118 & mutation_data$position < 4370, mutation_data$position - 4117,
  ifelse(mutation_data$position > 4370 & mutation_data$position < 4467, mutation_data$position - 4369,
    NA # Assign NA if none of the conditions match
     ))))))))))

# Function to replace digits in the column string with position_new
replace_digits_with_position <- function(column_name, position_new_col) {
  mutation_data[[column_name]] <- mapply(function(mutation_str, new_number) {
    gsub("\\d+", as.character(new_number), mutation_str)  # Replace digits with new_number
  }, mutation_data[[column_name]], mutation_data[[position_new_col]])
}

# List of columns to apply the function to
columns_to_process <- c("PB1", "PB2", "PA", "NP", "NA", "HA", "NS1", "NEP", "M1", "M2")

# Apply the mapply function to each column to replace digits with the corresponding position_new value
for (col in columns_to_process) {
  mutation_data[[col]] <- mapply(function(mutation_str, new_number) {
    gsub("\\d+", as.character(new_number), mutation_str)  # Replace digits with new_number
  }, mutation_data[[col]], mutation_data$position_new)
}


mutation_data<- mutation_data %>%
  mutate(mutation_final = coalesce(PA, PB1, PB2, NP, NA, HA, NS1, NEP, M1, M2))%>%
  mutate(
    segment = case_when(
      !is.na(PA) ~ "PA",
      !is.na(PB1) ~ "PB1",
      !is.na(PB2) ~ "PB2",
      !is.na(NP) ~ "NP",
      !is.na(NA) ~ "NA",
      !is.na(HA) ~ "HA",
      !is.na(NS1) ~ "NS1",
      !is.na(NEP) ~ "NEP",
      !is.na(M1) ~ "M1",
      !is.na(M2) ~ "M2",
      
      TRUE ~ NA_character_
    )
  )

write_csv(mutation_data, paste0(output_data_folder_path,"/mutation_data.csv"))



mutation_wide <- mutation_data %>% 
  select(query, mutation_final, segment) %>% 
  mutate(segment_mutation = paste0(segment, "_", mutation_final)) %>%  # Combine 'segment' and 'mutation_final'
  select(query, segment_mutation) %>%  # Use the new 'segment_mutation' column
  pivot_wider(id_cols = query, 
              values_from = segment_mutation,
              names_from = segment_mutation) %>% 
  select(-NA_NA)

#######
# Polymerase mutations of interest:

mutation_wide_final <- mutation_wide %>% 
  select(query,
         PA_E613K,
         PA_I13V,
         PA_K497R,
         PB2_E362G,
         PB2_M631L,
         PB2_D740N,
         PB2_Q591R,
         PB2_E677G) %>% 
  rename(tip.label = query) %>% 
  mutate(across(everything(), ~ gsub("^$", NA, .)))%>%
  mutate(across(everything(), ~ gsub("NULL", NA, .))) %>% 
  mutate(concat = apply(select(., -tip.label), 1, function(x) paste(na.omit(x), collapse = ","))) %>% 
mutate(concat = ifelse(tip.label== "EPI_ISL_19014396_A/Wild_Bird/Wyoming/24_003692_001/2024_2024_01_25", "Ref", concat))

polymerase_geno <- mutation_wide_final %>% 
  select(tip.label, concat)

ref_data <- polymerase_geno %>% filter(concat=="Ref")
polymerase_heat_map <- mutation_wide_final %>%
  mutate(across(c(-tip.label, -concat), ~ ifelse(!is.na(.) & . != "", "YES", "NO"))) %>% 
  mutate(group = ifelse(grepl("Ref", concat), 0, 1)) %>%  
  mutate(branch_colour = case_when(
    grepl("EPI_ISL_19014396", tip.label) ~ 1,
    grepl("EPI_ISL_19027114_A/Texas/37", tip.label) ~ 2,
    TRUE ~ 0
  ))

polymerase_heat_map_tbl <-data.frame(polymerase_heat_map) 
rownames(polymerase_heat_map_tbl) <- polymerase_heat_map$tip.label
polymerase_heat_map_tbl<- polymerase_heat_map_tbl %>% select(-tip.label, -concat, -group)


#######
ha_mutations <- mutation_data  %>%
filter(segment=="HA") %>% 
  select(query, mutation_final, segment) %>% 
  pivot_wider(
    id_cols = query, 
    names_from = mutation_final, 
    values_from = mutation_final
  ) %>% 
  rename(tip.label = query) 

ha_heatmap <-data.frame(ha_mutations) 
rownames(ha_heatmap) <- ha_mutations$tip.label
ha_heatmap<- ha_heatmap %>% select(-tip.label)%>%
  mutate(across(everything(), ~ na_if(., " "))) %>% 
  mutate(across(everything(),~ifelse(!is.na(.) & . != "", "YES", "NO"))) 



write_csv(ha_mutations, paste0(output_data_folder_path,"/ha_mutations.csv"))
# %>% 
#   mutate(group = ifelse(grepl("Wyoming", concat), 0, 1))

#Save out GISAID IDS:

gisaid_names <- alignment@ranges@NAMES %>% as.data.frame() %>% 
  mutate(gisaid_name = str_extract(., "EPI_ISL_\\d+")) %>% 
  select(gisaid_name)

write_csv(gisaid_names, paste0(output_data_folder_path,"/gisaid_names.csv"))
