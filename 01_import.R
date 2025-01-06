#Set directory where data is saved
raw_data_file <-"/Users/capelastegui.f/git/bovine_tree_figure/Data/01_raw_data/20250103"

#Read in raw cattle sequences
raw_cattle <- readAAStringSet(paste0(raw_data_file,"/","20241227_raw_dairy_cattle_pb1_pb2_na_ha_pa_np_ns1_ns2_m1_m2_nep.fasta"))

#Read in raw avian comparator sequences
raw_avian <- readAAStringSet(paste0(raw_data_file,"/","20241227_goose_wyoming_raw_pb1_pb2_pa_np_na_ha_ns1_ns2_m1_m2_nep.fasta"))

#Read in raw texas37  sequences
raw_texas37 <- readAAStringSet(paste0(raw_data_file,"/","20250102_texas37_pb1_pb2_na_ha_pa_np_ns1_ns2_m1_m2_nep.fasta"))

#raw proe 2024 

raw_animal_2023 <-readAAStringSet(paste0(raw_data_file,"/","20250103_avian_outgroup.fasta"))

#Separate the segments by name in the sequence name stem
raw_cattle_pb2 <- raw_cattle[grep("_PB2", names(raw_cattle))]
raw_cattle_pa <- raw_cattle[grep("_PA", names(raw_cattle))]
raw_cattle_pb1 <- raw_cattle[grep("_PB1", names(raw_cattle))]
raw_cattle_np <- raw_cattle[grep("_NP", names(raw_cattle))]
raw_cattle_na <- raw_cattle[grep("_NA", names(raw_cattle))]
raw_cattle_ha <- raw_cattle[grep("_HA", names(raw_cattle))]
raw_cattle_ns1 <- raw_cattle[grep("_NS1", names(raw_cattle))]
raw_cattle_nep <- raw_cattle[grep("_NEP", names(raw_cattle))]
raw_cattle_m1 <- raw_cattle[grep("_M1", names(raw_cattle))]
raw_cattle_m2 <- raw_cattle[grep("_M2", names(raw_cattle))]

#Separate the segments by name in the sequence name stem
raw_avian_pb2 <- raw_avian[grep("_PB2", names(raw_avian))]
raw_avian_pa <- raw_avian[grep("_PA", names(raw_avian))]
raw_avian_pb1 <- raw_avian[grep("_PB1", names(raw_avian))]
raw_avian_np <- raw_avian[grep("_NP", names(raw_avian))]
raw_avian_na <- raw_avian[grep("_NA", names(raw_avian))]
raw_avian_ha <- raw_avian[grep("_HA", names(raw_avian))]
raw_avian_ns1 <- raw_avian[grep("_NS1", names(raw_avian))]
raw_avian_nep <- raw_avian[grep("_NEP", names(raw_avian))]
raw_avian_m1 <- raw_avian[grep("_M1", names(raw_avian))]
raw_avian_m2 <- raw_avian[grep("_M2", names(raw_avian))]


raw_texas37_pb2 <- raw_texas37[grep("_PB2", names(raw_texas37))]
raw_texas37_pa <- raw_texas37[grep("_PA", names(raw_texas37))]
raw_texas37_pb1 <- raw_texas37[grep("_PB1", names(raw_texas37))]
raw_texas37_np <- raw_texas37[grep("_NP", names(raw_texas37))]
raw_texas37_na <- raw_texas37[grep("_NA", names(raw_texas37))]
raw_texas37_ha <- raw_texas37[grep("_HA", names(raw_texas37))]
raw_texas37_ns1 <- raw_texas37[grep("_NS1", names(raw_texas37))]
raw_texas37_nep <- raw_texas37[grep("_NEP", names(raw_texas37))]
raw_texas37_m1 <- raw_texas37[grep("_M1", names(raw_texas37))]
raw_texas37_m2 <- raw_texas37[grep("_M2", names(raw_texas37))]


raw_animal_2023_pb2 <- raw_animal_2023[grep("_PB2", names(raw_animal_2023))]
raw_animal_2023_pa <- raw_animal_2023[grep("_PA", names(raw_animal_2023))]
raw_animal_2023_pb1 <- raw_animal_2023[grep("_PB1", names(raw_animal_2023))]
raw_animal_2023_np <- raw_animal_2023[grep("_NP", names(raw_animal_2023))]
raw_animal_2023_na <- raw_animal_2023[grep("_NA", names(raw_animal_2023))]
raw_animal_2023_ha <- raw_animal_2023[grep("_HA", names(raw_animal_2023))]
raw_animal_2023_ns1 <- raw_animal_2023[grep("_NS1", names(raw_animal_2023))]
raw_animal_2023_nep <- raw_animal_2023[grep("_NEP", names(raw_animal_2023))]
raw_animal_2023_m1 <- raw_animal_2023[grep("_M1", names(raw_animal_2023))]
raw_animal_2023_m2 <- raw_animal_2023[grep("_M2", names(raw_animal_2023))]




#remove the segment identifier i. "_PA" so that the names can be joined

name_function <- function(filtered_sequences) {
  # Remove "_PA" and "_PB2" from the names
  names(filtered_sequences) <- gsub("_segment_PA", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_PB2", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_PB1", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_NP", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_NA", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_HA", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_NS1", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_NS2", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_NEP", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_M1", "", names(filtered_sequences))
  names(filtered_sequences) <- gsub("_segment_M2", "", names(filtered_sequences))
  
  names(filtered_sequences) <- gsub("-", "_", names(filtered_sequences))
  
  # Return the modified sequences
  return(filtered_sequences)
}

#Do this for avian and cattle:

raw_avian_pa<- name_function(raw_avian_pa)
raw_avian_pb2<- name_function(raw_avian_pb2)
raw_avian_pb1<- name_function(raw_avian_pb1)
raw_avian_np<- name_function(raw_avian_np)
raw_avian_na<- name_function(raw_avian_na)
raw_avian_ha<- name_function(raw_avian_ha)
raw_avian_ns1<- name_function(raw_avian_ns1)
raw_avian_nep<- name_function(raw_avian_nep)
raw_avian_m1<- name_function(raw_avian_m1)
raw_avian_m2<- name_function(raw_avian_m2)

raw_cattle_pa<- name_function(raw_cattle_pa)
raw_cattle_pb2<- name_function(raw_cattle_pb2)
raw_cattle_pb1<- name_function(raw_cattle_pb1)
raw_cattle_np<- name_function(raw_cattle_np)
raw_cattle_na<- name_function(raw_cattle_na)
raw_cattle_ha<- name_function(raw_cattle_ha)
raw_cattle_ns1<- name_function(raw_cattle_ns1)
raw_cattle_nep<- name_function(raw_cattle_nep)
raw_cattle_m1<- name_function(raw_cattle_m1)
raw_cattle_m2<- name_function(raw_cattle_m2)

raw_texas37_pa<- name_function(raw_texas37_pa)
raw_texas37_pb2<- name_function(raw_texas37_pb2)
raw_texas37_pb1<- name_function(raw_texas37_pb1)
raw_texas37_np<- name_function(raw_texas37_np)
raw_texas37_na<- name_function(raw_texas37_na)
raw_texas37_ha<- name_function(raw_texas37_ha)
raw_texas37_ns1<- name_function(raw_texas37_ns1)
raw_texas37_nep<- name_function(raw_texas37_nep)
raw_texas37_m1<- name_function(raw_texas37_m1)
raw_texas37_m2<- name_function(raw_texas37_m2)

raw_animal_2023_pa<- name_function(raw_animal_2023_pa)
raw_animal_2023_pb2<- name_function(raw_animal_2023_pb2)
raw_animal_2023_pb1<- name_function(raw_animal_2023_pb1)
raw_animal_2023_np<- name_function(raw_animal_2023_np)
raw_animal_2023_na<- name_function(raw_animal_2023_na)
raw_animal_2023_ha<- name_function(raw_animal_2023_ha)
raw_animal_2023_ns1<- name_function(raw_animal_2023_ns1)
raw_animal_2023_nep<- name_function(raw_animal_2023_nep)
raw_animal_2023_m1<- name_function(raw_animal_2023_m1)
raw_animal_2023_m2<- name_function(raw_animal_2023_m2)

#stack the sequences of the same segment (avian and cattle)
all_pa <- c(raw_animal_2023_pa,raw_cattle_pa, raw_texas37_pa) 
all_pb2 <- c(raw_animal_2023_pb2, raw_cattle_pb2, raw_texas37_pb2)
all_pb1 <- c(raw_animal_2023_pb1, raw_cattle_pb1, raw_texas37_pb1)
all_np <- c(raw_animal_2023_np, raw_cattle_np, raw_texas37_np)
all_na <- c(raw_animal_2023_na, raw_cattle_na, raw_texas37_na)
all_ha <- c(raw_animal_2023_ha, raw_cattle_ha, raw_texas37_ha)
all_ns1 <- c(raw_animal_2023_ns1, raw_cattle_ns1, raw_texas37_ns1)
all_nep <- c(raw_animal_2023_nep, raw_cattle_nep, raw_texas37_nep)
all_m1 <- c(raw_animal_2023_m1, raw_cattle_m1, raw_texas37_m1)
all_m2 <- c(raw_animal_2023_m2, raw_cattle_m2,raw_texas37_m2)

#Find sequence names that appear in all the sequences; this means whole genomes can be stuck together

common_names <- Reduce(intersect, list(
  names(all_pa), 
  names(all_pb2), 
  names(all_pb1),
  names(all_np),
  names(all_na),
  names(all_ha),
  names(all_ns1),
  names(all_nep),
  names(all_m1),
  names(all_m2)
))
concatenated_sequences <- lapply(common_names, function(name) {
  # Concatenate sequences by name
  c(all_pa[[name]], all_pb2[[name]],all_pb1[[name]],all_np[[name]],all_na[[name]],all_ha[[name]],all_ns1[[name]],all_nep[[name]],all_m1[[name]],all_m2[[name]])
})

# Step 4: Set names for the concatenated sequences
names(concatenated_sequences) <- common_names

concatenated_sequences<-AAStringSet(concatenated_sequences)


writeXStringSet(concatenated_sequences, filepath=paste0(raw_data_file,"/",tdate,"_raw_concat.fasta"))


