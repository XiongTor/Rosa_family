# !/usr/bin/Rscript
# Author: Tao Xiong
# Date: 2025.03.28
# This script is used to do the name standardization by WCVP

# *******  The input file shoule include the species name column named "Species"  *********

file<-commandArgs(TRUE)

# library packages
library(rWCVP)
library(rWCVPdata)
library(dplyr)
library(tidyverse)

rWCVP <- rWCVPdata::wcvp_names
rWCVP$plant_name_id <- as.double(rWCVP$plant_name_id)


sp <- read.csv(file[1], header = TRUE)
colnames(sp) <- "Species"

# get the accepted names from wcvp
globalWoodiness_matches <- wcvp_match_names(
  sp,
  name_col = "Species",
  fuzzy = FALSE,
  progress_bar = FALSE
)


# get the fuzzy name
fuzzy_matches <- globalWoodiness_matches %>%
  filter(str_detect(match_type, "Fuzzy")) %>%
  mutate(
    keep = case_when( #set up a keep column
      match_similarity < 0.9 ~ NA_real_, # fill with blank for dissimilar names
      # wcvp_author_edit_distance == 0 ~ 1, # fill with 1 if authors identical
      match_edit_distance == 1 ~ 1, # fill with 1 if only one letter different
    )
  )


if(nrow(fuzzy_matches)>0){
  write.csv(fuzzy_matches, "fuzzy_matches.csv",quote = F,row.names = FALSE)
}else {
  message("No fuzzy match data")
}


#Match with the wcvpdata
wcvp_matches <- globalWoodiness_matches %>%
  left_join(rWCVP, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  select(Species,wcvp_name,taxon_name,wcvp_status,taxon_status) %>% 
  distinct()

# get the accepted name
accepted_matches <- wcvp_matches %>%
  filter(taxon_status=="Accepted" & 
         wcvp_status != "Illegitimate" & 
         wcvp_status != "Misapplied" & 
         wcvp_status !="Invalid" & 
         wcvp_status !="NA" & 
         wcvp_status !="Unplaced") %>% 
  distinct()


accepted_matches2 <- accepted_matches %>% 
  group_by(Species) %>% 
  # 对每个组执行操作（无论是否重复）
  group_modify(~ {
    if (nrow(.x) == 1) { 
      # 无重复：直接保留 
      .x 
    } else {
      # 有重复：检查 Accepted 行
      accepted_rows <- .x %>% filter(wcvp_status == "Accepted")
      if (nrow(accepted_rows) >= 1) {
        # 存在 Accepted 行，保留所有 Accepted 行
        accepted_rows
      } else {
        synonym_codes <- c("Synonym", "Orthographic", "Artificial Hybrid", "Unplaced")
        synonyms <- .x %>%
          filter(wcvp_status %in% synonym_codes)
        if (nrow(synonyms) == 1)  {
          synonyms
        }else{
          .x %>% slice(0)
        }
        # 无 Accepted 行，保留第一行
      }
    }
  }) %>% 
  ungroup()

# get the no match name

df <- setdiff(wcvp_matches$Species,accepted_matches2$Species)

no_match <- wcvp_matches %>% 
  filter(Species %in% df)

write.csv(accepted_matches2,"Final_accepted_name.csv",quote = F,row.names = F)
write.csv(no_match,"No_matched_name.csv",quote = F,row.names = F)

