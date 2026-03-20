library(rWCVP)
library(rWCVPdata)
library(dplyr)
# This script is used to get the whole Rosaceae species list.

#Get the APG IV plant classified information
apg <- taxonomic_mapping

#Get the rosaceae sp list
rosa = wcvp_checklist(taxon = "Rosaceae", taxon_rank = "family") %>% 
  select(family,genus,taxon_name,accepted_name,taxon_status) %>% 
  filter(accepted_name!="NA" & taxon_status=="Accepted") %>% 
  distinct()

#Get genus list
rosa_genus <- rosa %>% 
  select(family,genus) %>% 
  distinct()

# write.csv(rosa_genus,"latest_25.03.16/WCVP_accepted_rosaceae_genus_list.csv",row.names = F)

#Get species list
rosa_species <- rosa %>% 
  select(family,accepted_name) %>% 
  distinct()

######################### Compare different name #############################
powo <- read.csv("latest_25.03.16/POWO_accepted_rosaceae_genus_list_25.03.15.csv")

wcvp_only <- setdiff(rosa_genus$genus,powo$Genus)

powo_only <- setdiff(powo$Genus,rosa_genus$genus)

#After checked the unique name ,make a new file named "Rosaceae_genus_accepted.csv"


######################### Add the family message ###############################

#genus level
rosa <- read.csv("latest_25.03.16/Rosaceae_genus_accepted.csv")

subfam <- read.csv("rosa_genus_subfam_inf.csv")

# remove all empty line
subfam <- subfam %>%
  mutate(across(where(is.character),~na_if(.,""))) %>% 
  filter(!if_all(everything(),is.na))

#get subfamily information
final <- merge(subfam,rosa,by.x = "Genera",by.y = "genus",all.y=T)

# write.csv(final,"latest_25.03.16/Rosaceae_genus_accepted.csv",row.names = F)
#add two genus "Taihangia" and "Oncostylus",remove two genus "Aremonia","Brachycaulos"




#species level
old_sp <- read.csv("rosa_all_sp_240521.csv")
wcvp_only <- setdiff(rosa_species$accepted_name,old_sp$taxon_name) %>% 
  as.data.frame()

rosa_species$genus <- sapply(strsplit(rosa_species$accepted_name," "),'[',1)

final_sp <- merge(subfam,rosa_species,by.x="Genera",by.y = "genus",all.y = T) %>% 
  select(-family) %>% 
  filter(rosa_species$genus!="Brachycaulos",
         rosa_species$accepted_name!="Aremonia agrimonoides subsp. agrimonoides")

final_sp$Genera[final_sp$Genera=="Aremonia"] <- "Agrimonia"

final_sp$accepted_name[final_sp$accepted_name=="Aremonia agrimonoides subsp. pouzarii"] <- "Agrimonia agrimonoides subsp. pouzarii"

final_sp$accepted_name[final_sp$accepted_name=="Aremonia agrimonoides"] <- "Agrimonia agrimonoides"



#add two genus "Taihangia" and "Oncostylus",change "Aremonia" to the accepted name
write.csv(final_sp,"latest_25.03.16/Rosaceae_sp_accepted.csv",row.names = F)



########################## Get the more species ################################
#read the new data
new_data <- read.csv("latest_25.03.16/Rosaceae_sp_accepted.csv")
old_data <- read.csv("rosa_all_sp_240521.csv")

df <- setdiff(new_data$accepted_name,old_data$taxon_name) %>% 
  as.data.frame()
colnames(df) <- "species_name"

df2 <- setdiff(old_data$taxon_name,new_data$accepted_name)

df3 <- setdiff(new_data$Genera,old_data$genus)

write.csv(df,"latest_25.03.16/Rosaceae_new_more_old_sp.csv",row.names = F)










