library(rWCVP)
library(rWCVPdata)
names = rWCVPdata::wcvp_names
head(taxonomic_mapping)

orders <- readLines("Orders.txt")  


checklist <- do.call(rbind, lapply(orders, function(order) {
  df <- wcvp_checklist(order, taxon_rank = "order")
  df$order <- order
  return(df)
}))

asterids_list <- checklist %>% 
  select(order, family, genus, taxon_rank, taxon_status, 
         taxon_name, accepted_name, accepted_plant_name_id) %>% 
  unique()

asterids_accepted_list <- asterids_list %>% 
  filter(taxon_status == "Accepted" & taxon_rank == "Species") %>% 
  arrange(accepted_plant_name_id) %>% 
  unique()

write.csv(asterids_accepted_list, "./asterids_accept_list.csv", row.names = FALSE)