# !/usr/bin/Rscript
# Author: Tao Xiong
# This script is used to classify the dispersal type of each species based on the dispersal terms provided in the dataset.
# And to explore the relashionship between dispersal type and distribution range.


library(dplyr)
library(stringr)

raw_data <- read.csv("Dispersal.csv")

# 定义传播方式到距离的映射规则（覆盖所有已知术语）
distance_mapping <- list(
  "极短距离（<1米）" = c(
    "Seeds Drop To The Ground Close To Or Beneath The Parent Plant",
    "Autochory", "Dry Seed", "Boleochory", "Germinule", "Generative Dispersule",
    "Semachor", "Hemerochor", "Blastochor", "Agochor", "Dummy", "Unassisted"
  ),
  "短距离（1-10米）" = c(
    "Adhesion", "Epizoochory", "Ant", "Snail", "Dormouse", "Hedgehog", 
    "Rainwash", "Omburochor", "Chamaechor", "Bythisochor", "Erosion Material",
    "Clothes And Footwear", "Seed Contamination", "Speirochor", "Pig", "Rabbit"
  ),
  "中距离（10米-1公里）" = c(
    "Mammals (non-bat)", "Ethelochor", "Meteorochor",
    "Mammals", "Fox", "Hare", "Marten", "Squirrel", "Mouse", "Sheep", "Cattle",
    "Deer", "Goat", "Wild Boar", "Badger", "Horse", "Pacarana", "Chamois",
    "Marmot", "Vertebrate", "Fleshy", "Shaken Fresh Water", "Standing Fresh Water",
    "Water Nautochor", "Nautochor", "Snail", 
    "Turdus Merula (Thrush)", "Sylvia Atricapilla (Blackcap)", "Vulpes Vulpes (Red Fox)",
    "Beech Marten", "Pine Marten", "Roe", "Diaspore Is Eaten Intentionally",
    "Seed Ingested Or Regurgitated","Zoochory"
  ),
  "长距离（>1公里）" = c(
    "Bird", "Birds","Hydrochor","Wind", "Anemochory","Endozoochory", "Commerce", "Water", "Anthropochory", "Hydrochorous",
    "Multi", "Air", "Dispersal Endozoochorous", "Dysochor", "Man", "Domestic Animal")
)

# 清洗函数：统一格式并标准化术语
clean_terms <- function(text) {
  text %>%
    # 替换所有分隔符为逗号
    str_replace_all("[;.,+（）()]", ",") %>%
    # 去除特殊符号和多余空格
    str_replace_all("\\*|\\?|\\s+", " ") %>%
    str_trim() %>%
    # 分割为独立术语
    str_split("\\s*,\\s*") %>%
    unlist() %>%
    # 标准化术语
    sapply(function(x) {
      case_when(
        str_detect(x, regex("zoochor(y|ous)", ignore_case = TRUE)) ~ "Zoochory",
        str_detect(x, regex("anemochor(y|ous)", ignore_case = TRUE)) ~ "Anemochory",
        str_detect(x, regex("epizoochor(y|ous)", ignore_case = TRUE)) ~ "Epizoochory",
        str_detect(x, regex("endozoochor(y|ous)", ignore_case = TRUE)) ~ "Endozoochory",
        str_detect(x, regex("hydrochor(y|ous)", ignore_case = TRUE)) ~ "Hydrochorous",
        str_detect(x, regex("dry seed", ignore_case = TRUE)) ~ "Dry Seed",
        str_detect(x, regex("air;?dry", ignore_case = TRUE)) ~ "Dry Seed",
        TRUE ~ str_to_title(x)  # 其余术语首字母大写
      )
    }) %>%
    unique() %>%  # 去重
    paste(collapse = ", ")  # 合并为字符串
}

# 分类函数：根据清洗后的术语确定最高优先级距离
classify_distance <- function(terms) {
  terms <- unlist(str_split(terms, ",\\s*"))
  
  # 按优先级检查每个术语
  if (any(terms %in% distance_mapping[["长距离（>1公里）"]])) {
    "长距离"
  } else if (any(terms %in% distance_mapping[["中距离（10米-1公里）"]])) {
    "中距离"
  } else if (any(terms %in% distance_mapping[["短距离（1-10米）"]])) {
    "短距离"
  } else if (any(terms %in% distance_mapping[["极短距离（<1米）"]])) {
    "极短距离"
  } else {
    "未分类"
  }
}

##############################属级水平#######################################

# 假设原始数据存储在数据框 raw_data 中，列名为 Dispersal
final_result <- raw_data %>%
  mutate(
    Cleaned_Terms = sapply(Dispersal, clean_terms),
    Distance_Category = sapply(Cleaned_Terms, classify_distance)
  ) %>%
  select(Genera,Dispersal, Cleaned_Terms, Distance_Category)



unclassified <- final_result %>%
  filter(Distance_Category == "未分类")

if (nrow(unclassified) > 0) {
  message("以下条目需检查分类规则：")
  print(unclassified$Cleaned_Terms)
} else {
  message("所有条目已成功分类！")
}

num <- read.csv("Final_genus_diff_num.csv")

final <- merge(num,final_result,by.x = "genus",by.y = "Genera",all.x = T)

write.csv(final,"Earth_num/Earth_num_dispersal.csv",fileEncoding = "GB18030",row.names = F)

#############################种级水平#######################################
data <- read.csv("TRY_rosa_sp_traits.csv")

data <- data %>% 
  select(AccSpeciesName,Dispersal.syndrome) %>% 
  distinct() %>% 
  group_by(AccSpeciesName) %>%
  summarise(across(everything(), ~paste(unique(.), collapse = ",")))

final_result <- data %>%
  mutate(
    Cleaned_Terms = sapply(Dispersal.syndrome, clean_terms),
    Distance_Category = sapply(Cleaned_Terms, classify_distance)
  ) %>%
  select(AccSpeciesName,Dispersal.syndrome, Cleaned_Terms, Distance_Category)



unclassified <- final_result %>%
  filter(Distance_Category == "未分类")

if (nrow(unclassified) > 0) {
  message("以下条目需检查分类规则：")
  print(unclassified$Cleaned_Terms)
} else {
  message("所有条目已成功分类！")
}
#=======================================================

num <- read.csv("Earth_num/Final_sp_diff_num.csv")

final <- merge(num,final_result,by.x = "species",by.y = "AccSpeciesName",all.x = T) %>% 
  select(species,diff_count,common_continents,diff_continents,Distance_Category) %>% 
  filter(Distance_Category!="NA")

final <- final %>% 
  filter(Distance_Category!="未分类")

final$Distance_Category <- factor(
  final$Distance_Category,
  levels = c("极短距离", "短距离", "中距离", "长距离"),
  ordered = TRUE
)
###########################与分布范围是否具有一定的相关性#######################
num2 <- read.csv("Earth_num/Final_sp_native_introduced.csv")

result <- num2 %>%
  group_by(taxon_name) %>%
  summarise(across(everything(), ~paste(unique(.), collapse = ",")),
            n_rows = n())  # n_rows用于统计每组的行数

final <- merge(result,final_result,by.x = "taxon_name",by.y = "AccSpeciesName",all.x = T) %>% 
  select(taxon_name,n_rows,Distance_Category) %>% 
  filter(Distance_Category!="NA")

final <- final %>% 
  filter(Distance_Category!="未分类")

final$Distance_Category <- factor(
  final$Distance_Category,
  levels = c("极短距离", "短距离", "中距离", "长距离"),
  ordered = TRUE
)

#============================


#转化为数据类型：
final$distance_rank <- as.numeric(final$Distance_Category)


final <- final[order(final$n_rows), ]
library(ggplot2)

ggplot(final, aes(x = reorder(taxon_name, n_rows))) + 
  geom_point(aes(y = n_rows, group = 1, color = "distribution_range")) +  # 使用 linewidth 替代 size
  geom_point(aes(y = distance_rank, group = 1, color = "dispersal ability"),alpha=0.3) +
  labs(x = "Species", y = "Value", title = "Relationship between distribution range and dispersal ability") +
  scale_color_manual(values = c("distribution_range" = "blue", "dispersal ability" = "red")) +  # 自定义颜色
  theme_classic() +
  theme(axis.text.x = element_blank())  # 调整x轴标签角度
