setwd("E:/Project/Rosacea/Rosaceae_fruit_type/character_data")
library(rWCVP)
library(rWCVPdata)

rosa_ga <- read.csv("rosa_genus_subfam_inf.csv")
#==================================================
###TRY
# install.packages('rtry')
# 加载rtry包
library(rtry)
library(dplyr)

# 1.1 使用rtry_import函数进行数据导入
mydata <- rtry_import("39259_14022025080608/39259.txt")

# 1.2 使用rtry_explore函数进行数据探索
mydata_explore <- rtry_explore(
  mydata,
  DataID,
  DataName,
  TraitID,
  TraitName,
  sortBy = TraitID,
  showOverview = FALSE
)
# View(mydata_explore)

# 1.4 现在数据集的信息可能有点冗余，要精简信息可以使用下面的方法
# 使用rtry_remove_col删除某一列
workdata <- rtry_remove_col(mydata, V29)
# 使用rtry_select_col函数选择需要的列
workdata <- rtry_select_col(
  workdata,
  ObsDataID,
  ObservationID,
  AccSpeciesID,
  AccSpeciesName,
  ValueKindName,
  TraitID,
  TraitName,
  DataID,
  DataName,
  OriglName,
  OrigValueStr,
  OrigUnitStr,
  StdValue,
  UnitName,
  OrigObsDataID,
  ErrorRisk,
  Comment
)
# 使用rtry_select_row函数选择所有性状记录和感兴趣的辅助数据
workdata <- rtry_select_row(workdata, TraitID > 0)

# 1.5 使用rtry_exclude函数排除数据，根据实际研究选择

# 根据错误风险排除离群值 （ErrorRisk），这里排除ErrorRisk大于等于3的值
workdata <- rtry_exclude(workdata, ErrorRisk >= 3, baseOn = ObsDataID)

# 1.6 使用rtry_remove_dup函数根据重复标识符 （OrigObsDataID） 删除重复项
workdata <- rtry_remove_dup(workdata)

# 1.7 使用rtry_trans_wider函数转换为宽表（长表格式是出于数据管理目的，宽表格式才是我们平时处理数据比较习惯的格式）
# 选择包含完整TraitID和StdValue的行
num_traits <- rtry_select_row(workdata, complete.cases(TraitID))
# 从筛选后的数据中选择指定的列
# num_traits <- rtry_select_col(num_traits, ObservationID, AccSpeciesID, AccSpeciesName, TraitID, TraitName, StdValue, UnitName)

#去除整列都为NA的列
num_traits <- num_traits %>%
  select(where(~ any(!is.na(.))))

#提取出属名
num_traits$AccGenusname <- sapply(
  strsplit(num_traits$AccSpeciesName, " "),
  '[',
  1
)

#提取出所属蔷薇科的属

try_rosa <- num_traits %>%
  filter(AccGenusname %in% rosa_ga$Genera)

#名称标准化与替换
#名称标准化
data1 <- try_rosa %>%
  select(AccSpeciesName) %>%
  distinct()

globalWoodiness_matches <- wcvp_match_names(
  data1,
  name_col = "AccSpeciesName",
  fuzzy = FALSE,
  progress_bar = FALSE
)

manually_resolved <- filter(
  globalWoodiness_matches,
  wcvp_status != "Illegitimate" &
    wcvp_status != "Misapplied" &
    wcvp_status != "Invalid" &
    wcvp_status != "NA"
)


rWCVP <- rWCVPdata::wcvp_names
rWCVP$plant_name_id <- as.double(rWCVP$plant_name_id)

accepted_matches <- manually_resolved %>%
  left_join(rWCVP, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  mutate(
    keep = case_when(
      taxon_status == "Accepted" & (wcvp_status != "Synonym" | wcvp_homotypic) ~
        "Matched to an accepted name",
      TRUE ~ "Not matched to an accepted name"
    )
  )


globalWoodiness_accepted <- accepted_matches %>%
  filter(accepted_plant_name_id %in% accepted_matches$wcvp_accepted_id) %>%
  group_by(genus) %>%
  filter(
    sum(wcvp_status == "Unplaced") == n() |
      wcvp_status != "Unplaced"
  ) %>%
  ungroup() %>%
  select(AccSpeciesName, taxon_name, wcvp_status) %>%
  distinct()


#与原始数据匹配
# 步骤 1：将两个数据框按 AccSpeciesName 列连接
try_rosa <- try_rosa %>%
  filter(AccSpeciesName %in% globalWoodiness_accepted$AccSpeciesName) %>%
  left_join(
    globalWoodiness_accepted %>% select(AccSpeciesName, taxon_name),
    by = "AccSpeciesName"
  )

# 步骤 2：替换 AccSpeciesName 列为 taxon_name（匹配时）
try_rosa <- try_rosa %>%
  mutate(
    AccSpeciesName = ifelse(
      is.na(taxon_name),
      AccSpeciesName, # 未匹配时保留原值
      taxon_name # 匹配时替换为 taxon_name
    )
  ) %>%
  select(-taxon_name) # 删除临时列


#再次去除哪些整列都为NA的列
final_ga_data <- try_rosa %>%
  select(where(~ any(!is.na(.))))

#提取出唯一的TraitID对应的形状
character_ID <- final_ga_data %>% select(TraitID, TraitName) %>% distinct()

# 将数据转换为宽表格式
num_traits_georef_wider <- rtry_trans_wider(
  final_ga_data,
  names_from = c(TraitID),
  values_from = c(OrigValueStr),
  values_fn = list(StdValue = mean)
)
#更换表头，筛选出所需列
num_traits_georef_wider <- num_traits_georef_wider %>%
  rename_with(~ coalesce(character_ID[[2]][match(.x, character_ID[[1]])], .x))

num_traits_georef_wider <- num_traits_georef_wider %>%
  select(
    AccSpeciesName,
    TraitName,
    DataName,
    OriglName,
    OrigUnitStr,
    StdValue,
    Comment,
    "Dispersal syndrome",
    "Flower color",
    "Fruit type",
    "Inflorescence type",
    "Flower petal number",
    "Flower number of pollen per ovule",
    "Flower morphology type",
    "Fruit length"
  )

#仅提取属
num_traits_georef_wider$genus <- sapply(
  strsplit(num_traits_georef_wider$AccSpeciesName, " "),
  '[',
  1
)

#去重
num_traits_genus <- num_traits_georef_wider %>%
  select(
    genus,
    TraitName,
    DataName,
    OriglName,
    OrigUnitStr,
    StdValue,
    Comment,
    "Dispersal syndrome",
    "Flower color",
    "Fruit type",
    "Inflorescence type",
    "Flower petal number",
    "Flower number of pollen per ovule",
    "Flower morphology type",
    "Fruit length"
  ) %>%
  distinct() %>%
  as.data.frame()

#将所有列均转化为字符串格式
num_traits_genus <- as.data.frame(lapply(num_traits_genus, as.character))

colnames(num_traits_genus) <- gsub(" ", "_", colnames(num_traits_genus))

# 1.8 使用rtry_export函数导出预处理后的 TRY 数据
rtry_export(num_traits_genus, "TRY_rosa_traits.csv")


##==============================================================================
#BIEN
# install.packages('BIEN')
library(RPostgreSQL)
library(DBI)
library(BIEN)
library(dplyr)
library(tidyr)

# 如果我们对某一性状感兴趣，第一步是检查该性状是否存在，并使用函数 BIEN_trait_list 验证拼写
BIEN_trait_list()

# 获取特定科的所有性状数据
family_traits <- BIEN_trait_family(family = "Rosaceae")

#读出整个蔷薇科的数据库
# write.csv(family_traits,"BIEN_rosaceae_trait_data.csv",row.names = F)
family_traits <- read.csv("BIEN_rosaceae_all_trait_data.csv")

head(family_traits)
# 获取多个科的所有性状数据，以下类似
# family_traits <- BIEN_trait_family(family = c("Poaceae","Orchidaceae"))
# View(family_traits)

# # 获取特定属的性状数据
# genus_traits <- BIEN_trait_genus(genus = "Acer")
# head(genus_traits)
#
# # 获取特定物种的性状数据
# species_traits <- BIEN_trait_species(species = "Poa annua")
# head(species_traits)
#
# # 获取特定性状的所有记录
# leaf_area_traits <- BIEN_trait_trait(trait = "leaf area")
# head(leaf_area_traits)
# # 获取多个性状的所有记录
# BIEN_traits <- BIEN_trait_trait(trait = c("whole plant height", "leaf dry mass per leaf fresh mass"))
# head(BIEN_traits)

# 估算物种平均性状值
# BIEN_trait_mean函数在没有物种水平数据的情况下，使用属或科水平数据估算给定性状的物种平均值
# BIEN_trait_mean(species=c("Poa annua","Juncus trifidus"),trait="leaf dry mass per leaf fresh mass")

# # 获取特定科和性状的数据
# family_leaf_area <- BIEN_trait_traitbyfamily(trait = "whole plant height", family = "Poaceae")
# head(family_leaf_area)

# # 获取特定属和性状的数据
# genus_wood_density <- BIEN_trait_traitbygenus(trait = "whole plant height", genus = "Carex")
# head(genus_wood_density)

# # 获取特定物种和性状的数据
# species_hl <- BIEN_trait_traitbyspecies(trait = c("whole plant height", "leaf area"), species = c("Carex capitata","Betula nana"))
# head(species_hl)

# BIEN_trait_species 提取物种国家的性状数据
# BIEN_trait_country("South Africa")
# BIEN_trait_country(country="South Africa",trait="whole plant growth form")

#挑选出需要的性状
family_traits2 <- family_traits %>%
  filter(
    trait_name %in%
      c(
        "flower color",
        "flower pollination syndrome",
        "maximum fruit length",
        "maximum whole plant height",
        "seed length",
        "seed mass",
        "whole plant dispersal syndrome",
        "whole plant height",
        "whole plant sexual system",
        "whole plant woodiness"
      )
  )

#提取出method和unit对应的性状解释
family_traits_um <- family_traits2 %>%
  select(trait_name, method, unit) %>%
  distinct()

write.csv(family_traits2, "BIEN_method.csv")

#挑选需要的列
family_traits3 <- family_traits2 %>%
  select(
    scrubbed_family,
    scrubbed_genus,
    scrubbed_species_binomial,
    trait_name,
    trait_value
  ) %>%
  filter(!is.na(trait_value))


#分别提取出value这一列为数字和字符的行
data_chara <- family_traits3 %>%
  mutate(is_numeric = !is.na(as.numeric(trait_value))) %>% # 判断是否能转换为数字
  filter(!is_numeric) %>% # 筛选出不能转换的行
  select(-is_numeric)

data_num <- family_traits3 %>%
  mutate(is_numeric = !is.na(as.numeric(trait_value))) %>% # 判断是否能转换为数字
  filter(is_numeric) %>% # 筛选出不能转换的行
  select(-is_numeric)

data_num$trait_value <- as.numeric(data_num$trait_value)

data_num <- data_num %>%
  filter(!is.na(data_num$scrubbed_species_binomial)) %>%
  distinct()


#长表格转换为宽表格
family_traits_wide_num <- data_num %>%
  pivot_wider(
    id_cols = c(
      "scrubbed_family",
      "scrubbed_genus",
      "scrubbed_species_binomial"
    ),
    names_from = trait_name,
    values_from = trait_value,
    values_fn = ~ median(.) # 对重复值取中位值
  )


family_traits_wide_character <- data_chara %>%
  pivot_wider(
    id_cols = c(
      "scrubbed_family",
      "scrubbed_genus",
      "scrubbed_species_binomial"
    ),
    names_from = trait_name,
    values_from = trait_value # 对重复值取中位值
  )

#挑选需要保留的行并转换成数据框格式
final_data_bien <- full_join(
  family_traits_wide_num,
  family_traits_wide_character,
  by = "scrubbed_species_binomial"
) %>%
  select(-c(scrubbed_family.y, scrubbed_genus.y)) %>%
  distinct() %>%
  as.data.frame()

#将所有的列都转化为字符串格式，防止部分列为list格式
final_data_bien_df <- data.frame(lapply(final_data_bien, as.character))

#######

#名称标准化
globalWoodiness_matches <- wcvp_match_names(
  final_data_bien_df,
  name_col = "scrubbed_species_binomial",
  fuzzy = FALSE,
  progress_bar = FALSE
)

manually_resolved <- filter(
  globalWoodiness_matches,
  wcvp_status != "Illegitimate" &
    wcvp_status != "Misapplied" &
    wcvp_status != "Invalid" &
    wcvp_status != "NA"
)


rWCVP <- rWCVPdata::wcvp_names
rWCVP$plant_name_id <- as.double(rWCVP$plant_name_id)

accepted_matches <- manually_resolved %>%
  left_join(rWCVP, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  mutate(
    keep = case_when(
      taxon_status == "Accepted" & (wcvp_status != "Synonym" | wcvp_homotypic) ~
        "Matched to an accepted name",
      TRUE ~ "Not matched to an accepted name"
    )
  )


globalWoodiness_accepted <- accepted_matches %>%
  filter(accepted_plant_name_id %in% accepted_matches$wcvp_accepted_id) %>%
  group_by(genus) %>%
  filter(
    sum(wcvp_status == "Unplaced") == n() |
      wcvp_status != "Unplaced"
  ) %>%
  ungroup() %>%
  select(scrubbed_species_binomial, taxon_name, wcvp_status, ) %>%
  distinct()


#与原始数据匹配
# 步骤 1：将两个数据框按 scrubbed_species_binomial 列连接
final_data_bien_df <- final_data_bien_df %>%
  filter(
    scrubbed_species_binomial %in%
      globalWoodiness_accepted$scrubbed_species_binomial
  ) %>%
  left_join(
    globalWoodiness_accepted %>% select(scrubbed_species_binomial, taxon_name),
    by = "scrubbed_species_binomial"
  )

# 步骤 2：替换 AccSpeciesName 列为 taxon_name（匹配时）
final_data_bien_df <- final_data_bien_df %>%
  mutate(
    scrubbed_species_binomial = ifelse(
      is.na(taxon_name),
      scrubbed_species_binomial, # 未匹配时保留原值
      taxon_name # 匹配时替换为 taxon_name
    )
  ) %>%
  select(-taxon_name) # 删除临时列


#再次去除哪些整列都为NA的列
final_data_bien_df <- final_data_bien_df %>%
  select(where(~ any(!is.na(.))))

#更新替换属名
final_data_bien_df$scrubbed_genus.x <- sapply(
  strsplit(final_data_bien_df$scrubbed_species_binomial, " "),
  '[',
  1
)


#读出数据
utils::write.csv(final_data_bien_df, "BIEN_rosaceae_character_data.csv")

#=================================================================================
#GIFT
# 安装并加载GIFT包
# install.packages("GIFT")
library("GIFT")
# 设置超时:查询需要很长时间才能完成，增加超时功能可以更容易地完成较大的下载
options(timeout = max(10000000, getOption("timeout")))

# 查看性状具体的ID，包括每个性状的类型和内容
trait_meta <- GIFT_traits_meta()
write.csv(trait_meta, "GIFT_trait_meta.csv")

# 使用GIFT_traits_tax函数在更高的分类学水平上检索性状
trait_tax <- GIFT_traits(trait_IDs = , bias_ref = FALSE, bias_deriv = FALSE)


#确认可以成功下载后开始使用循环下载数据

trait_tax_list <- list()

trait_name <- c(
  "1.2.2",
  "1.6.3",
  "2.4.1",
  "3.10.3",
  "3.11.3",
  "3.16.1",
  "3.21.1",
  "3.22.1",
  "3.23.1",
  "3.3.1",
  "3.3.2",
  "3.9.1",
  "5.1.1"
)

for (trait in trait_name) {
  trait_tax <- GIFT_traits(
    trait_IDs = trait,
    bias_ref = FALSE,
    bias_deriv = FALSE
  )
  trait_tax_list[[trait]] <- trait_tax
}

GIFT_data <- bind_rows(trait_tax_list)
GIFT_data$genus <- sapply(strsplit(GIFT_data$work_species, " "), '[', 1)

#挑选出仅蔷薇科有的属
GIFT_rosa_data <- merge(
  GIFT_data,
  rosa_ga,
  by.x = "genus",
  by.y = "Genera",
  all.y = T
)

# GIFT_rosa_data <- read.csv("GIFT_rosa_character_data.csv")

#名称标准化
globalWoodiness_matches <- wcvp_match_names(
  GIFT_rosa_data,
  name_col = "work_species",
  fuzzy = FALSE,
  progress_bar = FALSE
)

manually_resolved <- filter(
  globalWoodiness_matches,
  wcvp_status != "Illegitimate" &
    wcvp_status != "Misapplied" &
    wcvp_status != "Invalid" &
    wcvp_status != "NA"
)


rWCVP <- rWCVPdata::wcvp_names
rWCVP$plant_name_id <- as.double(rWCVP$plant_name_id)

accepted_matches <- manually_resolved %>%
  left_join(rWCVP, by = c("wcvp_accepted_id" = "plant_name_id")) %>%
  mutate(
    keep = case_when(
      taxon_status == "Accepted" & (wcvp_status != "Synonym" | wcvp_homotypic) ~
        "Matched to an accepted name",
      TRUE ~ "Not matched to an accepted name"
    )
  )


globalWoodiness_accepted <- accepted_matches %>%
  filter(accepted_plant_name_id %in% accepted_matches$wcvp_accepted_id) %>%
  group_by(genus.x) %>%
  filter(
    sum(wcvp_status == "Unplaced") == n() |
      wcvp_status != "Unplaced"
  ) %>%
  ungroup() %>%
  select(work_species, taxon_name, wcvp_status) %>%
  distinct()


#与原始数据匹配
# 步骤 1：将两个数据框按 work_species 列连接
GIFT_rosa_data <- GIFT_rosa_data %>%
  filter(work_species %in% globalWoodiness_accepted$work_species) %>%
  left_join(
    globalWoodiness_accepted %>% select(work_species, taxon_name),
    by = "work_species"
  )

# 步骤 2：替换 work_species 列为 taxon_name（匹配时）
GIFT_rosa_data <- GIFT_rosa_data %>%
  mutate(
    work_species = ifelse(
      is.na(taxon_name),
      work_species, # 未匹配时保留原值
      taxon_name # 匹配时替换为 taxon_name
    )
  ) %>%
  select(-taxon_name) # 删除临时列

#更新替换属名
GIFT_rosa_data$genus <- sapply(
  strsplit(GIFT_rosa_data$work_species, " "),
  '[',
  1
)


#再次去除哪些整列都为NA的列
GIFT_rosa_data <- GIFT_rosa_data %>%
  select(where(~ any(!is.na(.))))


#读出数据
write.csv(GIFT_rosa_data, "GIFT_rosa_character_data.csv")


#===============================================================================

#+==============================================================================
#合并TRY，BIEN，GIFT和之前自测收集的属级，种级别数据集合

#加载需要的分析包
library(dplyr)

#读取全部的性状数据

TRY <- read.csv("TRY_rosa_traits.csv")

#归并为一个属的类型
##用属名进行分组后，对每一组，选择每一列的唯一值，并使用“,”分隔开
TRY_genus <- TRY %>%
  group_by(genus) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%
  select(-c(X, TraitName, DataName, OriglName, OrigUnitStr, StdValue, Comment))

#给列名加上TRY的后缀
TRY_name <- data.frame(
  old <- c(
    "genus",
    "Family",
    "Subfamilies",
    "Tribes",
    "Dispersal_syndrome",
    "Flower_color",
    "Fruit_type",
    "Inflorescence_type",
    "Flower_petal_number",
    "Flower_number_of_pollen_per_ovule",
    "Flower_morphology_type",
    "Fruit_length"
  ),
  new <- c(
    "genus(TRY)",
    "Family(TRY)",
    "Subfamilies(TRY)",
    "Tribes(TRY)",
    "Dispersal_syndrome(TRY)",
    "Flower_color(TRY)",
    "Fruit_type(TRY)",
    "Inflorescence_type(TRY)",
    "Flower_petal_number(TRY)",
    "Flower_number_of_pollen_per_ovule(TRY)",
    "Flower_morphology_type(TRY)",
    "Fruit_length(TRY)"
  )
)

TRY_genus <- TRY_genus %>%
  rename_with(~ coalesce(TRY_name[[2]][match(.x, TRY_name[[1]])], .x))


#======================================

BIEN <- read.csv("BIEN_rosaceae_character_data.csv")
#归并为一个属
BIEN_genus <- BIEN %>%
  group_by(scrubbed_genus.x) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%
  select(-c(X))
#给列名加上BIEN的后缀
BIEN_name <- data.frame(
  old <- c(
    "scrubbed_genus.x",
    "scrubbed_family.x",
    "scrubbed_species_binomial",
    "whole.plant.height",
    "seed.mass",
    "seed.length",
    "maximum.whole.plant.height",
    "maximum.fruit.length",
    "whole.plant.dispersal.syndrome",
    "whole.plant.sexual.system",
    "flower.pollination.syndrome",
    "whole.plant.woodiness",
    "flower.color",
    "Subfamilies",
    "Tribes"
  ),
  new <- c(
    "genus(BIEN)",
    "family(BIEN)",
    "species(BIEN)",
    "whole_plant_height(BIEN)",
    "seed_mass(BIEN)",
    "seed_length(BIEN)",
    "maximum_whole_plant_height(BIEN)",
    "maximum_fruit_length(BIEN)",
    "whole_plant_dispersal_syndrome(BIEN)",
    "whole_plant_sexual_system(BIEN)",
    "flower_pollination_syndrome(BIEN)",
    "whole_plant_woodiness(BIEN)",
    "flower_color(BIEN)",
    "Subfamilies(BIEN)",
    "Tribes(BIEN)"
  )
)

BIEN_genus <- BIEN_genus %>%
  rename_with(~ coalesce(BIEN_name[[2]][match(.x, BIEN_name[[1]])], .x))

#===========================

GIFT <- read.csv("GIFT_rosa_character_data.csv")
#归并为一个属
GIFT_genus <- GIFT %>%
  group_by(genus) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", "))) %>%
  select(-c(X, X.1, X.2, work_ID, work_species, work_author)) %>%
  select(
    -c(
      "agreement_1.2.2",
      "n_1.2.2",
      "references_1.2.2",
      "cv_1.6.3",
      "n_1.6.3",
      "references_1.6.3",
      "agreement_2.4.1",
      "n_2.4.1",
      "references_2.4.1",
      "n_3.10.3",
      "references_3.10.3",
      "n_3.11.3",
      "references_3.11.3",
      "agreement_3.16.1",
      "n_3.16.1",
      "references_3.16.1",
      "agreement_3.21.1",
      "n_3.21.1",
      "references_3.21.1",
      "agreement_3.22.1",
      "n_3.22.1",
      "references_3.22.1",
      "agreement_3.23.1",
      "n_3.23.1",
      "references_3.23.1",
      "agreement_3.3.1",
      "n_3.3.1",
      "references_3.3.1",
      "agreement_3.3.2",
      "n_3.3.2",
      "references_3.3.2",
      "agreement_3.9.1",
      "n_3.9.1",
      "references_3.9.1"
    )
  )

#给列名加上GIFT的后缀
GIFT_name <- data.frame(
  old <- c(
    "genus",
    "trait_value_1.2.2",
    "trait_value_1.6.3",
    "trait_value_2.4.1",
    "trait_value_3.10.3",
    "trait_value_3.11.3",
    "trait_value_3.16.1",
    "trait_value_3.21.1",
    "trait_value_3.22.1",
    "trait_value_3.23.1",
    "trait_value_3.3.1",
    "trait_value_3.3.2",
    "trait_value_3.9.1",
    "Family",
    "Subfamilies",
    "Tribes"
  ),
  new <- c(
    "genus(GIFT)",
    "Growth_form(GIFT)",
    "Plant_height_mean(GIFT)",
    "Deciduousness(GIFT)",
    "Seed_length_mean(GIFT)",
    "Seed_width_mean(GIFT)",
    "Fruit_type(GIFT)",
    "Flower_colour (GIFT)",
    "Fruit colour(GIFT)",
    "Inflorescence(GIFT)",
    "Dispersal_syndrome_1(GIFT)",
    "Dispersal_syndrome_2(GIFT)",
    "Seeds_per_fruit(GIFT)",
    "Family(GIFT)",
    "Subfamilies(GIFT)",
    "Tribes(GIFT)"
  )
)

GIFT_genus <- GIFT_genus %>%
  rename_with(~ coalesce(GIFT_name[[2]][match(.x, GIFT_name[[1]])], .x))


#===================================
#读取自己手动收集的数据
rosa_genus <- read.csv("Rosaceae_ga_character_data.csv") %>%
  select(-"chinese")

#给列名加上rosa_genus的后缀
rosa_genus_name <- data.frame(
  old <- c(
    "family",
    "subfamily",
    "genus",
    "fruit_type",
    "inflorescence",
    "carpel",
    "Ovary.location",
    "life.form",
    "petal",
    "calyx"
  ),
  new <- c(
    "family(rosa_genus)",
    "subfamily(rosa_genus)",
    "genus(rosa_genus)",
    "fruit_type(rosa_genus)",
    "inflorescence(rosa_genus)",
    "carpel(rosa_genus)",
    "Ovary.location(rosa_genus)",
    "life.form(rosa_genus)",
    "petal(rosa_genus)",
    "calyx(rosa_genus)"
  )
)

rosa_genus <- rosa_genus %>%
  rename_with(
    ~ coalesce(rosa_genus_name[[2]][match(.x, rosa_genus_name[[1]])], .x)
  )


#===================================
rosa_sp <- read.csv("rosa_sp_character_data.csv")
rosa_sp$genus <- sapply(strsplit(rosa_sp$Species_name, " "), '[', 1)
#归并为一个属
rosa_sp_to_genus <- rosa_sp %>%
  group_by(genus) %>%
  summarise(across(everything(), ~ paste(unique(.), collapse = ", ")))

#给列名加上rosa_sp_to_genus的后缀
rosa_sp_name <- data.frame(
  old <- c(
    "genus",
    "Species_name",
    "Anthocyanidin_coloring_type",
    "Petal_color",
    "Petal_number",
    "Prick",
    "Inflorescence"
  ),
  new <- c(
    "genus(rosa_sp)",
    "Species_name(rosa_sp)",
    "Anthocyanidin_coloring_type(rosa_sp)",
    "Petal_color(rosa_sp)",
    "Petal_number(rosa_sp)",
    "Prick(rosa_sp)",
    "Inflorescence(rosa_sp)"
  )
)

rosa_sp_to_genus <- rosa_sp_to_genus %>%
  rename_with(~ coalesce(rosa_sp_name[[2]][match(.x, rosa_sp_name[[1]])], .x))


#=============================
#读取蔷薇科的属级和种级接受名单
rosa_accepted_genus <- read.csv("rosa_genus_subfam_inf.csv")
rosa_accepted_sp <- read.csv("rosa_species_subfam_inf.csv")

#========================
#读取染色体相关信息
rosa_chromosome <- read.csv("Rosaceae_genus_chromosome_num.csv")

#==========================
#开始进行数据合并
final <- rosa_accepted_genus %>%
  full_join(TRY_genus, by = c("Genera" = "genus(TRY)")) %>%
  full_join(BIEN_genus, by = c("Genera" = "genus(BIEN)")) %>%
  full_join(GIFT_genus, by = c("Genera" = "genus(GIFT)")) %>%
  full_join(rosa_genus, by = c("Genera" = "genus(rosa_genus)")) %>%
  full_join(rosa_sp_to_genus, by = c("Genera" = "genus(rosa_sp)")) %>%
  full_join(rosa_chromosome, by = c("Genera" = "genus"))

final <- final %>%
  select(
    -c("family(rosa_genus)", "subfamily(rosa_genus)", "Species_name(rosa_sp)")
  )

#给列名重排序
final <- final[, c(
  "Family",
  "Subfamilies",
  "Tribes",
  "Genera",
  "Anthocyanidin_coloring_type(rosa_sp)",
  "calyx(rosa_genus)",
  "carpel(rosa_genus)",
  "Deciduousness(GIFT)",
  "Dispersal_syndrome(TRY)",
  "Dispersal_syndrome_1(GIFT)",
  "Dispersal_syndrome_2(GIFT)",
  "flower_color(BIEN)",
  "Flower_color(TRY)",
  "Flower_colour (GIFT)",
  "Flower_morphology_type(TRY)",
  "Flower_number_of_pollen_per_ovule(TRY)",
  "Flower_petal_number(TRY)",
  "flower_pollination_syndrome(BIEN)",
  "Fruit colour(GIFT)",
  "Fruit_length(TRY)",
  "Fruit_type(GIFT)",
  "fruit_type(rosa_genus)",
  "Fruit_type(TRY)",
  "Growth_form(GIFT)",
  "Inflorescence(GIFT)",
  "inflorescence(rosa_genus)",
  "Inflorescence(rosa_sp)",
  "Inflorescence_type(TRY)",
  "life.form(rosa_genus)",
  "maximum_fruit_length(BIEN)",
  "maximum_whole_plant_height(BIEN)",
  "Ovary.location(rosa_genus)",
  "petal(rosa_genus)",
  "Petal_color(rosa_sp)",
  "Petal_number(rosa_sp)",
  "Plant_height_mean(GIFT)",
  "Prick(rosa_sp)",
  "seed_length(BIEN)",
  "Seed_length_mean(GIFT)",
  "seed_mass(BIEN)",
  "Seed_width_mean(GIFT)",
  "Seeds_per_fruit(GIFT)",
  "whole_plant_dispersal_syndrome(BIEN)",
  "whole_plant_height(BIEN)",
  "whole_plant_sexual_system(BIEN)",
  "whole_plant_woodiness(BIEN)",
  "gametophytic",
  "sporophytic",
  "parsed_n"
)]

write.csv(final, "Final_rosaceae_genus_character_data2.csv", row.names = F)

#更新处理一下数据
data <- read.csv("Final_rosaceae_genus_character_data2.csv")
data <- data[c(1:109), ]
data <- data %>%
  mutate(across(where(is.factor), as.character)) %>% # 先转字符
  mutate(across(where(is.character), ~ gsub("/", ";", ., fixed = TRUE)))
