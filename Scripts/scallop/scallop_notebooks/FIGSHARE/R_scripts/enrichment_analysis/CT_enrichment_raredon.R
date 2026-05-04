setwd("/media/amon/2TB-HD/Olga/scallop-notebooks/scallop-notebooks/scallop_biological/scallop_directory/R_enrichment_analysis/")
library(emmeans)
library(tidyverse)
path2dir = '/media/amon/2TB-HD/Olga/scallop-notebooks/scallop-notebooks/scallop_biological/scallop_directory/output/tables/cell_type_enrichment/'

dataset_name = 'Raredon'
filepath = paste0(path2dir, 'input/', dataset_name, '_input_for_enrichment.csv')
read_csv(filepath) -> cell_type_counts

cell_type_counts %>% group_by(Condition, Donor, Cell_type)
cell_type_counts %>% filter(Total > 0) -> df_a
model0 <- glm(
  formula = cbind(Count, Other) ~ Cell_type,
  family = binomial(link = 'logit'),
  data = df_a
)
emm0 <- emmeans(model0, specs = ~ Cell_type)
emm0 %>%
  summary(infer = TRUE, type = 'response') %>%
  arrange(prob) -> cell_type_probs
cell_type_probs %>% head

cell_type_probs %>% pull(prob) %>% sum

df_a %>% pull(Count) %>% sum -> N
df_a %>% 
  group_by(Cell_type) %>% 
  summarise(fraction = sum(Count) / N) %>% 
  arrange(fraction) %>% 
  as.data.frame -> cell_type_fractions

cell_type_fractions %>% left_join(cell_type_probs) -> cell_type_probs

df_a %>% filter(Condition %in% c('Young', 'Old')) -> df

formula = cbind(Count, Other) ~ Cell_type * Condition

model1 <- glm(formula = formula, family = 'binomial', data = df)

emm1 <- emmeans(model1, specs = revpairwise ~ Condition | Cell_type)
emm1$contrasts %>%
  summary(infer = TRUE, type = 'response') %>%
  rbind() %>%
  as.data.frame() -> c_results

c_results %>% arrange(desc(odds.ratio))

emm2 <- emmeans(model1, specs = ~ Cell_type)
emm2 %>%
  summary(type = 'response') %>%
  select(Cell_type, prob) -> mean_probs

c_results %>% left_join(mean_probs) -> m_results

write.table(m_results, paste0(path2dir, 'output/', dataset_name, '_enrichment_output.csv'), sep=',')


