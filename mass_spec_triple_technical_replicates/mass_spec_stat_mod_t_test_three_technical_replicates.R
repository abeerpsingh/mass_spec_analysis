#Import libraries
library(tidyverse)
library(infer)
library(tidymodels)
library(ggrepel)
library(htmlwidgets)
#Functions-------------------------------------------------------------------------
ord_t_test <- function(data, protein = NA){
  data <- data %>%
    select(name, protein) %>%
    rownames_to_column() %>%
    pivot_wider(names_from = name, values_from = protein)
  
  std_error <- data %>%
    summarise(std_error = sqrt((var(test, na.rm = T) * (rep-1) + var(control, na.rm = T) * (rep-1)) /
                                 ((rep + rep -2) * ((1 / rep) + (1 / rep)))))
  
  diff_mean <- data %>%
    summarise(dx = mean(test, na.rm = T) - mean(control, na.rm = T))
  t <- diff_mean / std_error 
  p <- pt(q = -abs(t[["dx"]]), df = rep + rep -2)
  tibble(protein = protein, std_error = std_error[["std_error"]],
         diff_mean = diff_mean[["dx"]] , t = t[["dx"]], p = p)
}
#--------------------------------------------------------------------------------

#Import proteingroups data file
proteingroups_file_path <- file.choose()
protein_groups <- read_tsv(file = proteingroups_file_path) %>%
  rownames_to_column(var = "rownum") %>%
  mutate(rownum = str_c("p",rownum, sep = "_"),
         'Gene names' = str_replace(string = .$`Gene names`, pattern = ";", replacement = "|"))

#Remove false positive protein hits
protein_groups <- protein_groups %>%
  filter(is.na(`Only identified by site`) & is.na(Reverse) &
           is.na(`Potential contaminant`))

#Select the relevant samples
test <- readline(prompt = "Enter name of test protein:")
test_1 <- readline(prompt = "Enter name of 1st test sample:")
test_2 <- readline(prompt = "Enter name of 2nd test sample:")
test_3 <- readline(prompt = "Enter name of 3rd test sample:")
control <- readline(prompt = "Enter name of control sample: ")
control_1 <- readline(prompt = "Enter name of 1st control sample: ")
control_2 <- readline(prompt = "Enter name of 2nd control sample: ")
control_3 <- readline(prompt = "Enter name of 3rd control sample: ")
rep <- readline(prompt = "Number of replicates for test or control: ") %>%
  parse_integer()
data <- protein_groups %>%
  select("rownum", "Gene names", contains(test, ignore.case = T),
         contains(control, ignore.case = T))
#-----------------------------------------------------------------------------
data_sequence_coverage <- data %>%
  select(contains("Sequence coverage")) %>%
  rowwise() %>%
  transmute(sequence_coverage_test_1 = max(across(contains(test_1)), na.rm = T),
         sequence_coverage_test_2 = max(across(contains(test_2)), na.rm = T),
         sequence_coverage_test_3 = max(across(contains(test_3)), na.rm = T),
         sequence_coverage_control_1 = max(across(contains(control_1)), na.rm = T),
         sequence_coverage_control_2 = max(across(contains(control_2)), na.rm = T),
         sequence_coverage_control_3 = max(across(contains(control_3)), na.rm = T)) %>%
  ungroup()

data_unique_peptides <- data %>%
  select(starts_with("Unique peptides")) %>%
  rowwise() %>%
  transmute(unique_peptides_test_1 = max(across(contains(test_1)), na.rm = T),
            unique_peptides_test_2 = max(across(contains(test_2)), na.rm = T),
            unique_peptides_test_3 = max(across(contains(test_3)), na.rm = T),
            unique_peptides_control_1 = max(across(contains(control_1)), na.rm = T),
            unique_peptides_control_2 = max(across(contains(control_2)), na.rm = T),
            unique_peptides_control_3 = max(across(contains(control_3)), na.rm = T)) %>%
  ungroup()

 data_ibaq <- data %>%
  select(contains("iBAQ")) %>%
  na_if(y = 0) %>%
  transmute(ibaq_test_1 = pmap_dbl(select(., contains(test_1)),
                                .f = function(...) mean(c(...), na.rm = T)),
            ibaq_test_2 = pmap_dbl(select(., contains(test_2)),
                                   .f = function(...) mean(c(...), na.rm = T)),
            ibaq_test_3 = pmap_dbl(select(., contains(test_3)),
                                   .f = function(...) mean(c(...), na.rm = T)),
            ibaq_control_1 = pmap_dbl(select(., contains(control_1)),
                                   .f = function(...) mean(c(...), na.rm = T)),
            ibaq_control_2 = pmap_dbl(select(., contains(control_2)),
                                   .f = function(...) mean(c(...), na.rm = T)),
            ibaq_control_3 = pmap_dbl(select(., contains(control_3)),
                                   .f = function(...) mean(c(...), na.rm = T)))
 data_ibaq_names <- colnames(data_ibaq)
 data_ibaq <- colnames(data_ibaq) %>%
       map_dfc(.f = ~replace_na(data_ibaq[[.x]], 0))
 names(data_ibaq) <- data_ibaq_names  
 data <- bind_cols(data %>% select(rownum, "Gene names"), data_ibaq, data_sequence_coverage, data_unique_peptides)
# Select proteins: sequence coverage > 5 for atleast 2 out of 3 samples in atleast
# one experimental group
names_sequence_coverage <- data %>%
  select(rownum, starts_with("sequence_coverage")) %>%
  mutate(across(.cols = -rownum, .fns = ~between(x = .x, left = 0, right = 5))) %>%
  mutate(across(.cols = -rownum, .fns = as.integer)) %>%
  rowwise() %>%
  mutate(test = sum(across(.cols = contains("test"))),
         control = sum(across(.cols = contains("control")))) %>%
  ungroup() %>%
  filter(test < 2 | control < 2) %>%
  select(rownum) %>%
  drop_na()

# select proteins with unique peptide > 1 in two out of three samples in atleast 
# one experimental group
names_unique_peptides <- data %>%
  select(rownum, starts_with(match = "unique")) %>%
  mutate(across(.cols = -rownum, .fns = ~between(x = .x, left = 0,
                                                 right = 1)),
         across(.cols = -rownum, .fns = as.integer)) %>%
  rowwise() %>%
  mutate(test = sum(across(.cols = contains("test"))),
         control = sum(across(.cols = contains("control")))) %>%
  ungroup() %>%
  filter(test < 2 | control < 2) %>%
  select(rownum) %>%
  drop_na()

# Select proteins with valid values in 2 out of 3 replicates in both experimental
# groups
valid_value_names <- data %>%
  select(rownum, starts_with(match = "ibaq")) %>%
  mutate(across(.cols = -rownum, .fns = ~near(x = .x, y = 0)),
         across(.cols = -rownum, .fns = as.integer)) %>%
  rowwise() %>%
  mutate(test = sum(across(.cols = contains("test"))),
         control = sum(across(.cols = contains("control")))) %>%
  ungroup() %>%
  filter(test < 2 & control < 2) %>%
  select(rownum) %>%
  drop_na()

#Inner join the selected names: right join
names_for_analysis <- valid_value_names %>%
  inner_join(names_sequence_coverage, by = "rownum") %>%
  inner_join(names_unique_peptides, by = "rownum")

#filter out the selected names for analysis from the data: right join
data_combined <- right_join(x = data, y = names_for_analysis, by = "rownum")
data_iBAQ <- data_combined %>%
  select(rownum, starts_with("ibaq"), 'Gene names')

data_transposed <- data_combined %>%
  select(rownum, starts_with("ibaq")) %>%
  pivot_longer(-rownum) %>%
  pivot_wider(names_from = rownum) %>%
  mutate(name = if_else(str_detect(name, "test"), true = "test", 
                        false = "control")) %>%
  mutate(across(contains("p"), .fns = ~na_if(x = .x, y = 0)))

#Preprocess data: roll imputation(median)
data_transpose_impute <- data_transposed %>%
  recipe() %>%
  step_impute_roll(all_numeric(), window = 3) %>%
  prep(data_transposed) %>%
  bake(data_transposed)

# Calculate ordinary t-statistic and the p_value
final_result <- names_for_analysis$rownum %>%
  map_dfr(.f = ~ord_t_test(data = data_transpose_impute,
                           protein = .x))

#Calculate moderated t-statistic and p_value
final_result <- final_result %>%
  mutate(mod_std_error = (std_error + median(std_error, na.rm = T)) / 2,
         t_mod = diff_mean / mod_std_error,
         p_mod = pt(q = -abs(t_mod), df = rep + rep -2))

#Determine significance
final_data <- data_iBAQ %>%
  right_join(y = final_result, by = c("rownum" = "protein")) %>%
  mutate(p = -log10(p),
         p_mod = -log10(p_mod),
         across(.cols = starts_with("ibaq"), .fns = ~log2(x = .x + 1))) %>%
  rowwise() %>%
  mutate(logFC = (sum(across(contains("test")), na.rm = T) /
                    sum(!near(across(contains("test")), y = 0))) - 
           (sum(across(contains("control")), na.rm = T) /
              sum(!near(across(contains("control")), y = 0)))) %>%
  ungroup() %>%
  mutate(significant_ord = (logFC > 1.5 | logFC < -1.5) & p > 1.30103,
         significant_mod = (logFC > 1.5 | logFC < -1.5) & p_mod > 1.30103)

#Select for mitochondrial proteins
mitocarta <- read_tsv("mitoCarta3.0.tsv") %>%
  as_tibble() %>%
  select(Symbol, Description, MitoCarta3.0_SubMitoLocalization, MitoCarta3.0_MitoPathways)

final_data <- final_data %>%
  inner_join(mitocarta, by = c('Gene names' = 'Symbol')) %>%
  rename(gene_names = 'Gene names')

#Plot graph
final_data %>%
  ggplot(aes(x = logFC, y = p_mod, group = gene_names)) +
  geom_point(aes(color = significant_mod), show.legend = F) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30103, linetype = "dashed") +
  geom_label_repel(aes(label = gene_names), box.padding = 0.5, max.overlaps = 10) +
  theme_classic() +
  labs(x = bquote(Log[2]*" fold change"), y = bquote(-Log[10]*" pvalue_mod"),
       title = str_c(c(test, control), collapse = " v/s "))
ggsave(filename = str_c(test, "_vs_", control, "_moderated_t_test.pdf"))
t <- ggplotly(p = final_data %>%
           ggplot(aes(x = logFC, y = p_mod, group = gene_names)) +
           geom_point(aes(color = significant_mod), show.legend = F) +
           geom_vline(xintercept = 0) +
           geom_vline(xintercept = 1.5, linetype = "dashed") +
           geom_hline(yintercept = 1.30103, linetype = "dashed") +
           theme_classic(), tooltip = c("gene_names", "logFC", "p_mod"))
saveWidget(widget = t, file = str_c(test, "_vs_", control, "_moderated_t_test.html"))

final_data %>%
  ggplot(aes(x = logFC, y = p)) +
  geom_point(aes(color = significant_ord), show.legend = F) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1.5, linetype = "dashed") +
  geom_hline(yintercept = 1.30103, linetype = "dashed") +
  labs(x = bquote(Log[2]*" fold change"), y = bquote(-Log[10]*" pvalue_ord"),
       title = str_c(c(test, control), collapse = " v/s ")) +
  geom_text_repel(aes(label = gene_names), box.padding = 0.5) +
  theme_classic()

ggsave(filename = str_c(test, "_vs_", control, "_ordinary_t_test.pdf"))
#print the graph and final_result file
write_tsv(x = final_data, file = str_c(test, "_vs_", control, "_final_data.tsv"))
#--------------------------------------------------------------------------------


