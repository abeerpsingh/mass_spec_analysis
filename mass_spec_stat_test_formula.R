#Import libraries
library(tidyverse)
library(infer)
library(tidymodels)
library(readxl)
library(ggrepel)

#Functions-------------------------------------------------------------------------
estimate_t_test <- function(data, protein = NA){
        data <- data %>%
                select(name, protein) %>%
                rename(resp = protein)
        
        point_estimate <- data %>%
                specify(response = resp, explanatory = name) %>%
                calculate(stat = "diff in means", order = c("test", "control"))
        
        if(point_estimate[["stat"]] > 0){
                tibble(protein = protein) %>%
                        bind_cols(t_test(x = data, formula = resp~name,
                                         order = c("test", "control"),
                                         alternative = "greater", conf_int = F))     
        }
        else {
                tibble(protein = protein) %>%
                        bind_cols(t_test(x = data, formula = resp~name,
                                         order = c("test", "control"),
                                         alternative = "less", conf_int = F))    
        }
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
test <- readline(prompt = "Enter name of test sample:")
control <- readline(prompt = "Enter name of control sample: ")
rep <- readline(prompt = "Number of replicates for test or control: ") %>%
        parse_integer()
data <- protein_groups %>%
        select("rownum", "Gene names", contains(test, ignore.case = T),
               contains(control, ignore.case = T))

# Select proteins: sequence coverage > 5 for atleast 2 out of 3 samples in atleast
# one experimental group
names_sequence_coverage <- data %>%
        select(rownum, starts_with("Sequence coverage")) %>%
        mutate(across(.cols = -rownum, .fns = ~between(x = .x, left = 0, right = 5))) %>%
        mutate(across(.cols = -rownum, .fns = as.integer)) %>%
        rowwise() %>%
        mutate(test = sum(across(.cols = contains(test))),
               control = sum(across(.cols = contains(control)))) %>%
        ungroup() %>%
        filter(test < 2 | control < 2) %>%
        select(rownum) %>%
        drop_na()

# select proteins with unique peptide > 1
names_unique_peptides <- data %>%
        select(rownum, matches(match = "^Unique")) %>%
        mutate(across(.cols = -rownum, .fns = ~between(x = .x, left = 0,
                                                       right = 1)),
               across(.cols = -rownum, .fns = as.integer)) %>%
        rowwise() %>%
        mutate(test = sum(across(.cols = contains(test))),
               control = sum(across(.cols = contains(control)))) %>%
        ungroup() %>%
        filter(test < 2 | control < 2) %>%
        select(rownum) %>%
        drop_na()

#Select proteins with valid values in 2 out of 3 replicates
valid_value_names <- data %>%
        select(rownum, matches(match = "^iBAQ")) %>%
        mutate(across(.cols = -rownum, .fns = ~near(x = .x, y = 0)),
               across(.cols = -rownum, .fns = as.integer)) %>%
        rowwise() %>%
        mutate(test = sum(across(.cols = contains(test))),
               control = sum(across(.cols = contains(control)))) %>%
        ungroup() %>%
        filter(test < 2 | control < 2) %>%
        select(rownum) %>%
        drop_na()

#Inner join the selected names: right join
names_for_analysis <- valid_value_names %>%
        inner_join(names_sequence_coverage, by = "rownum") %>%
        inner_join(names_unique_peptides, by = "rownum")

#filter out the selected names for analysis from the data: right join
data_combined <- right_join(x = data, y = names_for_analysis, by = "rownum")
data_iBAQ <- data_combined %>%
        select(rownum, contains("iBAQ"), 'Gene names')
data_transposed <- data_combined %>%
        select(rownum, contains("iBAQ")) %>%
        pivot_longer(-rownum) %>%
        pivot_wider(names_from = rownum) %>%
        mutate(name = if_else(str_detect(name, test), true = "test", 
                              false = "control"))

#Preprocess data: Imputation and scale
data_transpose_impute <- data_transposed %>%
        recipe() %>%
        step_impute_lower(all_numeric()) %>%
        prep(verbose = T, log_changes = T)

data_preprocessed_for_t_test <- data_transpose_impute[["template"]]

# Calculate p-value by standard t-test
final_result <- names_for_analysis$rownum %>%
        map_dfr(.f = ~estimate_t_test(data = data_preprocessed_for_t_test, protein = .x))

#Determine significance
final_result <- data_iBAQ %>%
        right_join(y = final_result, by = c("rownum" = "protein")) %>%
        mutate(p_value = -log10(p_value),
               across(.cols = contains("iBAQ"), .fns = ~log2(x = .x + 1))) %>%
        rowwise() %>%
  mutate(logFC = (sum(across(contains(test)), na.rm = T) /
                    sum(!near(across(contains(test)), y = 0))) - 
           (sum(across(contains(control)), na.rm = T) /
              sum(!near(across(contains(control)), y = 0)))) %>%
        ungroup() %>%
        mutate(significant = (logFC > 1.5 | logFC < -1.5) & p_value > 1.30103)

#Select for mitochondrial proteins
mitocarta <- read_xls(path = "Human.MitoCarta3.0.xls",
                      sheet = 2) %>%
        as_tibble() %>%
        select(Symbol, Description, MitoCarta3.0_SubMitoLocalization, MitoCarta3.0_MitoPathways)

final_result <- final_result %>%
        inner_join(mitocarta, by = c('Gene names' = 'Symbol')) %>%
        select(-statistic, t_df, estimate) %>%
  rename(gene_names = 'Gene names')

#Plot graph
final_result %>%
        ggplot(aes(x = logFC, y = p_value)) +
        geom_point(aes(color = significant), show.legend = F) +
        geom_vline(xintercept = 0) +
        geom_vline(xintercept = 1.5, linetype = "dashed") +
        geom_hline(yintercept = 1.30103, linetype = "dashed") +
        labs(x = bquote(Log[2]*" fold change"), y = bquote(-Log[10]*" pvalue"),
             title = str_c(c(test, control), collapse = " v/s ")) +
        geom_text_repel(aes(label = gene_names), box.padding = 0.5) +
        theme_classic()

#print the graph and final_result file
ggsave(filename = str_c(test, "_vs_", control, "_ordinary_t_test.pdf"))
write_tsv(x = final_result, file = str_c(test, "_vs_", control, "_final_data.tsv"))
#--------------------------------------------------------------------------------

