#Import libraries
library(tidyverse)
library(infer)
library(tidymodels)
library(readxl)
library(ggrepel)
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

# select proteins with unique peptide > 1 in two out of three samples in atleast 
# one experimental group
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

# Select proteins with valid values in 2 out of 3 replicates in both experimental
# groups
valid_value_names <- data %>%
        select(rownum, matches(match = "^iBAQ")) %>%
        mutate(across(.cols = -rownum, .fns = ~near(x = .x, y = 0)),
               across(.cols = -rownum, .fns = as.integer)) %>%
        rowwise() %>%
        mutate(test = sum(across(.cols = contains(test))),
               control = sum(across(.cols = contains(control)))) %>%
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

# Calculate ordinary t-statistic and the p_value
final_result <- names_for_analysis$rownum %>%
        map_dfr(.f = ~ord_t_test(data = data_preprocessed_for_t_test,
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
           across(.cols = contains("iBAQ"), .fns = ~log2(x = .x + 1))) %>%
    rowwise() %>%
    mutate(logFC = (sum(across(contains(test)), na.rm = T) /
                      sum(!near(across(contains(test)), y = 0))) - 
             (sum(across(contains(control)), na.rm = T) /
                sum(!near(across(contains(control)), y = 0)))) %>%
    ungroup() %>%
        mutate(significant_ord = (logFC > 1.5 | logFC < -1.5) & p > 1.30103,
               significant_mod = (logFC > 1.5 | logFC < -1.5) & p_mod > 1.30103)

#Select for mitochondrial proteins
mitocarta <- read_xls(path = "Human.MitoCarta3.0.xls",
                      sheet = 2) %>%
        as_tibble() %>%
        select(Symbol, Description, MitoCarta3.0_SubMitoLocalization, MitoCarta3.0_MitoPathways)

 final_data <- final_data %>%
        inner_join(mitocarta, by = c('Gene names' = 'Symbol')) %>%
         rename(gene_names = 'Gene names')

#Plot graph
final_data %>%
        ggplot(aes(x = logFC, y = p_mod)) +
        geom_point(aes(color = significant_mod), show.legend = F) +
        geom_vline(xintercept = 0) +
        geom_vline(xintercept = 1.5, linetype = "dashed") +
        geom_hline(yintercept = 1.30103, linetype = "dashed") +
        labs(x = bquote(Log[2]*" fold change"), y = bquote(-Log[10]*" pvalue"),
             title = str_c(c(test, control), collapse = " v/s ")) +
        geom_text_repel(aes(label = gene_names), box.padding = 0.5) +
        theme_classic()
ggsave(filename = str_c(test, "_vs_", control, "_moderated_t_test.pdf"))


final_data %>%
        ggplot(aes(x = logFC, y = p)) +
        geom_point(aes(color = significant_ord), show.legend = F) +
        geom_vline(xintercept = 0) +
        geom_vline(xintercept = 1.5, linetype = "dashed") +
        geom_hline(yintercept = 1.30103, linetype = "dashed") +
        labs(x = bquote(Log[2]*" fold change"), y = bquote(-Log[10]*" pvalue"),
             title = str_c(c(test, control), collapse = " v/s ")) +
        geom_text_repel(aes(label = gene_names), box.padding = 0.5) +
        theme_classic()

#print the graph and final_result file
ggsave(filename = str_c(test, "_vs_", control, "_ordinary_t_test.pdf"))
write_tsv(x = final_data, file = str_c(test, "_vs_", control, "_final_data.tsv"))
#--------------------------------------------------------------------------------


