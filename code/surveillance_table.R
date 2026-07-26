# this snippet has been pinched from secure_code/pathogen_p2_sort_surveillance_data.R
# (some locations were recoverable from that script so won't be able to publish it)
library(xtable)
library(tidyverse)

# should hopefully be able to publish this dataset?
combined_surveillance_df <- read.csv('data/secure/phase_2/clean/surveillance/combined_surveillance_p2p2_dec23.csv')

# for summary of detections in Introduction:
combined_surveillance_df %>%
  filter(virus_detection == "JEV positive") %>%
  group_by(host_type) %>%
  summarise(n = n(), nevents = n_distinct(longitude, latitude, date))


# const used to toggle promotion of sero-positives where -ve PCR tests were recorded
# (also used into sort_surveillance_data script)
PROMOTE_SEROS <- TRUE

tab <- combined_surveillance_df %>%
  filter(virus_detection == "JEV positive" | virus_detection == "JEV negative") %>% 
  mutate(outcome = case_when(virus_detection == "JEV negative" ~ "JEV negative",
                             virus_detection == "JEV positive" & diagnostic_type == "Serology" ~ "JEV antibody positive",
                             virus_detection == "JEV positive" & diagnostic_type == "PCR" ~ "JEV PCR positive",
                             virus_detection == "JEV positive" & diagnostic_type == "Unassigned" & host_type == "chicken" ~ "JEV antibody positive"),
         outcome = factor(outcome,
                          levels = c("JEV negative", "JEV antibody positive", "JEV PCR positive"))) %>%
  # had to single out some South Australian chickens ... may need to get back to those at some point if they're still unassigned ...
  group_by(host_type, state, outcome) %>%
  summarise(event = n_distinct(longitude, latitude, date)) %>%
  pivot_wider(names_from = state,
              values_from = event) %>%
  arrange(factor(host_type, 
                 levels = c("human", "mosquito", "domestic_pig", "chicken", "feral_pig")),
          outcome) %>%
  dplyr::select(c("host_type", "outcome", "NSW", "Vic", "Qld", "SA", "WA", "Tas", "NT")) %>%
  mutate(across(NSW:NT, function(x){
    format(ifelse(is.na(x), 0, x), # chunky handling of NAs
           big.mark = ",",
           trim = TRUE)
  })) %>%
  mutate(outcome = as.character(outcome),
         outcome = ifelse(host_type == "human", "JEV positive", outcome),
         NT = ifelse(host_type == "mosquito" & outcome == "JEV PCR positive", 
                     paste0(NT, "$^\\dag$"), NT)) %>%
  #mutate(outcome = paste0('\  \  \  \  ', outcome)) %>%
  #rename(Host = outcome) %>%
  mutate(Qld = case_when(Qld == "0" ~ "-",
                         TRUE ~ Qld),
         Tas = case_when(Tas == "0" & host_type %in% c("mosquito", "chicken") ~ "-",
                         TRUE ~ Tas),
         Vic = case_when(Vic == "0" & host_type == "chicken" ~ "-",
                         TRUE ~ Vic)) %>%
  rename(Host = host_type,
         " "= outcome) %>%
  ungroup() %>%
  mutate(Host = "\ \ \ \ ") %>%
  xtable::xtable()

host_multicols <- c("\\textbf{Humans} - PCR and serological test data",
                    "\\textbf{Mosquitoes} - PCR test data",
                    "\\textbf{Domestic pigs} - PCR and serological test data",
                    "\\textbf{Sentinel chicken flocks} - Serological test data",
                    "\\textbf{Feral pigs} - PCR and serological test data")
addtorow <- list()
addtorow$pos <- list(0, 1, 3, 6, 8)
addtorow$command <- sapply(host_multicols, function(x){
  paste0('\\midrule\\addlinespace[2.5pt] \\multicolumn{9}{l}{', x, '}\\\\', collapse='')})

print(tab,
      add.to.row = addtorow,
      include.rownames=FALSE,
      sanitize.colnames.function = function(x){paste0("\\textbf{", x, "}")},
      sanitize.text.function = function(x){x},
      hline.after = NULL,
      file = paste0("surveillance_table/sero_promo_", ifelse(PROMOTE_SEROS, "true", "false"), ".tex"))

