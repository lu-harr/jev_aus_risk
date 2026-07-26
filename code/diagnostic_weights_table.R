diagnostic_weights_table <- tribble(
  ~diagnostic_type, ~diagnostic_endpoint, ~virus_detection, ~diagnostic_weight,
  "PCR","PCR", "JEV positive", 1,
  "PCR","PCR", "JEV negative", 1,
  "Unassigned", "Unassigned", "JEV negative", 0.8,
  "Serology","panflavi ELISA", "JEV negative", 0.5,
  "Serology","panflavi ELISA", "Inconclusive", 0,
  "Serology", "Unassigned", "JEV negative", 0.25,
  "Serology", "Unassigned", "JEV positive", 0.5, # probably PRNT
  "Serology","ELISA", "JEV negative", 0.25,
  "Serology","ELISA", "JEV positive", 0.25,
  "Serology","ELISA", "Inconclusive", 0,
  "Serology", "Unassigned", 'Inconclusive', 0,
  "Serology","PRNT", "JEV negative", 0.35,
  "Serology","PRNT", "JEV positive", 0.35,
  "Serology","PRNT", "Inconclusive", 0
) %>%
  expand_grid(
    tibble(
      host_type = c("human",
                    "mosquito",
                    "domestic_pig",
                    "feral_pig",
                    "chicken",
                    "other")
    )
  ) %>%
  relocate(
    diagnostic_weight,
    .after = everything()
  )