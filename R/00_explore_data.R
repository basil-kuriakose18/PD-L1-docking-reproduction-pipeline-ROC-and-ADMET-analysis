
library(readr)
library(dplyr)

table1<-read_csv("data/published/table1_published_affinities.csv")
head(table1)
nrow(table1)
summary(table1$Affinity_kcal_mol)
table1 %>% group_by(Subgroup) %>% summarize(mean_affinity=mean(Affinity_kcal_mol))

library(webchem) 
test_cid <- get_cid("117-39-5", from = "xref/RegistryID", domain = "compound")
print(test_cid)

results_list <- list()
for (i in seq_len(nrow(table1))) {
  cas <- table1$CAS_No[i]
  name <- table1$Compound[i]
  result <- tryCatch({
    cid_lookup <- get_cid(cas, from = "xref/RegistryID", domain = "compound")
    cid <- cid_lookup$cid[1]
    if (is.na(cid)) {
      data.frame(Compound = name, CAS_No = cas, CID = NA, SMILES = NA)
    } else {
      smiles_lookup <- pc_prop(cid, properties = "ConnectivitySMILES")
      data.frame(Compound = name, CAS_No = cas, CID = cid, SMILES = smiles_lookup$ConnectivitySMILES[1])
    }
  }, error = function(e) {
    data.frame(Compound = name, CAS_No = cas, CID = NA, SMILES = NA)
  })
  results_list[[i]] <- result
  Sys.sleep(0.3)
}
structures <- bind_rows(results_list)

failed <- c("Delphinidin", "Malvidin", "Pelargonidin", "Rosinidin")
for (compound_name in failed) {
  cid_lookup <- get_cid(compound_name, from = "name", domain = "compound")
  cid <- cid_lookup$cid[1]
  if (!is.na(cid)) {
    smiles_lookup <- pc_prop(cid, properties = "ConnectivitySMILES")
    structures$CID[structures$Compound == compound_name] <- cid
    structures$SMILES[structures$Compound == compound_name] <- smiles_lookup$ConnectivitySMILES[1]
  }
  Sys.sleep(0.3)
}
write_csv(structures, "data/raw/ligand_structures.csv")
write.table(structures[, c("SMILES", "CAS_No")], "data/raw/ligands.smi", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

structures$pdbqt_file <- paste0("ligand_", seq_len(nrow(structures)), ".pdbqt")
write_csv(structures, "data/raw/ligand_structures.csv")

download.file("https://files.rcsb.org/download/7DY7.pdb", "data/raw/7DY7.pdb")

library(bio3d)
pdb <- read.pdb("data/raw/7DY7.pdb")
table(pdb$atom$resid[pdb$atom$type == "HETATM"])

protein_inds <- atom.select(pdb, "protein")
protein_pdb <- trim.pdb(pdb, protein_inds)
write.pdb(protein_pdb, file = "data/raw/7DY7_receptor_clean.pdb")

ligand_inds <- atom.select(pdb, resid = "HOU")
ligand_pdb <- trim.pdb(pdb, ligand_inds)
write.pdb(ligand_pdb, file = "data/raw/7DY7_native_ligand.pdb")

dir.create("docking/receptor_pdbqt", recursive = TRUE)

dir.create("docking/config", recursive = TRUE)
config_lines <- c(
  "receptor = docking/receptor_pdbqt/7DY7_receptor.pdbqt",
  "ligand = docking/ligands_pdbqt/native_ligand.pdbqt",
  "center_x = 144.2",
  "center_y = -13.2",
  "center_z = 19.8",
  "size_x = 20",
  "size_y = 20",
  "size_z = 20",
  "exhaustiveness = 8",
  "out = results/native_ligand_redock.pdbqt",
  "log = results/native_ligand_redock.log"
)
writeLines(config_lines, "docking/config/redock_config.txt")

library(bio3d)
crystal <- read.pdb("data/raw/7DY7_native_ligand.pdb")
redocked <- read.pdb("results/native_ligand_redock_mode1.pdb")
nrow(crystal$atom)
nrow(redocked$atom)

test <- read.pdb("results/native_ligand_roundtrip.pdb")
nrow(test$atom)

sum(substr(trimws(crystal$atom$elety), 1, 1) == "H")

crystal_heavy_inds <- atom.select(crystal, "noh")
crystal_heavy <- trim.pdb(crystal, crystal_heavy_inds)
redocked_heavy_inds <- atom.select(redocked, "noh")
redocked_heavy <- trim.pdb(redocked, redocked_heavy_inds)
nrow(crystal_heavy$atom)
nrow(redocked_heavy$atom)

crystal_elements <- gsub("[0-9]", "", trimws(crystal_heavy$atom$elety))
redocked_elements <- gsub("[0-9]", "", trimws(redocked_heavy$atom$elety))
identical(crystal_elements, redocked_elements)

table(crystal_elements)
table(redocked_elements)

vina_path <- "C:/Program Files (x86)/The Scripps Research Institute/Vina/Vina.exe"
dock_results <- data.frame(pdbqt_file = character(), best_affinity = numeric(), stringsAsFactors = FALSE)
for (i in 1:3) {
  ligand_file <- paste0("docking/ligands_pdbqt/ligand_", i, ".pdbqt")
  out_file <- paste0("results/ligand_", i, "_docked.pdbqt")
  log_file <- paste0("results/ligand_", i, "_docked.log")
  cmd <- sprintf('"%s" --receptor docking/receptor_pdbqt/7DY7_receptor.pdbqt --ligand %s --center_x 144.2 --center_y -13.2 --center_z 19.8 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 8 --out %s --log %s', vina_path, ligand_file, out_file, log_file)
  system(cmd)
  out_lines <- readLines(out_file)
  result_line <- grep("VINA RESULT", out_lines, value = TRUE)[1]
  affinity <- as.numeric(strsplit(trimws(result_line), "\\s+")[[1]][4])
  dock_results <- rbind(dock_results, data.frame(pdbqt_file = paste0("ligand_", i, ".pdbqt"), best_affinity = affinity))
}
dock_results

dock_results <- data.frame(pdbqt_file = character(), best_affinity = numeric(), stringsAsFactors = FALSE)
for (i in 1:62) {
  cat("Docking ligand", i, "of 62\n")
  ligand_file <- paste0("docking/ligands_pdbqt/ligand_", i, ".pdbqt")
  out_file <- paste0("results/ligand_", i, "_docked.pdbqt")
  log_file <- paste0("results/ligand_", i, "_docked.log")
  cmd <- sprintf('"%s" --receptor docking/receptor_pdbqt/7DY7_receptor.pdbqt --ligand %s --center_x 144.2 --center_y -13.2 --center_z 19.8 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 8 --out %s --log %s', vina_path, ligand_file, out_file, log_file)
  system(cmd)
  out_lines <- readLines(out_file)
  result_line <- grep("VINA RESULT", out_lines, value = TRUE)[1]
  affinity <- as.numeric(strsplit(trimws(result_line), "\\s+")[[1]][4])
  dock_results <- rbind(dock_results, data.frame(pdbqt_file = paste0("ligand_", i, ".pdbqt"), best_affinity = affinity))
}
write_csv(dock_results, "results/docking_scores.csv")

write_csv(structures, "data/raw/ligand_structures.csv")

writeLines(paste(structures$SMILES[45], structures$CAS_No[45]), "data/raw/fix45.smi")
writeLines(paste(structures$SMILES[46], structures$CAS_No[46]), "data/raw/fix46.smi")
writeLines(paste(structures$SMILES[49], structures$CAS_No[49]), "data/raw/fix49.smi")
writeLines(paste(structures$SMILES[50], structures$CAS_No[50]), "data/raw/fix50.smi")

for (n in c(45, 46, 49, 50)) {
  test <- read.pdb(paste0("docking/ligands_pdbqt/ligand_", n, ".pdbqt"))
  cat("ligand_", n, ": ", nrow(test$atom), "atoms\n")
}
for (i in 45:62) {
  cat("Docking ligand", i, "of 62\n")
  ligand_file <- paste0("docking/ligands_pdbqt/ligand_", i, ".pdbqt")
  out_file <- paste0("results/ligand_", i, "_docked.pdbqt")
  log_file <- paste0("results/ligand_", i, "_docked.log")
  cmd <- sprintf('"%s" --receptor docking/receptor_pdbqt/7DY7_receptor.pdbqt --ligand %s --center_x 144.2 --center_y -13.2 --center_z 19.8 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 8 --out %s --log %s', vina_path, ligand_file, out_file, log_file)
  system(cmd)
  out_lines <- readLines(out_file)
  result_line <- grep("VINA RESULT", out_lines, value = TRUE)[1]
  affinity <- as.numeric(strsplit(trimws(result_line), "\\s+")[[1]][4])
  dock_results <- rbind(dock_results, data.frame(pdbqt_file = paste0("ligand_", i, ".pdbqt"), best_affinity = affinity))
  if (i %% 5 == 0) write_csv(dock_results, "results/docking_scores_partial.csv")
}
write_csv(dock_results, "results/docking_scores.csv")

comparison <- dock_results %>%
  left_join(structures %>% select(pdbqt_file, Compound, CAS_No), by = "pdbqt_file") %>%
  left_join(table1 %>% select(CAS_No, Affinity_kcal_mol), by = "CAS_No") %>%
  rename(our_affinity = best_affinity, published_affinity = Affinity_kcal_mol)

library(ggplot2)
ggplot(comparison, aes(x = published_affinity, y = our_affinity)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Published affinity (kcal/mol)", y = "Our affinity (kcal/mol)",
       title = "Docking score reproduction: our Vina runs vs. published Table 1")

comparison <- comparison %>%
  left_join(table1 %>% select(CAS_No, Subgroup), by = "CAS_No")

comparison %>%
  group_by(Subgroup) %>%
  summarise(mean_diff = mean(our_affinity - published_affinity), n = n()) %>%
  arrange(desc(abs(mean_diff)))

library(ggrepel)

top5 <- comparison[order(-comparison$diff), ][1:5, ]

ggplot(comparison, aes(x = published_affinity, y = our_affinity, color = Subgroup)) +
  geom_point(size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_text_repel(data = top5, aes(label = Compound), size = 3, show.legend = FALSE) +
  labs(x = "Published affinity (kcal/mol)", y = "Our affinity (kcal/mol)",
       title = "Docking score reproduction: our Vina runs vs. published Table 1",
       subtitle = "Colored by chemical subgroup; top 5 largest deviations labeled")

dir.create("results/figures", recursive = TRUE)
ggsave("results/figures/affinity_comparison.png", width = 8, height = 6, dpi = 300)

dude_input <- structures[structures$Compound %in% c("BMS-202", "BMS-1166"), c("SMILES", "CAS_No")]
writeLines(paste(dude_input$SMILES, dude_input$CAS_No), "data/raw/dude_actives.smi")

extract_decoys <- function(file_lines) {
  decoy_lines <- file_lines[-1]
  sapply(strsplit(decoy_lines, "\t"), function(x) x[1])
}

decoys1 <- extract_decoys(decoy_file1)
decoys2 <- extract_decoys(decoy_file2)
all_decoys <- unique(c(decoys1, decoys2))
length(all_decoys)

writeLines(paste(all_decoys, paste0("decoy_", seq_along(all_decoys))), "data/raw/decoys.smi")

decoy_dock_results <- data.frame(pdbqt_file = character(), best_affinity = numeric(), stringsAsFactors = FALSE)
for (i in 1:100) {
  cat("Docking decoy", i, "of 100\n")
  result <- tryCatch({
    ligand_file <- paste0("docking/ligands_pdbqt/decoy_", i, ".pdbqt")
    out_file <- paste0("results/decoy_", i, "_docked.pdbqt")
    log_file <- paste0("results/decoy_", i, "_docked.log")
    cmd <- sprintf('"%s" --receptor docking/receptor_pdbqt/7DY7_receptor.pdbqt --ligand %s --center_x 144.2 --center_y -13.2 --center_z 19.8 --size_x 20 --size_y 20 --size_z 20 --exhaustiveness 8 --out %s --log %s', vina_path, ligand_file, out_file, log_file)
    system(cmd)
    out_lines <- readLines(out_file)
    result_line <- grep("VINA RESULT", out_lines, value = TRUE)[1]
    affinity <- as.numeric(strsplit(trimws(result_line), "\\s+")[[1]][4])
    data.frame(pdbqt_file = paste0("decoy_", i, ".pdbqt"), best_affinity = affinity)
  }, error = function(e) {
    data.frame(pdbqt_file = paste0("decoy_", i, ".pdbqt"), best_affinity = NA)
  })
  decoy_dock_results <- rbind(decoy_dock_results, result)
  if (i %% 5 == 0) write_csv(decoy_dock_results, "results/decoy_docking_scores_partial.csv")
}
write_csv(decoy_dock_results, "results/decoy_docking_scores.csv")

library(pROC)

active_scores <- comparison$our_affinity
decoy_scores <- decoy_dock_results$best_affinity[!is.na(decoy_dock_results$best_affinity)]

roc_data <- data.frame(
  score = c(active_scores, decoy_scores),
  label = c(rep(1, length(active_scores)), rep(0, length(decoy_scores)))
)

roc_result <- roc(roc_data$label, -roc_data$score)
auc(roc_result)
plot(roc_result, main = "ROC: Vina score discrimination of PD-L1 actives vs. decoys")

bms_scores <- comparison$our_affinity[comparison$Compound %in% c("BMS-202", "BMS-1166")]
roc_data_bms <- data.frame(
  score = c(bms_scores, decoy_scores),
  label = c(rep(1, length(bms_scores)), rep(0, length(decoy_scores)))
)
roc_result_bms <- roc(roc_data_bms$label, -roc_data_bms$score)
auc(roc_result_bms)

png("results/figures/roc_bms_matched.png", width = 800, height = 800)
plot(roc_result_bms, main = "ROC: Vina discrimination of confirmed PD-L1 inhibitors vs. matched decoys")
dev.off()

admet_subset <- structures[structures$Compound %in% c("Ginkgetin", "Theaflavin", "Diosmin", "Formononetin", "Idaein", "Neohesperidin"), c("Compound", "SMILES")]
admet_subset

writeLines(paste(admet_subset$Compound, admet_subset$SMILES, sep = ": "), "data/raw/admet_smiles_for_submission.txt")

pkcsm_batch <- read_tsv("data/raw/adme_1788681311.64.csv")
swiss <- read_csv("data/raw/swissadme.csv")
ginkgetin_solo <- read_csv("data/raw/ginkgetin_pkcsm_solo.csv")

swiss$Compound <- c("Ginkgetin", "Theaflavin", "Diosmin", "Formononetin", "Idaein", "Neohesperidin")

mw_lookup <- swiss %>% select(Compound, MW) %>% mutate(MW_round = round(MW, 0))
pkcsm_batch <- pkcsm_batch %>% mutate(MW_round = round(MOL_WEIGHT, 0)) %>%
  left_join(mw_lookup %>% select(Compound, MW_round), by = "MW_round")

pkcsm_batch %>% select(Compound, MOL_WEIGHT)

published_admet <- data.frame(
  Compound = c("Ginkgetin", "Theaflavin", "Diosmin", "Formononetin", "Idaein", "Neohesperidin"),
  pub_GI = c(95.38, 65.08, 29.32, 96.12, 45.39, 20.65),
  pub_BBB = c(-1.884, -1.729, -1.795, 0.157, -1.713, -1.720),
  pub_CYP3A4 = c("No", "No", "No", "Yes", "No", "No"),
  pub_OCT2 = c("No", "No", "No", "No", "No", "No"),
  pub_AMES = c("No", "No", "No", "No", "No", "No"),
  pub_hERG1 = c("No", "No", "No", "No", "No", "No"),
  pub_Lipinski = c(0, 3, 3, 0, 2, 3),
  pub_PAINS = c(0, 1, 0, 0, 1, 0)
)

our_pkcsm <- bind_rows(
  ginkgetin_solo %>% select(Compound, GI = Intestinal_absorption, BBB = BBB_permeability, CYP3A4 = CYP3A4_substrate, OCT2 = Renal_OCT2_substrate, AMES = AMES_toxicity, hERG1 = hERG_I_inhibitor),
  pkcsm_batch %>% select(Compound, GI = `Intestinal absorption (human)`, BBB = `BBB permeability`, CYP3A4 = `CYP3A4 substrate`, OCT2 = `Renal OCT2 substrate`, AMES = `AMES toxicity`, hERG1 = `hERG I inhibitor`)
)
our_swiss <- swiss %>% select(Compound, Lipinski = `Lipinski #violations`, PAINS = `PAINS #alerts`)
our_admet <- our_pkcsm %>% left_join(our_swiss, by = "Compound")

admet_comparison <- published_admet %>% left_join(our_admet, by = "Compound")
write_csv(admet_comparison, "results/admet_comparison.csv")
admet_comparison