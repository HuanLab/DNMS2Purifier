######################################
# Structural similarity calculation # 
######################################
library(rJava)
library(rcdk)
library(fingerprint)
#------------main function: two SMILES ----------------------
structureSmi <- function (smiles_A, smiles_B) {
        # SMILES vector
        smiles_v <- c(smiles_A, smiles_B)
        mols <- parse.smiles(smiles_v)
        # fingerprint types, commonly used: PubChem, ECFP4, ECFP6, etc.
        type <- c("standard", "extended", "graph", "hybridization",
                  "maccs", "estate",
                  "pubchem",
                  "kr", "shortestpath","signature", "circular","substructure")
        fps <- lapply(mols, get.fingerprint, type=type[7])
        output <- fp.sim.matrix(fps, method='tanimoto')
        # output: structural similarity score matrix
        return (output[2,1])
}