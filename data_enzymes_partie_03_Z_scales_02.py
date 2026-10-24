from Bio import SeqIO
import pandas as pd

fichier_alignement = "Sequences_alignes.fasta"
id_reference = "C7RAM1"               

vraies_positions_reference = [42, 76, 77, 102, 153, 173, 196, 239, 240]

def trouver_colonnes_alignement(fichier, id_ref, vraies_positions):
    colonnes_trouvees = []
    
    for record in SeqIO.parse(fichier, "fasta"):
        if id_ref in record.id:
            seq_alignee = str(record.seq)
            compteur_vrais_aa = 0
            
            for index_colonne, char in enumerate(seq_alignee):
                if char != '-':
                    compteur_vrais_aa += 1
                    
                    if compteur_vrais_aa in vraies_positions:
                        colonnes_trouvees.append(index_colonne)
                        
            return colonnes_trouvees
            
    print("Erreur : ID de référence introuvable !")
    return []


colonnes_cles = trouver_colonnes_alignement(fichier_alignement, id_reference, vraies_positions_reference)
print(f"Les positions {vraies_positions_reference} de {id_reference} correspondent aux colonnes {colonnes_cles} de l'alignement.")
print(f"Nombre de colonnes trouvées : {len(colonnes_cles)} (doit être égal à 9 !)")
print(f"Colonnes dans l'alignement : {colonnes_cles}")

Z_SCALES = {
    'A': (0.24, -2.32, 0.60), 'R': (3.52, 2.50, -3.50), 'N': (3.05, 1.62, 1.04),
    'D': (3.98, 0.93, 1.93), 'C': (0.84, -1.67, 3.71), 'Q': (1.75, 0.50, -1.44),
    'E': (3.11, 0.26, -0.11), 'G': (2.05, -4.06, 0.36), 'H': (2.47, 1.95, 0.26),
    'I': (-3.89, -1.73, -1.71), 'L': (-4.28, -1.30, -1.49), 'K': (2.29, 0.89, -2.49),
    'M': (-2.85, -0.22, 0.47), 'F': (-4.22, 1.94, 1.06), 'P': (-1.66, 0.27, 1.84),
    'S': (2.39, -1.07, 1.15), 'T': (0.75, -2.18, -1.12), 'W': (-4.36, 3.94, 0.59),
    'Y': (-2.54, 2.44, 0.43), 'V': (-2.59, -2.64, -1.54),
    '-': (0.00, -5.00, 0.00), 'X': (0.00, 0.00, 0.00)
}

resultats = []

for record in SeqIO.parse(fichier_alignement, "fasta"):
    features = {"Acc. Number": record.id}
    seq_alignee = str(record.seq).upper()
    
    for col in colonnes_cles:
        aa = seq_alignee[col]
        z1, z2, z3 = Z_SCALES.get(aa, (0.0, 0.0, 0.0))
        
        features[f"Col{col}_Z1_Pol"] = z1
        features[f"Col{col}_Z2_Vol"] = z2
        features[f"Col{col}_Z3_Elec"] = z3
        
    resultats.append(features)

df_zscales = pd.DataFrame(resultats)
print(f"Extraction terminée pour {len(df_zscales)} enzymes. Voici un aperçu :")
print(df_zscales.head())

df_precedent = pd.read_csv("data_exploitable_python_02 (ESM2).csv", sep=";")
df_final = pd.merge(df_precedent, df_zscales, on="Acc. Number", how="left")
df_final.to_csv("data_exploitable_python_03 (ESM2, Z-Scales).csv", index=False, sep=";")
print("Fichier final sauvegardé !")