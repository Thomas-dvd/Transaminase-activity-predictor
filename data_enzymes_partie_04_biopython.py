import pandas as pd
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis

nom_fichier = "data_exploitable_python_03 (ESM2, Z-Scales).csv"
df = pd.read_csv(nom_fichier, sep=";")
df.columns = df.columns.str.strip()

nom_colonne_sequence = "Sequence"

def calculer_biopython_4(seq_brute):
    seq_propre = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', str(seq_brute).upper())
    longueur = len(seq_propre)
    
    if longueur == 0:
        return pd.Series([0.0, 0.0, 0.0, 0.0])
        
    pa = ProteinAnalysis(seq_propre)
    
    aa_counts = pa.count_amino_acids()
    
    aliphatic_idx = (
        (aa_counts.get('A', 0) / longueur) + 
        2.9 * (aa_counts.get('V', 0) / longueur) + 
        3.9 * ((aa_counts.get('I', 0) + aa_counts.get('L', 0)) / longueur)
    ) * 100
    
    return pd.Series([
        pa.isoelectric_point(),
        pa.gravy(),
        pa.instability_index(),
        aliphatic_idx
    ])


print("Calcul des 4 descripteurs Biopython...")
colonnes_bio = ["Bio_pI", "Bio_GRAVY", "Bio_Instability", "Bio_Aliphatic"]
df[colonnes_bio] = df[nom_colonne_sequence].apply(calculer_biopython_4)

df.to_csv("data_exploitable_python_04 (ESM2, Z-Scales, Biopython).csv", index=False, sep=";")
print("✅ Terminé ! Les 4 colonnes Biopython ont été ajoutées avec succès.")
print(df[colonnes_bio].head())