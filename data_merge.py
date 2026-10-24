import pandas as pd

df_enzymes = pd.read_csv("data_exploitable_python_04 (ESM2, Z-Scales, Biopython).csv", sep=";")
df_substrats = pd.read_csv("descripteurs_substrats_01.csv", sep=";") 

substrats = [
    "2OATA", "Pyr", "2OPAA", "2OHAA", "4MT2OB", "MOPA", "tMePyr", "2O2PAA", "2O4PBA",
    "HPA", "3H2OBA", "BOHiPr", "BOHPh", "4H2OA", "4H4PA", "2KDGA", "3.3diMe"
]

colonnes_id_vars = [col for col in df_enzymes.columns if col not in substrats]

df_long = pd.melt(
    df_enzymes,
    id_vars=colonnes_id_vars,
    value_vars=substrats,
    var_name="Substrat",
    value_name="Activite_Brute"
)

df_substrats = df_substrats.rename(columns={"Substrate": "Substrat", "Nom": "Substrat"}, errors="ignore")

# Croisement (Merge)
df_final = pd.merge(df_long, df_substrats, on="Substrat", how="inner")

print("\n--- Bilan du dataset brut pour le Machine Learning ---")
print(f"Lignes totales : {len(df_final)} ({len(df_enzymes)} enzymes x 17 substrats)")
print(f"Colonnes : {len(df_final.columns)}")

df_final.to_csv("data_exploitable_python_05 (merged).csv", index=False, sep=";")
print("\n✅ Fichier 'data_exploitable_python_05 (merged).csv' généré avec succès !")