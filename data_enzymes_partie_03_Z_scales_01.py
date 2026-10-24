import pandas as pd

df = pd.read_csv("data_exploitable_python_01.csv", sep=";") 

with open("mes_enzymes_a_aligner.fasta", "w") as f:
    for _, row in df.iterrows():
        f.write(f">{row['Acc. Number']}\n{row['Sequence']}\n")
print("Fichier FASTA créé avec succès !")