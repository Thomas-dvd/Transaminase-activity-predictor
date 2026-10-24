import re
import pandas as pd



def charger_fasta(chemin_fasta):
    sequences = {}
    acc_actuel = None
    seq_actuelle = []

    with open(chemin_fasta, "r", encoding="utf-8") as f:
        for ligne in f:
            ligne = ligne.strip()
            if not ligne:
                continue

            if ligne.startswith(">"):
                if acc_actuel is not None:
                    sequences[acc_actuel] = "".join(seq_actuelle)

                parties = ligne.split("|")
                if len(parties) >= 2:
                    acc_actuel = parties[1].strip()
                else:
                    acc_actuel = ligne[1:].split()[0]

                seq_actuelle = []

            else:
                if re.match(r"^[A-Z\*]+$", ligne):
                    seq_actuelle.append(ligne)

        if acc_actuel is not None:
            sequences[acc_actuel] = "".join(seq_actuelle)

    return sequences

fichier_fasta = "FASTAs_TAMA.txt"
fichier_csv = "Dataset.csv"
fichier_sortie = "data_exploitable_python_01.csv"

dico_sequences = charger_fasta(fichier_fasta)
print(f"{len(dico_sequences)} séquences chargées depuis le fichier FASTA.")

df = pd.read_csv(fichier_csv, sep=";", decimal=",")

df["Acc. Number"] = df["Acc. Number"].astype(str).str.strip()

df["Sequence"] = df["Acc. Number"].map(dico_sequences)

manquantes = df["Sequence"].isna().sum()
if manquantes > 0:
    print(
        f"Attention : {manquantes} protéine(s) du CSV n'ont pas été trouvées"
        " dans le FASTA."
    )
else:
    print("Toutes les séquences ont été associées avec succès !")

df.to_csv(fichier_sortie, sep=";", index=False, encoding="utf-8")
print(f"Fichier sauvegardé sous : {fichier_sortie}")