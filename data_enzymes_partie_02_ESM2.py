import pandas as pd
import torch
import numpy as np
from transformers import AutoTokenizer, AutoModel
from sklearn.decomposition import PCA
from tqdm import tqdm # Pour la barre de progression



# PARAMÈTRES ET CHARGEMENT DU MODÈLE
model_name = "facebook/esm2_t12_35M_UR50D"

print("Chargement du modèle ESM-2 (cela peut prendre quelques secondes)...")
tokenizer = AutoTokenizer.from_pretrained(model_name)
model = AutoModel.from_pretrained(model_name)
model.eval()



# FONCTION D'EXTRACTION (MEAN POOLING)
def get_esm_embedding(sequence):
    # Nettoyer la séquence au cas où
    sequence = sequence.replace("*", "").upper()
    
    inputs = tokenizer(sequence, return_tensors="pt")
    
    with torch.no_grad(): 
        outputs = model(**inputs)
    
    last_hidden_states = outputs.last_hidden_state
    

    sequence_embedding = last_hidden_states[0, 1:-1, :].mean(dim=0).numpy()
    
    return sequence_embedding



# LECTURE DU CSV ET CALCUL DES VECTEURS
df_enzymes = pd.read_csv("data_exploitable_python_01.csv", sep=";")
nom_colonne_sequence = "Sequence"

print(f"Calcul des embeddings pour {len(df_enzymes)} enzymes...")
embeddings = []

# For fun
for seq in tqdm(df_enzymes[nom_colonne_sequence]):
    vec = get_esm_embedding(seq)
    embeddings.append(vec)

X_embeddings = np.array(embeddings) 



# PCA 480 -> n_composantes
print("\nPCA...")

n_composantes = 10
pca = PCA(n_components=n_composantes)
X_pca = pca.fit_transform(X_embeddings)

variance_conservee = sum(pca.explained_variance_ratio_) * 100
print(f"La PCA à {n_composantes} composantes conserve {variance_conservee:.2f}% de l'information évolutive !")



# AJOUT AU DATASET ACTUEL ET SAUVEGARDE
pca_cols = [f"ESM_PC{i+1}" for i in range(n_composantes)]
df_pca = pd.DataFrame(X_pca, columns=pca_cols)

df_resultat = pd.concat([df_enzymes, df_pca], axis=1)

df_resultat.to_csv("data_exploitable_python_02 (ESM2).csv", index=False, sep=";")
print(f"Le fichier 'data_exploitable_python_02 (ESM2).csv' contient désormais {n_composantes} colonnes ESM.")