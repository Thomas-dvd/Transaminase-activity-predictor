import pubchempy as pcp
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, GraphDescriptors, Crippen, MolSurf, Lipinski
import json

noms_molecules = {
    "2OATA": "Oxoacetic Acid",
    "Pyr": "2-Oxopropanoic Acid",
    "2OPAA": "2-OxoPentanoic Acid",
    "2OHAA": "2-OxoHexanoic Acid",
    "4MT2OB": "4-(Methylthio)-2-oxobutanoic Acid",
    "MOPA": "4-Methyl-2-oxopentanoic Acid",
    "tMePyr": "3,3-Dimethyl-2-oxobutanoic acid",
    "2O2PAA": "2-Oxo-2-phenylacetic Acid",
    "2O4PBA": "2-Oxo-4-phenylbutanoic Acid",
    "HPA": "3-Hydroxy-2-oxopropanoic Acid",
    "3H2OBA": "3-Hydroxy-2-oxobutanoic Acid",
    "BOHiPr": "3-Hydroxy-4-methyl-2-oxopentanoic Acid",
    "BOHPh": "3-Hydroxy-2-oxo-3-phenylpropanoic Acid",
    "4H2OA": "4-Hydroxy-2-oxopentanoic Acid",
    "4H4PA": "4-Hydroxy-2-oxo-4-phenylbutanoic Acid",
    "2KDGA": "3,4,5,6-Tetrahydroxy-2-oxohexanoic Acid",
    "3.3diMe": "3,3-Dimethyl-2-oxopentanedioic Acid"
}

dictionnaire_smiles = {}

print("Interrogation de PubChem en cours...\n")

for abbr, vrai_nom in noms_molecules.items():
    try:
        composé = pcp.get_compounds(vrai_nom, 'name')[0]
        dictionnaire_smiles[abbr] = composé.canonical_smiles
        print(f"'{abbr}': '{composé.canonical_smiles}',")
    except IndexError:
        print(f"Erreur 1 sur '{vrai_nom}'.")
    except Exception as e:
        print(f"Erreur 2 sur '{vrai_nom}': {e}")

smarts_tete_keto = Chem.MolFromSmarts("[CX3](=O)[CX3](=O)[OX2H1,OX1-]")
smarts_acide_carboxylique = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
smarts_beta_oh = Chem.MolFromSmarts("[CX4]([OX2H1])[CX3](=O)[CX3](=O)[OX2H1,OX1-]")
smarts_gamma_oh = Chem.MolFromSmarts("[CX4]([OX2H1])[CX4][CX3](=O)[CX3](=O)[OX2H1,OX1-]")

donnees_dream_team = []


print("Calcul des 14 descripteurs de la 'Dream Team'...")

for abbr, smiles in dictionnaire_smiles.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Impossible de lire le SMILES pour {abbr}: {smiles}")

    mol_3d = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol_3d, params) == 0:
        AllChem.MMFFOptimizeMolecule(mol_3d)
        mol_volume = AllChem.ComputeMolVolume(mol_3d)
    else:
        AllChem.Compute2DCoords(mol_3d)
        mol_volume = AllChem.ComputeMolVolume(mol_3d)

    AllChem.ComputeGasteigerCharges(mol)
    
    nb_acides = len(mol.GetSubstructMatches(smarts_acide_carboxylique))
    net_charge_ph7 = -1 * nb_acides

    match_keto = mol.GetSubstructMatch(smarts_tete_keto)
    if match_keto:
        idx_c2 = match_keto[0]
        idx_c1 = match_keto[2]
        atome_c2 = mol.GetAtomWithIdx(idx_c2)
        
        electrophilie_c2 = float(atome_c2.GetProp('_GasteigerCharge'))
        
        voisins_c3 = [
            n for n in atome_c2.GetNeighbors() 
            if n.GetIdx() != idx_c1 and n.GetAtomicNum() == 6
        ]
        
        if len(voisins_c3) == 1:
            atome_c3 = voisins_c3[0]

            degre_sub_c3 = atome_c3.GetDegree() - 1
        else:
            degre_sub_c3 = 0
    else:
        electrophilie_c2 = np.nan
        degre_sub_c3 = 0

    donnees_dream_team.append({
        "Substrat": abbr,
        
        "MolVolume": round(mol_volume, 3),
        "HallKierAlpha": round(GraphDescriptors.HallKierAlpha(mol), 4),
        "Kappa2": round(GraphDescriptors.Kappa2(mol), 4),
        "Kappa3": round(GraphDescriptors.Kappa3(mol), 4),
        
        "MolLogP": round(Crippen.MolLogP(mol), 4),
        "TPSA": round(MolSurf.TPSA(mol), 3),
        "NumHDonors": Lipinski.NumHDonors(mol),
        "NetCharge_pH7": net_charge_ph7,
        
        "NumRotatableBonds": Lipinski.NumRotatableBonds(mol),
        "FractionCSP3": round(Lipinski.FractionCSP3(mol), 4),
        
        "Degre_Sub_C3": degre_sub_c3,
        "Electrophilie_C2": round(electrophilie_c2, 5),
        "Motif_Beta_OH": int(mol.HasSubstructMatch(smarts_beta_oh)),
        "Motif_Gamma_OH": int(mol.HasSubstructMatch(smarts_gamma_oh))
    })

df_dream_team = pd.DataFrame(donnees_dream_team)

nom_fichier = "descripteurs_substrats_01.csv"
df_dream_team.to_csv(nom_fichier, sep=";", decimal=".", index=False)

print(f"\nSuccès ! Fichier '{nom_fichier}' généré avec {df_dream_team.shape[0]} substrats et {df_dream_team.shape[1]-1} descripteurs.")


