"""Cage module."""

import numpy as np
import rdkit.Chem.AllChem as rdkit
import stk
from rdkit import Chem


class Cage:
    def __init__(self, file_path=None, molecule=None):
        """Initialize a Cage object either from a file or directly from an stk molecule.

        Args:
        ----
        - file_path (str): Path to the file from which to load the cage.
        - molecule (stk.Molecule): Directly provided stk molecule.

        """
        if file_path:
            self.cage = stk.BuildingBlock.init_from_file(file_path)
        elif molecule:
            self.cage = molecule
        else:
            msg = "Either file_path or molecule must be provided."
            raise ValueError(msg)

    def calculate_cage_vector(self, arene_smile):
        """Calculate the vector representation of a cage molecule.

        Parameters
        ----------
        - stk_mol: An STK molecular object representing the cage.
        - arene_smile: SMILES string for the arene part of the molecule.

        Returns
        -------
        - A vector representation of the cage molecule.

        """
        rdkit_mol = self.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        arenes = [
            self.get_centroid(atom_ids=atom_ids)
            for atom_ids in rdkit_mol.GetSubstructMatches(
                query=rdkit.MolFromSmarts(arene_smile)
            )
        ]
        arenes = np.asarray(arenes)
        centroid_mol = self.get_centroid()
        vectors = np.array(
            [
                (arene - centroid_mol) / np.linalg.norm(arene - centroid_mol)
                for arene in arenes
            ]
        )
        return vectors[0]

    def cage_size(self, arene_smile):
        """Calculate the average distance between the centroid of a cage and its arenes.

        Parameters
        ----------
        - cage: The stk cage molecule.
        - arene_smile: SMILES string representing the arene component.

        Returns
        -------
        - The average distance between the centroid of the cage and its arenes.

        """
        rdkit_mol = self.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        arenes = [
            self.get_centroid(atom_ids=atom_ids)
            for atom_ids in rdkit_mol.GetSubstructMatches(
                query=rdkit.MolFromSmarts(arene_smile)
            )
        ]
        centroid_mol = np.mean(arenes, axis=0)
        distances = [np.linalg.norm(arene - centroid_mol) for arene in arenes]
        return np.mean(distances)

    def fix_atom_set(self, diamine_smiles, metal_atom=None):
        """Return atom IDs of carbons in diamine adjacent to amine nitrogens.

        Args:
        ----
        - mol_file_path (str): Path to the .mol file.
        - diamine_smiles (str): SMILES string of the diamine.
        Default is "NC1CCCCC1N".

        Returns:
        -------
        - List[int]: List of atom IDs.

        """
        # Use the SMILES string of the diamine to identify its structure
        diamine = Chem.MolFromSmiles(diamine_smiles)
        rdkit_mol = self.to_rdkit_mol()
        # Get the substructure matches for the diamine
        matches = rdkit_mol.GetSubstructMatches(diamine)
        # List to store carbon atom ids
        fixed_atom_ids = []
        # Iterate over the matches
        if metal_atom is not None:  # Checks if metal_ids is not empty
            for atom in rdkit_mol.GetAtoms():
                # Check if the atom symbol matches metal_atom
                if atom.GetSymbol() == metal_atom:
                    # Append the atom's index if it matches
                    fixed_atom_ids.append(atom.GetIdx())
        else:
            for match in matches:
                for atom_id in match:
                    atom = rdkit_mol.GetAtomWithIdx(atom_id)
                    # Check if the atom is carbon and adjacent to nitrogen
                    if atom.GetSymbol() == "C" and any(
                        neighbor.GetSymbol() == "N"
                        for neighbor in atom.GetNeighbors()
                    ):
                        fixed_atom_ids.append(atom_id)
        return fixed_atom_ids
