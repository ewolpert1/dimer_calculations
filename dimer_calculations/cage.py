# Assuming stk and other necessary imports are already handled in this file
import itertools

import numpy as np
import rdkit.Chem.AllChem as rdkit
import stk
from rdkit import Chem

from . import utils


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
            raise ValueError("Either file_path or molecule must be provided.")

    def construct(self, bb1_smiles, bb2_smiles):
        """Given the SMILES strings of two building blocks,
        optimize both building blocks and construct a cage with FourPlusSix topology.
        The resulting STK molecule of the cage is centered at the origin.
        """
        bb1 = stk.BuildingBlock(
            smiles=bb1_smiles,
            functional_groups=[
                stk.PrimaryAminoFactory(),
                stk.AldehydeFactory(),
            ],
        )
        rdkit_bb1 = bb1.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_bb1)
        rdkit.MMFFOptimizeMolecule(rdkit_bb1)
        bb1 = bb1.with_position_matrix(
            position_matrix=rdkit_bb1.GetConformer().GetPositions(),
        )

        bb2 = stk.BuildingBlock(
            smiles=bb2_smiles,
            functional_groups=[
                stk.AldehydeFactory(),
                stk.PrimaryAminoFactory(),
            ],
        )
        rdkit_bb2 = bb2.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_bb2)
        rdkit.MMFFOptimizeMolecule(rdkit_bb2)
        bb2 = bb2.with_position_matrix(
            position_matrix=rdkit_bb2.GetConformer().GetPositions(),
        )

        # Construct the cage
        self.cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
            optimizer=stk.MCHammer(),
        )
        self.cage = self.cage.with_centroid([0, 0, 0])

    # Add other methods as necessary, like optimize, calculate vectors, etc.

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

    def generate_new_cages(
        self,
        guest_bb_aa,
        guest_bb_ww,
        guest_bb_wa,
        guest_bb_wv,
        direction_vector_w,
        direction_vector_a,
        rotated_vectors_w,
        rotated_vectors_a,
        displacement_distance,
    ):
        """Generate new cage complexes with various packings and rotations.

        Parameters
        ----------
        - cage: The base cage molecule.
        - guest_bb_aa: Guest molecule with arene-to-arene packing.
        - guest_bb_ww: Guest molecule with window-to-window packing.
        - guest_bb_wa: Guest molecule with window-to-arene packing.
        - direction_vector: Direction vector for displacement.
        - perpendicular_vector: Vector perpendicular to the direction vector.
        - rotated_vectors: A list of vectors representing all the rotation directions.
        - displacement_distance: The distance for displacement.

        Returns
        -------
        - Lists of constructed complexes with different packings and rotations.
        - A list of centroids for all displaced cages.

        """
        perpendicular_vector_w = utils.calculate_perpendicular_vector(
            direction_vector_w
        )
        perpendicular_vector_a = utils.calculate_perpendicular_vector(
            direction_vector_a
        )
        list_wa, list_ww, list_aa, list_wv = [], [], [], []
        for i in np.arange(0, 7, 1):
            displaced_centers_w = utils.find_integer_points(
                direction_vector_w,
                displacement_distance * 2 - 2 + i,
                int(displacement_distance) + 1,
            )
            displaced_centers_a = utils.find_integer_points(
                direction_vector_a,
                displacement_distance * 2 - 2 + i,
                int(displacement_distance) + 1,
            )
            list_wa_cent, list_ww_cent, list_aa_cent, list_wv_cent = (
                [],
                [],
                [],
                [],
            )
            for center in range(len(displaced_centers_w)):
                list_wa_rot, list_ww_rot, list_aa_rot, list_wv_rot = (
                    [],
                    [],
                    [],
                    [],
                )
                for rotation_vector in range(len(rotated_vectors_w)):
                    guest_aa = utils.create_rotated_guest(
                        guest_bb_aa,
                        perpendicular_vector_a,
                        rotated_vectors_a[rotation_vector],
                        displaced_centers_a[center],
                    )
                    guest_wa = utils.create_rotated_guest(
                        guest_bb_wa,
                        perpendicular_vector_w,
                        rotated_vectors_w[rotation_vector],
                        displaced_centers_w[center],
                    )
                    guest_wv = utils.create_rotated_guest(
                        guest_bb_wv,
                        perpendicular_vector_w,
                        rotated_vectors_w[rotation_vector],
                        displaced_centers_w[center],
                    )
                    guest_ww = utils.create_rotated_guest(
                        guest_bb_ww,
                        perpendicular_vector_w,
                        rotated_vectors_w[rotation_vector],
                        displaced_centers_w[center] * -1,
                    )
                    complex_wa = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(
                            host=self, guests=guest_wa
                        )
                    )
                    complex_ww = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(
                            host=self, guests=guest_ww
                        )
                    )
                    complex_aa = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(
                            host=self, guests=guest_aa
                        )
                    )
                    complex_wv = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(
                            host=self, guests=guest_wv
                        )
                    )
                    list_wa_rot.append(complex_wa)
                    list_ww_rot.append(complex_ww)
                    list_aa_rot.append(complex_aa)
                    list_wv_rot.append(complex_wv)
                list_wa_cent.append(list_wa_rot)
                list_ww_cent.append(list_ww_rot)
                list_aa_cent.append(list_aa_rot)
                list_wv_cent.append(list_wv_rot)
            list_wa.append(list_wa_cent)
            list_ww.append(list_ww_cent)
            list_aa.append(list_aa_cent)
            list_wv.append(list_wv_cent)
        return list_wa, list_ww, list_aa, list_wv

    def generate_ww_cages(
        cage, guest_bb_ww, vec, vec_perpendicular, rotated_vectors, dist
    ):
        """Args:
        ----
        - cage: constructed molecule of the cage
        - guest_bb_aa: guest molecule with arene to arene packing
        - guest_bb_ww: guest molecule with window to window packing
        - guest_bb_wa: guest molecule with window to arene packing
        - vec: direction vector of displacement
        - vec_perpendicular: orthogonal vector of vec
        - rotated_vectors: a list of rotated vectors
        - dist: magnitude of minimum displacement
        Returns:
        - lists of constructed dimers with packings of ww, wa, aa and rotations
        - list_centroids: a list of vectors representing the centroids of all the displaced cages

        """

    def fix_atom_set(self, diamine_smiles, metal_atom=None):
        """Return the atom IDs of the carbon atoms from the diamine that are adjacent to the amine nitrogen atoms.

        Args:
        ----
        - mol_file_path (str): Path to the .mol file.
        - diamine_smiles (str): SMILES string of the diamine. Default is "NC1CCCCC1N".

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
                        [
                            neighbor.GetSymbol() == "N"
                            for neighbor in atom.GetNeighbors()
                        ]
                    ):
                        fixed_atom_ids.append(atom_id)
        return fixed_atom_ids

    def cage_vertex_vec(self, amine_smile):
        """Returns the vector of average distance between the centroid of the cage and arenes"""
        # pdb_file = "%s.pdb"%(cagename)
        rdkit_mol = self.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        arenes = []
        for atom_ids in rdkit_mol.GetSubstructMatches(
            query=rdkit.MolFromSmarts(amine_smile),
        ):
            centroid = self.get_centroid(atom_ids=atom_ids)
            arenes.append(centroid)
        arenes = np.asarray(arenes)
        centroid_mol = self.get_centroid()
        vec = np.zeros((len(arenes), 3))
        for i in range(len(arenes)):
            vec[i] = arenes[i] - centroid_mol
            vec[i] = vec[i] / (np.sqrt(np.dot(vec[i], vec[i])))
        return vec

    def find_midpoint_of_facet(self, vectors, sym, tolerance=0.1):
        # Calculate all distances and find the minimum non-zero distance
        all_distances = [
            utils.distance(v1, v2)
            for v1, v2 in itertools.combinations(vectors, 2)
        ]
        min_distance = min(filter(lambda d: d > 0, all_distances))
        # Iterate through all combinations of three vectors
        if sym == "oct":
            for v1, v2, v3 in itertools.combinations(vectors, 3):
                # Calculate distances between the vectors
                d1, d2, d3 = (
                    utils.distance(v1, v2),
                    utils.distance(v2, v3),
                    utils.distance(v3, v1),
                )
                # Check if the distances are close to each other within the defined tolerance
                if (
                    max(d1, d2, d3) - min(d1, d2, d3)
                    <= tolerance * min_distance
                ):
                    # Calculate the midpoint of the facet
                    midpoint = [
                        (v1[0] + v2[0] + v3[0]) / 3,
                        (v1[1] + v2[1] + v3[1]) / 3,
                        (v1[2] + v2[2] + v3[2]) / 3,
                    ]
                    return midpoint
        elif sym == "trunc_oct":
            for v1, v2, v3, v4 in itertools.combinations(vectors, 4):
                # Calculate distances between the vectors
                d1, d2, d3, d4 = (
                    utils.distance(v1, v2),
                    utils.distance(v2, v3),
                    utils.distance(v3, v4),
                    utils.distance(v4, v1),
                )
                # Check if the distances are close to each other within the defined tolerance
                if (
                    max(d1, d2, d3, d4) - min(d1, d2, d3, d4)
                    <= tolerance * min_distance
                ):
                    # Calculate the midpoint of the facet
                    midpoint = [
                        (v1[0] + v2[0] + v3[0] + v4[0]) / 4,
                        (v1[1] + v2[1] + v3[1] + v4[1]) / 4,
                        (v1[2] + v2[2] + v3[2] + v4[2]) / 4,
                    ]
                    return midpoint

    def displace_cages(self, arene_smile, diamine_smile, sym, metal_atom):
        dist = Cage.cage_size(self, arene_smile)
        vecs_vertex = Cage.cage_vertex_vec(
            self, diamine_smile
        )  # this has all the vectors of the arene to centroid
        vec_vertex = vecs_vertex[0]
        if sym == "4_6":
            vec_w = -Cage.calculate_cage_vector(self, arene_smile)
            vec_a = Cage.calculate_cage_vector(self, arene_smile)
        elif sym == "oct" or sym == "trunc_oct":
            vec_w = Cage.find_midpoint_of_facet(self, vecs_vertex, sym)
            vec_a = Cage.calculate_cage_vector(self, arene_smile)
        if sym == "4_6":
            guest_bb_wa = stk.BuildingBlock.init_from_molecule(self)
            guest_bb_aa = stk.BuildingBlock.init_from_molecule(self)
            origin = guest_bb_wa.get_centroid()
            guest_bb_ww = guest_bb_wa.with_rotation_between_vectors(
                vec_w, vec_w * -1, origin
            )
            guest_bb_wv = guest_bb_wa.with_rotation_between_vectors(
                vec_vertex, vec_w, origin
            )
        elif sym == "oct":
            guest_bb_ww = stk.BuildingBlock.init_from_molecule(self)
            origin = guest_bb_ww.get_centroid()
            guest_bb_wa = guest_bb_ww.with_rotation_between_vectors(
                vec_w, vec_a, origin
            )
            guest_bb_aa = guest_bb_ww
            guest_bb_wv = guest_bb_wa.with_rotation_between_vectors(
                vec_vertex, vec_w, origin
            )
        elif sym == "trunc_oct":
            guest_bb_ww = stk.BuildingBlock.init_from_molecule(self)
            origin = guest_bb_ww.get_centroid()
            guest_bb_wa = guest_bb_ww.with_rotation_between_vectors(
                vec_w, vec_a, origin
            )
            guest_bb_aa = guest_bb_ww
            guest_bb_wv = guest_bb_ww.with_rotation_between_vectors(
                vec_vertex, vec_w, origin
            )
        origin = guest_bb_wa.get_centroid()
        rotated_vectors_w = utils.generate_rotated_vectors(vec_w, 4, 30)
        rotated_vectors_a = utils.generate_rotated_vectors(vec_a, 4, 30)
        list_wa, list_ww, list_aa, list_wv = Cage.generate_new_cages(
            self,
            guest_bb_aa,
            guest_bb_ww,
            guest_bb_wa,
            guest_bb_wv,
            vec_w,
            vec_a,
            rotated_vectors_w,
            rotated_vectors_a,
            dist,
        )
        # list_wv = CageOperations.generate_ww_cages(guest_bb_wa, guest_bb_wv, vec, vec_perpendicular, rotated_vectors, dist)
        rdkit_mol = list_wa[-1][-1][-1].to_rdkit_mol()
        fixed_atom_set = Cage.fix_atom_set(
            rdkit_mol, diamine_smile, metal_atom=metal_atom
        )
        return list_wa, list_ww, list_aa, list_wv, fixed_atom_set
