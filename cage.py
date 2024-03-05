# Assuming stk and other necessary imports are already handled in this file
import stk
import rdkit.Chem.AllChem as rdkit
import numpy as np
import utils

class Cage:
    def __init__(self, file_path=None, molecule=None):
        """
        Initialize a Cage object either from a file or directly from an stk molecule.

        Args:
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
        """
        Given the SMILES strings of two building blocks,
        optimize both building blocks and construct a cage with FourPlusSix topology.
        The resulting STK molecule of the cage is centered at the origin.
        """
        bb1 = stk.BuildingBlock(
            smiles=bb1_smiles,
            functional_groups=[stk.PrimaryAminoFactory(), stk.AldehydeFactory()],
        )
        rdkit_bb1 = bb1.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_bb1)
        rdkit.MMFFOptimizeMolecule(rdkit_bb1)
        bb1 = bb1.with_position_matrix(
            position_matrix=rdkit_bb1.GetConformer().GetPositions(),
        )

        bb2 = stk.BuildingBlock(
            smiles=bb2_smiles,
            functional_groups=[stk.AldehydeFactory(), stk.PrimaryAminoFactory()],
        )
        rdkit_bb2 = bb2.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_bb2)
        rdkit.MMFFOptimizeMolecule(rdkit_bb2)
        bb2 = bb2.with_position_matrix(
            position_matrix=rdkit_bb2.GetConformer().GetPositions(),
        )

        # Verify Kekulization
        for bond in bb1.to_rdkit_mol().GetBonds():
            assert bond.GetBondTypeAsDouble() != 1.5
        for bond in bb2.to_rdkit_mol().GetBonds():
            assert bond.GetBondTypeAsDouble() != 1.5

        # Construct the cage
        self.cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusSix((bb1, bb2))
        )
        self.cage = self.cage.with_centroid([0, 0, 0])

    # Add other methods as necessary, like optimize, calculate vectors, etc.


class CageOperations:
    def __init__(self, cage):
        self.cage = cage

    def calculate_cage_vector(self, arene_smile):
        """
        Calculate the vector representation of a cage molecule.

        Parameters:
        - stk_mol: An STK molecular object representing the cage.
        - arene_smile: SMILES string for the arene part of the molecule.

        Returns:
        - A vector representation of the cage molecule.
        """

        rdkit_mol = self.to_rdkit_mol()
        rdkit.SanitizeMol(rdkit_mol)
        arenes = [self.get_centroid(atom_ids=atom_ids) for atom_ids in rdkit_mol.GetSubstructMatches(query=rdkit.MolFromSmarts(arene_smile))]
        arenes = np.asarray(arenes)
        centroid_mol = self.get_centroid()
        vectors = np.array([(arene - centroid_mol) / np.linalg.norm(arene - centroid_mol) for arene in arenes])
        return vectors[0]

    def generate_new_cages(self, guest_bb_aa, guest_bb_ww, guest_bb_wa, direction_vector, perpendicular_vector, rotated_vectors, displacement_distance):
        """
        Generate new cage complexes with various packings and rotations.

        Parameters:
        - cage: The base cage molecule.
        - guest_bb_aa: Guest molecule with arene-to-arene packing.
        - guest_bb_ww: Guest molecule with window-to-window packing.
        - guest_bb_wa: Guest molecule with window-to-arene packing.
        - direction_vector: Direction vector for displacement.
        - perpendicular_vector: Vector perpendicular to the direction vector.
        - rotated_vectors: A list of vectors representing all the rotation directions.
        - displacement_distance: The distance for displacement.

        Returns:
        - Lists of constructed complexes with different packings and rotations.
        - A list of centroids for all displaced cages.
        """
        list_wa, list_ww, list_aa, list_centroids = [], [], [], []
        for i in np.arange(0, 7, 1):
            list_wa_run,list_ww_run,list_aa_run  = [],[],[]
            displaced_centers =utils.find_integer_points(direction_vector, displacement_distance*2-2+i, int(displacement_distance) + 1)

            for center in displaced_centers:
                for rotation_vector in rotated_vectors:
                    list_centroids.append(center)
                    guest_aa = utils.create_rotated_guest(guest_bb_aa, perpendicular_vector, rotation_vector, center)
                    guest_wa = utils.create_rotated_guest(guest_bb_wa, perpendicular_vector, rotation_vector, center)
                    guest_ww = utils.create_rotated_guest(guest_bb_ww, perpendicular_vector, rotation_vector, center * -1)

                    complex_wa = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(host=self, guests=guest_wa)
                    )
                    complex_ww = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(host=self, guests=guest_ww)
                    )
                    complex_aa = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(host=self, guests=guest_aa)
                    )

                    list_wa_run.append(complex_wa)
                    list_ww_run.append(complex_ww)
                    list_aa_run.append(complex_aa)
            list_wa.append(list_wa_run)
            list_ww.append(list_ww_run)
            list_aa.append(list_aa_run)
        return list_wa, list_ww, list_aa, list_centroids
    
    def generate_ww_cages(cage,guest_bb_ww, vec, vec_perpendicular, rotated_vectors, dist):  
        """
        Args:
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

        list_ww = []
        list_centroids = []
        for i in np.arange(0, 7, 1):
            list_ww_run = []
            displaced_centres=utils.find_integer_points(vec, dist*2-2+i,int(dist)+1)

            for k in range(len(displaced_centres)):
                for j in range(len(rotated_vectors)):
                    vec_new = displaced_centres[k]
                    #vec_new = displace_vector(vec, dist)
                    list_centroids.append(vec_new)
                    #print((rotated_vectors[j]),k*len(rotated_vectors)+j,i)
                    guest_ww = utils.create_rotated_guest(guest_bb_ww, vec_perpendicular, rotated_vectors[j], vec_new*-1)
                    complex_ww = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(
                            host=cage,
                            guests=guest_ww,
                        ),
                    )
                    list_ww_run.append(complex_ww)
            list_ww.append(list_ww_run)
        return  list_ww, list_centroids