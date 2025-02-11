"""Module for generating dimers."""

import numpy as np
import stk

from dimer_calculations import utils


class DimerGenerator:
    """Class to generate dimers."""

    def __init__(
        self,
        displacement: float = 7,
        displacement_step_size: float = 1,
        rotation_limit: float = 120,
        rotation_step_size: float = 30,
        overlap_tolerance: float = 0.1,
        slide: bool = False,
        radius: float = 1,
    ):
        self._displacement = displacement
        self._displacement_step_size = displacement_step_size
        self._rotation_limit = rotation_limit
        self._rotation_step_size = rotation_step_size
        self._overlap_tolerance = overlap_tolerance
        self._slide = slide
        self._radius = radius

    def generate(
        self,
        molecule: stk.Molecule,
        molecule_2: stk.Molecule,
        axes_1: np.ndarray,
        axes_2: np.ndarray,
        displacement_distance: float,
    ):
        origin_2 = molecule_2.get_centroid()
        guest_cage = molecule_2.with_rotation_between_vectors(
            axes_2, -axes_1, origin_2
        )
        rotated_vectors = utils.generate_rotated_vectors(
            axes_1, self._rotation_limit / self._rotation_step_size, 30
        )
        perpendicular_vector = utils.calculate_perpendicular_vector(axes_1)

        dimer_list = []
        for i in range(int(self._displacement / self._displacement_step_size)):
            if self._slide:
                displaced_centers = utils.find_integer_points(
                    axes_1, displacement_distance + i, self._radius + 1
                )
            else:
                displaced_centers = [(displacement_distance + i) * axes_1]
            slide_up = 0
            for center in displaced_centers:
                rot_by = 0
                for vector in rotated_vectors:
                    rotated_guest = utils.create_rotated_guest(
                        guest_cage, perpendicular_vector, vector, center
                    )
                    dimer = stk.ConstructedMolecule(
                        topology_graph=stk.host_guest.Complex(
                            host=molecule,
                            guests=rotated_guest,
                        )
                    )
                    mol = dimer.to_rdkit_mol()
                    overlaps = utils.check_overlaps(
                        mol, self._overlap_tolerance
                    )
                    if overlaps:
                        rot_by = rot_by + 1
                        continue
                    dimer_list.append(
                        {
                            "Displacement shell": i,
                            "Slide": slide_up,
                            "Rotation": rot_by * self._rotation_step_size,
                            "Displacement centroid": center,
                            "Dimer": dimer,
                        }
                    )
                    rot_by = rot_by + 1
                slide_up = slide_up + 1
        return dimer_list
