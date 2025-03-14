import h5py
import awkward as ak
import numpy as np
from pathlib import Path


def read_depletion_results(filename: str | Path = 'depletion_results.h5') -> ak.Array:
    """
    Reads depletion results from an OpenMC depletion_results.h5 HDF5 file and
    returns a memory efficient awkward array of atom quantities.

    Args:
        filename: Path to the HDF5 file containing depletion results.

    Returns:
        awkward.Array: An awkward array containing atom quantities for each timestep.
            - "timestep" (float): The end time of the timestep in seconds.
            - "material_id" (int): The ID of the material.
            - "nuclide" (str): The name of the nuclide.
            - "atoms" (float): The quantity of the nuclides in the material at that timestep in atoms.
    """

    all_atom_quantities = ak.Array([])

    with h5py.File(str(filename), "r") as handle:

        mat_id_map = {}
        nuclide_map = {}

        for mat, mat_handle in handle["/materials"].items():
            vol = mat_handle.attrs["volume"]
            ind = mat_handle.attrs["index"]

            mat_id_map[int(ind)] = int(mat)

        for nuc, nuc_handle in handle["/nuclides"].items():
            ind_atom = nuc_handle.attrs["atom number index"]
            nuclide_map[int(ind_atom)] = nuc

        # Get number of results stored
        n = handle["number"][...].shape[0]

        for step in range(n):

            number_dset = handle["/number"]
            data = number_dset[step, :, :, :]

            timestep_data = data[0, :, :]
            materials, nuclides = np.nonzero(timestep_data)
            atoms = timestep_data[materials, nuclides]
            # time returns a tuple with start and end time, end time is the one which corresponds to the atoms

            timestep_awk = ak.Array(
                [
                    {
                        "timestep": step,
                        "material_id": int(mat_id_map[m]),
                        "nuclide": nuclide_map[n],
                        "atoms": a,
                    }
                    for m, n, a in zip(materials, nuclides, atoms)
                ]
            )

            all_atom_quantities = ak.concatenate([all_atom_quantities, timestep_awk])

    return all_atom_quantities
