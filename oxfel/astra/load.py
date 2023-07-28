from pathlib import Path
import pandas as pd
import numpy as np
import sys


COLUMNS = ["x", "y", "z", "px", "py", "pz", "clock", "charge", "index", "flag"]

def are_parallel(vector_1, vector_2, tolerance=1e-10):
    """
    Check if vector vector_1 is parallel to vector vector_2 down to
    some tolerance.

    :param vector_1: First input vector
    :type vector_1: array
    :param vector_2: Second input vector
    :type vector_2: array
    :param tolerance: Tolerance for calculation
    :type tolerance: float
    :returns: if vectors are parallel
    :rtype: bool

    """
    return (np.dot(vector_1, vector_2) / (np.linalg.norm(vector_1)
                                           * np.linalg.norm(vector_2))
            > 1 - tolerance)

def are_anti_parallel(vector_1, vector_2, tolerance=1e-10):
    """
    Check if vector vector_1 is parallel to vector vector_2 down to
    some tolerance.

    :param vector_1: First input vector
    :type vector_1: array
    :param vector_2: Second input vector
    :type vector_2: array
    :param tolerance: Tolerance for calculation
    :type tolerance: float
    :returns: if vectors are antiparallel
    :rtype: bool
    """
    return (np.dot(vector_1, vector_2) / (np.linalg.norm(vector_1)
                                          * np.linalg.norm(vector_2))
            < -1 + tolerance)

def matrix_from(v_from, v_to):
    """
    Returns the rotation matrix that rotates v_from to parallel to
    v_to.

    v_from and v_to should be array-like three-vectors.

    """
    if are_parallel(v_from, v_to):
        return np.eye(3)
    elif are_anti_parallel(v_from, v_to):
        return _rodrigues_anti_parallel(v_from, v_to)

    axis = (np.cross(v_from,
                      v_to)
            / np.linalg.norm(np.cross(v_from,
                                        v_to)))
    angle = np.arccos(np.dot(v_from,
                               v_to)
                       / (np.linalg.norm(v_from)
                          * np.linalg.norm(v_to)))

    # Construct the skew-symmetric cross product matrix.
    cross_matrix = np.array([[0,       -axis[2],  axis[1]],
                             [axis[2],       0,  -axis[0]],
                             [-axis[1], axis[0],        0]])
    # Rodrigues' rotation formula.
    rotation_matrix = (np.eye(3)
                       + (np.sin(angle) * (cross_matrix))
                       + ((1 - np.cos(angle))
                          * (cross_matrix).dot(cross_matrix)))

    assert are_parallel(v_to, rotation_matrix.dot(v_from)), (
        "not parallel")
    assert np.allclose(rotation_matrix.T.dot(rotation_matrix), np.eye(3)), (
        "Rotation matrix not orthogonal")
    return rotation_matrix


def new_load(fname, outdir):
    df = pd.read_csv(fname, delim_whitespace=True, names=COLUMNS)

    # Get the reference particle
    reference_particle = df.iloc[0]
    # And the rest of the beam sans reference particle
    beam = df.iloc[1:]

    # Update coords that are relative to the reference particle
    beam.loc[:, "z"] += reference_particle.z
    beam.loc[:, "pz"] += reference_particle.pz
    beam.loc[:, "clock"] += reference_particle.clock

    # Average momentum and position of the beam
    average_momentum = beam[["px", "py", "pz"]].mean()
    average_position = beam[["x", "y", "z"]].mean()

    # Get the rotation matrix that maps the average momentum onto the z-axis
    rot = matrix_from(average_momentum , [0, 0, 1])

    # Rotate the momenta
    new_momenta = np.einsum("ij,hj", rot, beam[["px", "py", "pz"]])
    # Rotate the positions after centering them.
    new_positions = np.einsum("ij,hj", rot, beam[["x", "y", "z"]])

    # Set the new positions and momenta in the dataframe
    beam.loc[:, ["x", "y", "z"]] = new_positions
    beam.loc[:, ["px", "py", "pz"]] = new_momenta

    # Now update the reference particle so that it is at the mean
    # position, momenta and clock
    reference_particle = beam.iloc[0].copy()
    reference_particle[["x", "y", "z"]] = beam[["x", "y", "z"]].mean()
    reference_particle[["px", "py", "pz"]] = beam[["px", "py", "pz"]].mean()
    reference_particle["clock"] = beam["clock"].mean()

    # Now subtract the reference particle relative longitudinal coordinates
    beam.loc[:, ["z", "pz", "clock"]] -= reference_particle[["z", "pz", "clock"]]

    # Put the reference particle back again at the head of the
    # dataframe and sort by index so that it's at the top
    beam.loc[0] = reference_particle
    beam = beam.sort_index()

    # Write.
    outdir = Path(outdir)
    np.savetxt(outdir / fname, beam)


def main(fname):
    new_load(fname, "test")


if __name__ == '__main__':
    main(sys.argv[1])
