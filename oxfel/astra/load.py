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
    assert not are_parallel(v_from, v_to)
    assert not are_anti_parallel(v_from, v_to)

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
    new_positions = np.einsum("ij,hj", rot, beam[["x", "y", "z"]] - average_position)

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


def RRv(u1, u2):
    # % x-y-x
    th1 = np.arctan2(u1[1],u1[2]);
    M1x = RRx(-th1);
    th2 = np.arctan2(u2[1],u2[2]);
    M2x = RRx(-th2);
    v1 = M1x @ u1;
    v2 = M2x @ u2;
    b = np.arctan2(v2[2],v2[0]);
    a = np.arctan2(v1[2],v1[0]);
    My = RRy(b-a);
    R = M2x.T @ My @ M1x;
    return R

def RRx(theta):
    # % ccw rotation with rotation axis out of the page
    R = np.array([[1,  0,           0],
                  [0,  np.cos(theta), -np.sin(theta)],
                  [0,  np.sin(theta),  np.cos(theta)]])
    return R


def RRy(theta):
    # % ccw rotation with rotation axis out of the page
    R = np.array([[np.cos(theta), 0,  np.sin(theta)],
                  [0,           1,  0],
                  [-np.sin(theta), 0,  np.cos(theta)]])
    return R

# % % read in phase space
# % phasespace = dlmread('E:\EG\Work\DESY_2023\XFEL_IBS2\Simulations\Astra_spch_withdipoles\euxfel-injector.2451.001');
# phasespace = dlmread('E:\EG\Work\DESY_2023\XFEL_IBS2\Simulations\Reptil_spch\beam-42.5.astra');
# phasespace = phasespace(:,1:10);

# % translations and scalings
# phasespace(2:end,3) = phasespace(1,3) + phasespace(2:end,3); % z
# phasespace(2:end,6) = phasespace(1,6) + phasespace(2:end,6); % pz
# %phasespace(2:end,7) = phasespace(1,7) + phasespace(2:end,7); % clock
# phasespace(:,8) = phasespace(:,8)*1E-9; % q
# % phasespace(:,3) = phasespace(:,3) - mean(phasespace(:,3));
# % phasespace(:,7) = phasespace(:,7) - mean(phasespace(:,7));
              
# % rotate beam
# r0 = [mean(phasespace(:,1)) mean(phasespace(:,2)) mean(phasespace(:,3))]';
# u0 = [mean(phasespace(:,4)) mean(phasespace(:,5)) mean(phasespace(:,6))]';

# if (Rotate_Beam)
#     nz = [0 0 1]';
#     R = RRv(nz,u0);
#     for i=1:Np
#         r1 = [phasespace(i,1) phasespace(i,2) phasespace(i,3)]';
#         u1 = [phasespace(i,4) phasespace(i,5) phasespace(i,6)]';
#         r2 = R*r1; %R*(r1-r0);
#         u2 = R*u1; %R*(u1-u0);
#         phasespace(i,1) = r2(1);
#         phasespace(i,2) = r2(2);
#         phasespace(i,3) = r2(3);
#         phasespace(i,4) = u2(1);
#         phasespace(i,5) = u2(2);
#         phasespace(i,6) = u2(3);
#     end
# end

# % test rotation
# rr = [mean(phasespace(:,1)) mean(phasespace(:,2)) mean(phasespace(:,3))]';
# ur = [mean(phasespace(:,4)) mean(phasespace(:,5)) mean(phasespace(:,6))]';

# % parameters
# param.N_slices=300;
# param.Np_plot=size(phasespace(:,1),1);
# param.L_no_zcorr=true;
# param.X_hor=1;
# param.X_ver=2;
# param.X_long=3;
# param.X_phase=4;
# param.Slice_Analysis=5;
# param.X_space = 6;
# param.Np_plot = 100000;

# out_postpro = postpro(phasespace(:,[1 4 2 5 3 6 8]),param);

# %-eg dlmwrite( 'E:\EG\Work\DESY_2023\PSI2\Reptil\Ibs\Ref_200k\beam-66.txt', out_postpro.slice, 'delimiter','\t','precision','%e');

# % Plot slice emittance
#  %-egfigure(3);
#  %-egplot(out_postpro.slice(:,1),out_postpro.slice(:,3))

# % Plot slice energy spread
# %-egfigure(4);
# %-egplot(out_postpro.slice(:,1),out_postpro.slice(:,6))
    
    
if __name__ == '__main__':
    main(sys.argv[1])
