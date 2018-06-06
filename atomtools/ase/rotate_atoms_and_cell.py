import numpy as np

def rotate_atoms_and_cell( atoms=None, target_direction=(0,0,1) ):
    """
    This function rotates an atoms object such that third axis of the unitcell
    is parallel to the target_z_direction
    """
    cell = atoms.cell.T

    # Normalize the vectors
    orig_z = cell[:,2]
    orig_z /= np.sqrt( np.sum(orig_z**2) )
    target_direction = np.array(target_direction).astype(np.float64)
    target_direction /= np.sqrt( np.sum(target_direction**2) )

    if ( np.allclose(orig_z,target_direction) ):
        return

    rot_axis = np.cross(target_direction,orig_z)
    rot_axis /= np.sqrt( np.sum(rot_axis**2) )
    angle = np.arccos( orig_z.dot(target_direction) )
    atoms.rotate( -angle*180.0/np.pi, v=rot_axis, rotate_cell=True )

    # Make sure that the rotation worked as it should
    cell = atoms.cell.T
    new_z_direction = atoms.cell.T[:,2]
    new_z_direction /= np.sqrt(np.sum(new_z_direction**2))

    dotprod = new_z_direction.dot(target_direction)
    if ( np.abs(dotprod-1.0) > 1E-6 ):
        angle = np.arccos(dotprod)*180.0/np.pi
        raise RuntimeError( "Rotation failed. Final angle between target_z_direction and resulting {} deg".format(angle) )
    return atoms

def rotate_crd_system( atoms=None, target_z_direction=(0,0,1) ):
    """
    Rotate the coordinate system of the atoms object such that its z-axis
    points along target_z_direction.
    """

    target_z_direction = np.array( target_z_direction ).astype(np.float64)
    target_z_direction /= np.sqrt( np.sum(target_z_direction**2) )
    orig_z_direction = np.array([0,0,1])

    if ( np.allclose(orig_z_direction,target_z_direction) ):
        return

    rot_axis = np.cross(target_z_direction,orig_z_direction)
    rot_axis /= np.sqrt( np.sum(rot_axis**2) )
    angle = np.arccos( orig_z_direction.dot(target_z_direction) )
    atoms.rotate( -angle*180.0/np.pi, v=rot_axis, rotate_cell=True )
    return atoms

def align_direction_with_z( atoms=None, direction=(0,0,1) ):
    """
    Rotates the atoms cell such that the direction specified is oriented along the z axis
    """
    cell = atoms.get_cell().T
    dir_to_align = direction[0]*cell[:,0] + direction[1]*cell[:,1] + direction[2]*cell[:,2]
    dir_to_align /= np.sqrt( np.sum(dir_to_align**2) )

    zhat = np.array( [0,0,1] )
    rot_axis = np.cross( zhat, dir_to_align )
    rot_axis /= np.sqrt( np.sum(rot_axis**2) )
    angle = np.arccos( zhat.dot(rot_axis) )
    atoms.rotate( -angle*180.0/np.pi, v=rot_axis, rotate_cell=True )
    return atoms
