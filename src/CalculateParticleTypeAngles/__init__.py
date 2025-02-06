from ovito.data import CutoffNeighborFinder, DataCollection, DataTable
from ovito.pipeline import ModifierInterface
from traits.api import Int, Float
import numpy as np 

class CalculateParticleTypeAngles(ModifierInterface):
    donor_type = Int()
    accep_type = Int()
    center_type = Int()
    HD_rcut = Float()
    HA_rcut = Float()
    DA_rcut = Float()
    DHA_acut = Float()
    

    def __init__(self, donor_type, accep_type, center_type, HD_rcut, HA_rcut, DA_rcut, DHA_acut):
        self.donor_type = donor_type
        self.accep_type = accep_type
        self.center_type = center_type
        self.HD_rcut = HD_rcut
        self.HA_rcut = HA_rcut
        self.DA_rcut = DA_rcut
        self.DHA_acut = DHA_acut
    
    # Function to calculate angle (donor-hydrogen-acceptor) with PBC handling
    def Calculate_angle(self, vec_hd, vec_ha):
        """Calculate the angle (in degrees) between two vectors."""
        cos_ang = np.dot(vec_hd, vec_ha) / (np.linalg.norm(vec_hd) * np.linalg.norm(vec_ha))
        return np.round(np.degrees(np.arccos(np.clip(cos_ang, -1.0, 1.0))), 3)
    
    def modify(self, data, **kwargs):
        """Identify hydrogen bonds based on distance and angle criteria."""
        # Initialize neighbor finders
        donor_finder = CutoffNeighborFinder(cutoff=self.HD_rcut, data_collection=data)
        acceptor_finder = CutoffNeighborFinder(cutoff=self.HA_rcut, data_collection=data)

        # Get hydrogen atom indices (type 6)
        ptypes = data.particles.particle_types
        center_indices = np.where(ptypes == self.center_type)[0]

        donor_list, center_list, acceptor_list, angles_list = [], [], [], []

        # Process each center atom
        for H_idx in center_indices:
            donor_data = next((neigh for neigh in donor_finder.find(H_idx) if ptypes[neigh.index] == self.donor_type), None)
            if not donor_data:
                continue  # Skip if no donor found
    
            vec_hd, donor_idx, dist_hd = donor_data.delta, donor_data.index, donor_data.distance

            # Find all valid acceptors
            for A_neigh in acceptor_finder.find(H_idx):
                if A_neigh.distance > self.HD_rcut and ptypes[A_neigh.index] == self.accep_type:
                    vec_ha, acceptor_idx, dist_ha = A_neigh.delta, A_neigh.index, A_neigh.distance
                    
                    angle = self.Calculate_angle(vec_hd, vec_ha)
                    dist_da = np.sqrt(dist_hd**2 + dist_ha**2 - 2 * dist_hd * dist_ha * np.cos(np.radians(angle)))
    
                    if angle >= self.DHA_acut and dist_da <= self.DA_rcut:
                        angles_list.append(angle)
                        donor_list.append(donor_idx)
                        center_list.append(H_idx)
                        acceptor_list.append(acceptor_idx)
    
        # Convert index to particle Identifier
        donor_list = np.array(data.particles.identifiers[donor_list])
        center_list = np.array(data.particles.identifiers[center_list])
        acceptor_list = np.array(data.particles.identifiers[acceptor_list])

        # Store results in a table
        table = data.tables.create(identifier="angle-triplet", title="Bond Angles of Particle", plot_mode=DataTable.PlotMode.NoPlot)
        table.y = table.create_property("Angle", data=angles_list)
        table.x = table.create_property("Particle", data=np.column_stack([donor_list, center_list, acceptor_list]), components=["A", "B", "C"])
        data.attributes["No-angles"] = len(angles_list)
