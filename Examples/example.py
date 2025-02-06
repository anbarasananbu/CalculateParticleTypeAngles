input_file = "lammps.data"
pipeline = import_file(input_file)

pipeline.modifiers.append(CalculateParticleTypeAngles(donor_type= 5, accep_type=5,  center_type=6, DA_rcut=3.5, DHA_acut=150, HD_rcut=1.05, HA_rcut=3.0))

export_file(pipeline, "21.test/test.dat", format="txt/table", key="angle-triplet")
data = pipeline.compute()
print("total_angles", data.attributes["No-angles"])