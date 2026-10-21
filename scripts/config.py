import MDAnalysis as mda

mol2_file=snakemake.input.mol2
receptor_file=snakemake.input.receptor
ligand_file=snakemake.input.ligand
config_output=snakemake.output[0] #file path

def write_vina_config(mol2_file, receptor_file, ligand_file, config_path):
    #load mol2 file
    f = mda.Universe(mol2_file)

    #calculate center of geometry of atom group
    center_x, center_y, center_z = f.atoms.center_of_geometry()

    #approximate search space size
    #typically between 20 and 30 angstroms
    size_x, size_y, size_z = 20.0, 20.0, 20.0

    #write info to output file
    with open(config_path, 'w') as out:
        out.write(f'receptor = {receptor_file}\n')
        out.write(f'ligand = {ligand_file}\n')
        out.write(f'center_x = {center_x:.3f}\n')
        out.write(f'center_y = {center_y:.3f}\n')
        out.write(f'center_z = {center_z:.3f}\n')
        out.write(f'size_x = {size_x}\n')
        out.write(f'size_y = {size_y}\n')
        out.write(f'size_z = {size_z}\n')

write_vina_config(mol2_file, receptor_file, ligand_file, config_output)