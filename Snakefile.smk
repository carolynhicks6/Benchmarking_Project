IDS = ['1bzc', '1c5z']

rule all:
    input:
        expand("sample_data/{id}/{id}_proteinH.pdb", id=IDS), 
        expand("sample_data/{id}/{id}_protein.pqr", id=IDS),
        expand("sample_data/{id}/{id}_protein.pdbqt", id=IDS),
        expand("sample_data/{id}/{id}_ligand.pdbqt", id=IDS),
        expand("sample_data/{id}/{id}_config.txt", id=IDS),
        expand("sample_data/{id}/{id}_vina_out.pdbqt", id=IDS),
        expand("sample_data/{id}/{id}_proteinH.mol2", id=IDS),
        expand("sample_data/{id}/{id}_ligand_fixed.sdf", id=IDS),
        expand("sample_data/{id}/{id}_protein.prm", id=IDS),
        expand("sample_data/{id}/{id}_protein.as", id=IDS),
        expand("sample_data/{id}/{id}_rdock_out.sd", id=IDS)

#remove unwanted water molecules, add heavy atoms and hydrogens where needed to new pdb file
#pqr file contains additional info about atomic radii and atomic charges
rule fix_pdb:
    input:
        "sample_data/{id}/{id}_protein.pdb"
    output:
        pdb = "sample_data/{id}/{id}_proteinH.pdb",
        pqr = "sample_data/{id}/{id}_protein.pqr"
    shell:
        """
        pdb2pqr --pdb-output={output.pdb} --pH=7.4 {input} {output.pqr} --whitespace
        """

#output of previous rule to pdbqt
# -xr preserve residues and treat as rigid
# -r retain hydrogens
rule protein_pdbqt:
    input:
        "sample_data/{id}/{id}_protein.pqr"
    output:
        "sample_data/{id}/{id}_protein.pdbqt"
    shell:
        """
        obabel {input} -O {output} -xr -r
        """

#ligand preparation using meeko - specialized for small molecules
rule ligand_pdbqt:
    input:
        "sample_data/{id}/{id}_ligand_opt.mol2"
    output:
        "sample_data/{id}/{id}_ligand.pdbqt"
    shell:
        """
        mk_prepare_ligand.py -i {input} -o {output}
        """

#parse coordinates of crystal ligand from optimized mol2 and write config.txt
#run with --use-conda
rule config:
    input:
        mol2 = "sample_data/{id}/{id}_ligand_opt.mol2",
        receptor = "sample_data/{id}/{id}_protein.pdbqt",
        ligand = "sample_data/{id}/{id}_ligand.pdbqt"
    output:
        "sample_data/{id}/{id}_config.txt"
    conda:
        "docking_pipeline"
    script:
        "scripts/config.py"

#run vina using config file
rule run_vina:
    input:
        "sample_data/{id}/{id}_config.txt"
    output:
        "sample_data/{id}/{id}_vina_out.pdbqt"
    shell:
        """
        vina --config {input} --out {output}
        """

#prepare receptor - add hydrogens
rule fix_mol2:
    input:
        "sample_data/{id}/{id}_protein.mol2"
    output:
        "sample_data/{id}/{id}_proteinH.mol2"
    shell:
        """
        obabel {input} -O {output} -h
        """       

#next, fix ligand by making 3D and adding hydrogens
rule fix_sdf:
    input:
        "sample_data/{id}/{id}_ligand.sdf"
    output:
        "sample_data/{id}/{id}_ligand_fixed.sdf"
    shell:
        """
        obabel {input} -O {output} --gen3d -p 7.4 --rebuild-bonds --partialcharge gasteiger
        """

#another cleanup option
#mol2 to sdf - often cleaner
rule fix:
    input:
        mol="sample_data/{id}/{id}_ligand.mol2"
    output:
        fixed="sample_data/{id}/{id}_ligand_fixed.sdf"
    shell:
        """
        obabel {input.mol} -O {output.fixed} -h
        """
    
#first step of rdock is system definition - write .prm file 
rule write_prm:
    input:
        receptor="sample_data/{id}/{id}_proteinH.mol2",
        ligand="sample_data/{id}/{id}_ligand_fixed.sdf"
    output:
        "sample_data/{id}/{id}_protein.prm"
    shell:
        """
        cat << 'EOF' > {output} 
RBT_PARAMETER_FILE_V1.00
TITLE {wildcards.id}

RECEPTOR_FILE {input.receptor}
RECEPTOR_FLEX 3.0

SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {input.ligand}
    RADIUS 6.0
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
    VOL_INCR 0.0
    GRIDSTEP 0.5
END_SECTION

SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION
EOF
        """

#next rdock requires .as file which defines region of receptor that ligand should bind into
#this is called cavity mapping
rule cavity_mapping:
    input:
        "sample_data/{id}/{id}_protein.prm"
    output:
        "sample_data/{id}/{id}_protein.as"
    shell:
        """
        rbcavity -r {input} -W
        """

#last step is to actually run rdock
#output is a .sd file
rule run_rdock:
    input:
        ligand="sample_data/{id}/{id}_ligand_fixed.sdf",
        receptor="sample_data/{id}/{id}_protein.prm"
    output:
        docked="sample_data/{id}/{id}_rdock_out.sd"
    shell:
        """
        OUT_PREFIX="{output.docked}"
        rbdock -i {input.ligand} -o "$${{OUT_PREFIX%.sd}}" -r {input.receptor} -p dock.prm
        """