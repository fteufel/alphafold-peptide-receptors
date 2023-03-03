import subprocess
import tempfile
import os

ROSETTA_PATH= '/isdata/winthergrp/zpf738/rosetta_bin_linux_2021.16.61629_bundle/main/source/bin/'
DATABASE = '/isdata/winthergrp/zpf738/rosetta_src_2021.16.61629_bundle/main/database'
SCORING_XML = '2_score_filters.xml'
ROSETTA_RELAXED_DIR = '/isdata/winthergrp/zpf738/af/predictions_pdb_benchmark_relaxed'

ROSETTA_LOG_DIR = '/isdata/winthergrp/zpf738/af/rosetta_logs'

def score_complex(
        pdb_file: str, 
        peptide_chain = 'C',
        receptor_chain='B',
        rosetta_path = None,
        database_path = None,
        scoring_xml=None,
        relaxed_dir=None,
        log_dir=None
        ) -> dict:
    '''
    Run rosetta to relax and score the complex.
    '''


    if rosetta_path is None:
        rosetta_path = ROSETTA_PATH
    if database_path is None:
        database_path = DATABASE
    if scoring_xml is None:
        scoring_xml = SCORING_XML
    if relaxed_dir is None:
        relaxed_dir = ROSETTA_RELAXED_DIR
    if log_dir is None:
        log_dir = ROSETTA_LOG_DIR

    out_dir = os.path.join(relaxed_dir, os.path.normpath(pdb_file).split(os.sep)[-2])
    os.makedirs(out_dir, exist_ok=True)
    out_pdb_name = os.path.basename(pdb_file).split('.')[0] + '_0001.pdb'
        
    # 1. prepare the structure
    if not os.path.exists(os.path.join(out_dir, out_pdb_name)):
        command = [
            os.path.join(rosetta_path, 'relax.static.linuxgccrelease'),
            '-database', database_path,
            '-relax:constrain_relax_to_start_coords',
            '-relax:coord_constrain_sidechains',
            '-relax:ramp_constraints', 'false',
            '-ignore_zero_occupancy', 'false',
            '-ex1',
            '-ex2',
            '-use_input_sc',
            '-flip_HNQ',
            '-no_optH', 'false',
            '-ignore_waters', 'false',
            '-auto_setup_metals', 'true',
            '-score:weights', 'ref2015',
            '-in:file:s', pdb_file,
            '-out:path:all', out_dir,
        
        ]


        subprocess.run(command, stdout=open(os.path.join(log_dir, pdb_file.replace('/', '_')+'_relax.log'), 'w'))
        # ranked_0.pdb --> out_dir/ranked_0_0001.pdb

    with tempfile.TemporaryDirectory() as tmpdir:

        # 2. score
        command = [
            os.path.join(rosetta_path, 'rosetta_scripts.static.linuxgccrelease'),
            '-database', database_path,
            '-in:file:s', os.path.join(out_dir, out_pdb_name),
            '-parser:protocol', scoring_xml,
            '-parser:script_vars', f'ReceptorChain={receptor_chain}', f'PeptideChain={peptide_chain}',
            '-out:path:all', tmpdir,
            '-out:user_tag',  pdb_file,
            '-out:file:silent', 'ScoredComplexes.out'

        ]

        subprocess.run(command, stdout=open(os.path.join(log_dir, pdb_file.replace('/', '_')+'_scoring.log'), 'w'))

        with open(os.path.join(tmpdir, 'ScoredComplexes.out'), 'r') as f:
            lines = f.readlines()

            keys = lines[1].split(':')[1].split()
            values = lines[3].split(':')[1].split()

            result = dict(zip(keys, values))


        return result

        #SEQUENCE: MPLVVAVIFFSLWVFALGQLEQPEISISRPANKSAHISWKASIQGFSSKIIHWYWQKPNKGLEYLLHVFLTISAQDCSGGKTKKLEVSKNAHTSTSTLKIKFLEKEDEVVYHCACWIRHCYIQNCPLG
        #SCORE:     score     fa_atr     fa_rep     fa_sol    fa_intra_rep    fa_intra_sol_xover4    lk_ball_wtd    fa_elec    pro_close    hbond_sr_bb    hbond_lr_bb    hbond_bb_sc    hbond_sc    dslf_fa13    atom_pair_constraint    coordinate_constraint    angle_constraint    dihedral_constraint      omega     fa_dun    p_aa_pp    yhh_planarity        ref    rama_prepro    contact_MS    ddg    ddg_norepack    interface_buried_sasa    interface_hydrophobic_sasa    interface_polar_sasa    peptide_score    shape_complementarity    total_energy    shape_complementarity_int_area       time                                                          user_tag        description
        #REMARK BINARY SILENTFILE
        #SCORE:   -21.839   -612.932    144.809    351.125           1.438                 18.609        -13.349   -160.880        0.698        -10.666        -35.542         -7.339      -5.173        0.000                   0.000                    0.000               0.000                  0.000     35.076    197.487    -18.436            0.399     66.460         26.376       211.604 -13.422         -14.135                 1175.085                       814.151                 360.934           18.710                    0.488         -21.839                           591.617      9.000 ../predicti ranked_0_0001_0001

