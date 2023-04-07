# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2022-12-29 23:04:12
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-04-05 17:34:22

from noise2read.config import Config
import sys, getopt
from noise2read.data_generation import DataGneration
from noise2read.error_orrection import ErrorCorrection
from noise2read.data_analysis import DataAnalysis
import os
from noise2read.simulation import Simulation
from noise2read.utils import custom_logger, usage
from noise2read.data_preprocessing import DataProcessing
from contextlib import redirect_stdout

def main():
    argv = sys.argv[1:]
    # create logger
    logger = custom_logger("noise2read", debug_mode=False)

    ##############################################################
    try:
        # opts, args = getopt.getopt(argv, "m:c:i:u:t:r:d:p:a:g:o:l:h:", ["module=", "config=", "input=", "umi_file", "true", "rectification", "directory", "parallel", "high_ambiguous", "tree_method", "over_sampling", "libray_layout", "help"]) 
        opts, args = getopt.getopt(argv, "m:c:i:u:t:r:d:p:a:g:l:h", ["module=", "config=", "input=", "umi_file", "true", "rectification", "directory", "parallel", "high_ambiguous", "tree_method", "libray_layout", "help"]) 
        script_name = sys.argv[0]
        input_commands = [script_name]

        for opt, arg in opts:
            input_commands.append(opt)
            if arg:
                input_commands.append(arg)

        for arg in args:
            input_commands.append(arg)
        logger.info("Commands: " + " ".join(input_commands))
        
        if opts:
            opts_dict = dict(opts)
            opts_keys = list(opts_dict.keys())
            tar_set = set(opts_keys)

            h_lst = list({'-h', '--help'}.intersection(tar_set))
            if h_lst:
                usage()           
                sys.exit()  
            else: 
                m_lst = list({"-m", "--module"}.intersection(tar_set))
                c_lst = list({"-c", "--config"}.intersection(tar_set)) 
                i_lst = list({"-i", "--input"}.intersection(tar_set)) 
                u_lst = list({"-u", "--umi_file"}.intersection(tar_set)) 
                t_lst = list({"-t", "--true"}.intersection(tar_set)) 
                r_lst = list({"-r", "--rectification"}.intersection(tar_set))   
                d_lst = list({"-d", "--directory"}.intersection(tar_set))
                p_lst = list({"-p", "--parallel"}.intersection(tar_set))
                a_lst = list({"-a", "--high_ambiguous"}.intersection(tar_set))
                g_lst = list({"-g", "--tree_method"}.intersection(tar_set))
                # o_lst = list({"-o", "--over_sampling"}.intersection(tar_set))

                # l_lst = list({"-l", "--libray_layout"}.intersection(tar_set))
                # ref_lst = list({"-f", "--reference_in"}.intersection(tar_set))
                # r0_lst = list({"-0", "--read"}.intersection(tar_set))
                # r1_lst = list({"-1", "--read1"}.intersection(tar_set))
                # r2_lst = list({"-2", "--read1"}.intersection(tar_set))
                # align_lst = list({"-A", "--Alignment"}.intersection(tar_set))
                    
                if m_lst:
                    module_arg = opts_dict[m_lst[0]]
                        
                else:
                    logger.error("Not select any module, using command noise2read -h/--help for usage.")
                    sys.exit()
                available_cpu_cores = os.cpu_count()

############################################################################################################################
                if module_arg == "correction":
                    # try: 
                    if c_lst:
                        config = Config(opts_dict[c_lst[0]], logger) 
                        if i_lst:
                            config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Must set sequencing dataset in configuration.')
                            sys.exit()
                        if t_lst:
                            config.ground_truth_data = opts_dict[t_lst[0]]
                            if not os.path.exists(config.ground_truth_data):
                                logger.error('Ground truth data does not exsit.')
                    elif i_lst:
                        config = Config(None, logger)  
                        config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Input file does not exsit.')
                            sys.exit()
                    else:
                        logger.error('Must input configuration file or sequencing dataset.')
                        sys.exit()
                    ##############################################################
                    if t_lst:
                        config.ground_truth_data = opts_dict[t_lst[0]] 
                        if not os.path.exists(config.ground_truth_data):
                            logger.exception('Ground truth data does not exsit.')
                            raise
                    if d_lst:
                        config.result_dir = opts_dict[d_lst[0]] 
                    if p_lst:
                        config.num_workers = int(opts_dict[p_lst[0]])      
                    if a_lst:
                        config.high_ambiguous = opts_dict[a_lst[0]]
                    if g_lst:
                        config.tree_method = opts_dict[g_lst[0]]
                    # if o_lst:
                    #     if opts_dict[o_lst[0]] == 'False' or 'false':
                    #         config.over_sampling = False
                    #     elif opts_dict[o_lst[0]] == 'True' or 'true':
                    #         config.over_sampling = True
                    #     else:
                    #         logger.exception("Wrongly set over_sampling")
                    #         raise
                    if config.num_workers <= 0:
                        config.num_workers = available_cpu_cores - 2
                    if config.num_workers > available_cpu_cores:
                        logger.error(f"Only {available_cpu_cores} available to use.") 
                        config.num_workers = available_cpu_cores - 2
                    ##############################################################
                    DG = DataGneration(logger, config)
                    if config.high_ambiguous:
                        isolates_file, non_isolates_file, unique_seqs, read_max_len, read_min_len, genuine_df, negative_df, ambiguous_df, high_ambiguous_df = DG.data_files(edit_dis=1)
                    else:
                        isolates_file, non_isolates_file, unique_seqs, read_max_len, read_min_len, genuine_df, negative_df, ambiguous_df = DG.data_files(edit_dis=1)      
                    config.read_max_len = read_max_len
                    ###############################################################
                    EC = ErrorCorrection(logger, config)
                    ## one model to predict
                    if config.high_ambiguous:
                        # corrected_file, no_high_correct, new_negative_df = EC.all_in_one_correction(isolates_file, non_isolates_file, unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df)
                        corrected_file, new_negative_df = EC.all_in_one_correction(isolates_file, non_isolates_file, unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df)
                    else:
                        corrected_file, new_negative_df = EC.all_in_one_correction(isolates_file, non_isolates_file, unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df =None) 

                    if read_min_len > config.min_read_len:
                        genuine_df2, negative_df2, ambiguous_df2, unique_seqs2 = DG.extract_ed2_errors(corrected_file)
                        config.correct_data = EC.all_in_one_ed2_correction(corrected_file, unique_seqs2, genuine_df2, negative_df2, ambiguous_df2)
                    else:
                        config.correct_data = corrected_file

                    DataAnalysis(logger, config).evaluation()  
                    # delete bcool result
                    bcool_dir = os.path.join(config.result_dir, 'bcool/')
                    if os.path.exists(bcool_dir):
                        os.system("rm -rf %s" % bcool_dir)
                    ###############################################################
                    # if config.high_ambiguous:
                    #     if read_min_len > config.min_read_len:
                    #         genuine_df3, negative_df3, ambiguous_df3, unique_seqs3 = DG.extract_ed2_errors(no_high_correct)
                    #         config.correct_data = EC.all_in_one_ed2_correction(no_high_correct, unique_seqs3, genuine_df3, negative_df3, ambiguous_df3)
                    #     else:
                    #         config.correct_data = no_high_correct
                    # DataAnalysis(logger, config).evaluation()                    
############################################################################################################################
                elif module_arg == "amplicon_correction": 
                    if c_lst:
                        config = Config(opts_dict[c_lst[0]], logger) 
                        if i_lst:
                            config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Must set sequencing dataset in configuration.')
                            sys.exit()
                        if t_lst:
                            config.ground_truth_data = opts_dict[t_lst[0]] 
                            if not os.path.exists(config.ground_truth_data):
                                logger.exception('Ground truth data does not exsit.')
                                raise
                    elif i_lst:
                        config = Config(None, logger)  
                        config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Input file does not exsit.')
                            sys.exit()
                    else:
                        logger.error('Must input configuration file or sequencing dataset.')
                        sys.exit()

                    if t_lst:
                        config.ground_truth_data = opts_dict[t_lst[0]] 
                        if not os.path.exists(config.ground_truth_data):
                            logger.exception('Ground truth data does not exsit.')
                            raise
                    if d_lst:
                        config.result_dir = opts_dict[d_lst[0]] 
                    if p_lst:
                        config.num_workers = int(opts_dict[p_lst[0]])    
                    if a_lst:
                        config.high_ambiguous = opts_dict[a_lst[0]] 
                    if g_lst:
                        config.tree_method = opts_dict[g_lst[0]] 
                    # if o_lst:
                    #     if opts_dict[o_lst[0]] == 'False' or 'false':
                    #         config.over_sampling = False
                    #     elif opts_dict[o_lst[0]] == 'True' or 'true':
                    #         config.over_sampling = True
                    #     else:
                    #         logger.exception("Wrongly set over_sampling")
                    #         raise
                    if not os.path.exists(config.input_file):
                        logger.exception('Input file does not exsit.')
                        raise
                    if config.num_workers <= 0:
                        config.num_workers = available_cpu_cores - 2
                    if config.num_workers > available_cpu_cores:
                        logger.error(f"Only {available_cpu_cores} available to use.") 
                        config.num_workers = available_cpu_cores - 2
                    ##############################################################
                    DG = DataGneration(logger, config)
                    if config.high_ambiguous:
                        isolates_file, non_isolates_file, unique_seqs, read_max_len, read_min_len, genuine_df, negative_df, ambiguous_df, high_ambiguous_df = DG.data_files(edit_dis=1)
                    else:
                        isolates_file, non_isolates_file, unique_seqs, read_max_len, read_min_len, genuine_df, negative_df, ambiguous_df = DG.data_files(edit_dis=1)      
                    config.read_max_len = read_max_len
                    ###############################################################
                    EC = ErrorCorrection(logger, config)
                    ## one model to predict
                    if config.high_ambiguous:
                        corrected_file, new_negative_df = EC.all_in_one_correction(isolates_file, non_isolates_file, unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df)
                    else:
                        corrected_file, new_negative_df = EC.all_in_one_correction(isolates_file, non_isolates_file, unique_seqs, genuine_df, negative_df, ambiguous_df, high_ambiguous_df =None) 
                    if read_min_len > config.min_read_len:
                        genuine_df2, negative_df2, ambiguous_df2, unique_seqs2 = DG.extract_ed2_errors(corrected_file)
                        mid_result = EC.all_in_one_ed2_correction(corrected_file, unique_seqs2, genuine_df2, negative_df2, ambiguous_df2)
                    else:
                        mid_result = corrected_file

                    amplicon_df = DG.extract_amplicon_err_samples(mid_result, config.amplicon_low_freq, config.amplicon_high_freq)
                    config.correct_data = EC.correct_amplicon_err(mid_result, genuine_df, negative_df, new_negative_df, amplicon_df, config.amplicon_threshold_proba)
                
                    DataAnalysis(logger, config).evaluation()
                    # delete bcool result
                    bcool_dir = os.path.join(config.result_dir, 'bcool/')
                    if os.path.exists(bcool_dir):
                        os.system("rm -rf %s" % bcool_dir)
############################################################################################################################
                elif module_arg == "umi_correction": 
                    if c_lst:
                        config = Config(opts_dict[c_lst[0]], logger) 
                        if i_lst:
                            config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Must set sequencing dataset in configuration.')
                            sys.exit()
                        if t_lst:
                            config.ground_truth_data = opts_dict[t_lst[0]] 
                            if not os.path.exists(config.ground_truth_data):
                                logger.exception('Ground truth data does not exsit.')
                                raise
                    elif i_lst:
                        config = Config(None, logger)  
                        config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Input file does not exsit.')
                            sys.exit()
                    else:
                        logger.error('Must input configuration file or sequencing dataset.')
                        sys.exit()
                        
                    if t_lst:
                        config.ground_truth_data = opts_dict[t_lst[0]] 
                        if not os.path.exists(config.ground_truth_data):
                            logger.exception('Ground truth data does not exsit.')
                            raise
                    if d_lst:
                        config.result_dir = opts_dict[d_lst[0]] 
                    if p_lst:
                        config.num_workers = int(opts_dict[p_lst[0]])  
                    if a_lst:
                        config.high_ambiguous = opts_dict[a_lst[0]] 

                    if config.num_workers <= 0:
                        config.num_workers = available_cpu_cores - 2
                    if config.num_workers > available_cpu_cores:
                        logger.error(f"Only {available_cpu_cores} available to use.") 
                        config.num_workers = available_cpu_cores - 2
                    ##############################################################
                    config.high_ambiguous=False
                    DG = DataGneration(logger, config)
                    genuine_df = DG.extract_umi_genuine_errs()
                    # ##############################################################
                    EC = ErrorCorrection(logger, config)
                    config.correct_data = EC.umi_correction(config.input_file, genuine_df)
                    DataAnalysis(logger, config).evaluation()
############################################################################################################################
                elif module_arg == "mimic_umi":   
                    if i_lst and t_lst:
                        config = Config(None, logger)  
                        config.input_file = opts_dict[i_lst[0]]
                        config.ground_truth_data = opts_dict[t_lst[0]] 
                        if not os.path.exists(config.input_file):
                            logger.exception('Input file does not exsit.')
                            sys.exit()
                        if not os.path.exists(config.ground_truth_data):
                            logger.exception('Ground truth data does not exsit.')
                            sys.exit()
                    else:
                        logger.error("Must input raw and true datasets simultaneously or set them in the configuration file.")
                        sys.exit()

                    if d_lst:
                        config.result_dir = opts_dict[d_lst[0]] 

                    if config.num_workers <= 0:
                        config.num_workers = available_cpu_cores - 2
                    if config.num_workers > available_cpu_cores:
                        logger.error(f"Only {available_cpu_cores} available to use.") 
                        config.num_workers = available_cpu_cores - 2
                #############################################################
                    DP = DataProcessing( 
                        logger,
                        config.num_workers,
                        config.result_dir,
                        config.umi_start,
                        config.umi_end,
                        config.non_umi_start)
                    DP.write_mimic_umis(config.input_file, config.ground_truth_data)
############################################################################################################################
                elif module_arg == "real_umi":   
                    if i_lst:
                        config = Config(None, logger)  
                        config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Input file does not exsit.')
                            sys.exit()
                    else:
                        logger.error("Must input umi-based sequencing dataset or configuration file.")
                        sys.exit()

                    if u_lst:
                        config.umi_file = opts_dict[u_lst[0]] 
                    if d_lst:
                        config.result_dir = opts_dict[d_lst[0]] 
                    if config.num_workers <= 0:
                        config.num_workers = available_cpu_cores - 2
                    if config.num_workers > available_cpu_cores:
                        logger.error(f"Only {available_cpu_cores} available to use.") 
                        config.num_workers = available_cpu_cores - 2
                    DP = DataProcessing( 
                        logger,
                        config.num_workers,
                        config.result_dir,
                        config.umi_start,
                        config.umi_end,
                        config.non_umi_start)
                    DP.real_umi_data(config.input_file)                  
############################################################################################################################
                elif module_arg == "evaluation":   
                    if i_lst and t_lst and r_lst:
                        config = Config(None, logger)  
                        config.input_file = opts_dict[i_lst[0]]
                        config.ground_truth_data = opts_dict[t_lst[0]]
                        config.correct_data = opts_dict[r_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Input file does not exsit.')
                            sys.exit()
                        if not os.path.exists(config.ground_truth_data):
                            logger.exception('Ground truth data does not exsit.')
                            sys.exit()
                        if not os.path.exists(config.correct_data):
                            logger.exception('Corrected data does not exsit.')
                            sys.exit()      

                        if d_lst:
                            config.result_dir = opts_dict[d_lst[0]]
                    elif i_lst and r_lst:
                        config = Config(None, logger)  
                        config.input_file = opts_dict[i_lst[0]]
                        config.correct_data = opts_dict[r_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Input file does not exsit.')
                            sys.exit()
                        if not os.path.exists(config.correct_data):
                            logger.exception('Corrected data does not exsit.')
                            sys.exit()   
                        if d_lst:
                            config.result_dir = opts_dict[d_lst[0]]
                    else:
                        logger.error("Must input raw, true and rectified datasets or raw and rectified datasets simultaneously or set them in the configuration file.")
                        sys.exit()

                    if config.num_workers <= 0:
                        config.num_workers = available_cpu_cores - 2
                    if config.num_workers > available_cpu_cores:
                        logger.error(f"Only {available_cpu_cores} available to use.") 
                        config.num_workers = available_cpu_cores - 2
                #############################################################
                    DataAnalysis(logger, config).evaluation()
############################################################################################################################
                elif module_arg == "simulation":  
                    if c_lst:
                        config = Config(opts_dict[c_lst[0]], logger) 
                        if i_lst:
                            config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Must set sequencing dataset in configuration.')
                            sys.exit()
                    elif i_lst:
                        config = Config(None, logger)  
                        config.input_file = opts_dict[i_lst[0]]
                        if not os.path.exists(config.input_file):
                            logger.exception('Input file does not exsit.')
                            sys.exit()
                    else:
                        logger.error('Must input configuration file or sequencing dataset.')
                        sys.exit()

                    if d_lst:
                        config.result_dir = opts_dict[d_lst[0]] 
                    if a_lst:
                        config.high_ambiguous = opts_dict[a_lst[0]]
                    if g_lst:
                        config.tree_method = opts_dict[g_lst[0]]
                    # if o_lst:
                    #     if opts_dict[o_lst[0]] == 'False' or 'false':
                    #         config.over_sampling = False
                    #     elif opts_dict[o_lst[0]] == 'True' or 'true':
                    #         config.over_sampling = True
                    #     else:
                    #         logger.exception("Wrongly set over_sampling")
                    #         raise
                    if config.num_workers <= 0:
                        config.num_workers = available_cpu_cores - 2
                    if config.num_workers > available_cpu_cores:
                        logger.error(f"Only {available_cpu_cores} available to use.") 
                        config.num_workers = available_cpu_cores - 2
                    Simulation(logger, config).simulation()
############################################################################################################################
                # elif module_arg == "extract_isolates": 

                else:
                    # logger.error("Invalid module name, please check.")
                    logger.error("Input wrong module name, using command noise2read -h/--help for usage.")
                    sys.exit()
        
        else:
            logger.error("No valid arguments input, using command noise2read -h/--help for usage")
    except getopt.GetoptError as err:
        logger.error(err)
