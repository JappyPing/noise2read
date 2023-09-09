# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-01-19 10:56:38
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-08-14 13:02:11

import configparser
import os

class Config(object):
    def __init__(self, config_file, logger):
        conf = configparser.ConfigParser()
        base_path = os.getcwd()
        if config_file:
            conf.read(config_file)
            # path
            if conf.has_option("Paths", "result_dir"):
                self.result_dir = conf.get("Paths", "result_dir")
            else:
                self.result_dir = os.path.join(base_path, 'result/')
                
            # input
            if conf.has_option("SourceInputData", "input_file"):
                self.input_file = conf.get("SourceInputData", "input_file")      
            if conf.has_option("SourceInputData", "umi_file"):
                self.umi_file = conf.get("SourceInputData", "umi_file")  

            if conf.has_option("SourceInputData", "ground_truth_data"):
                self.ground_truth_data = conf.get("SourceInputData", "ground_truth_data")
            else:
                self.ground_truth_data = None   
                
            if conf.has_option("SourceInputData", "correct_data"):
                self.correct_data = conf.get("SourceInputData", "correct_data")

            # setting model parameters
            if conf.has_option("General", "num_workers"):
                self.num_workers = conf.getint("General", "num_workers")
            else:
                self.num_workers = -1
            if conf.has_option("General", "chunks_num"):
                self.chunks_num = conf.getint("General", "chunks_num")
            else:   
                self.chunks_num = 100
            if conf.has_option("General", "verbose"):
                self.verbose = conf.getboolean("General", "verbose")
            else:
                self.verbose = False
            if conf.has_option("General", "min_iters"):
                self.min_iters = conf.getint("General", "min_iters")
            else:
                self.min_iters = 1000
            if conf.has_option("General", "iso_change_detail"):
                self.iso_change_detail = conf.getboolean("General", "iso_change_detail")
            else:
                self.iso_change_detail = False
            if conf.has_option("General", "top_n"):
                self.top_n = conf.getint("General", "top_n")
            else:
                self.top_n = 10 
            if conf.has_option("General", "min_read_len"):
                self.min_read_len = conf.getint("General", "min_read_len")
            else:
                self.min_read_len = 30 
            # when the number of negative samples larger than preseting threshold, noise2read will downsample negative samples for training
            if conf.has_option("General", "negative_sample_num"):
                self.negative_sample_num = conf.getint("General", "negative_sample_num")
            else:
                self.negative_sample_num = 300000 # when 
            # if conf.has_option("General", "over_sampling"):
            #     self.over_sampling = conf.getboolean("General", "over_sampling")
            # else:
            #     self.over_sampling = True
            # GraphSetup
            if conf.has_option("GraphSetup", "high_freq_thre"):
                self.high_freq_thre = conf.getint("GraphSetup", "high_freq_thre")
            else:
                self.high_freq_thre = 4

            if conf.has_option("GraphSetup", "max_error_freq"):
                self.max_error_freq = conf.getint("GraphSetup", "max_error_freq")
            else:
                self.max_error_freq = 4

            if conf.has_option("GraphSetup", "save_graph"):
                self.save_graph = conf.getboolean("GraphSetup", "save_graph")
            else:
                self.save_graph = False

            if conf.has_option("GraphSetup", "graph_visualization"):
                self.graph_visualization = conf.getboolean("GraphSetup", "graph_visualization")
            else:
                self.graph_visualization = False

            if conf.has_option("GraphSetup", "drawing_graph_num"):
                self.drawing_graph_num = conf.getint("GraphSetup", "drawing_graph_num")
            else:
                self.drawing_graph_num = 50   
            # EmbeddingSetup
            if conf.has_option("EmbeddingSetup", "entropy_kmer"):
                self.entropy_kmer = conf.getint("EmbeddingSetup", "entropy_kmer")
            else:
                self.entropy_kmer = 3  
            if conf.has_option("EmbeddingSetup", "entropy_q"):
                self.entropy_q = conf.getint("EmbeddingSetup", "entropy_q")
            else:
                self.entropy_q = 2  
            if conf.has_option("EmbeddingSetup", "kmer_freq"):
                self.kmer_freq = conf.getint("EmbeddingSetup", "kmer_freq")
            else:
                self.kmer_freq = 3  
            if conf.has_option("EmbeddingSetup", "read_type"):
                self.read_type = conf.get("EmbeddingSetup", "read_type")
            else:
                self.read_type = 'DNA'

            # AmbiguousSetup
            # if conf.has_option("AmbiguousSetup", "ambiguous_error_node_degree"):
            #     self.ambiguous_error_node_degree = conf.getint("AmbiguousSetup", "ambiguous_error_node_degree")
            # else:
            #     self.ambiguous_error_node_degree = 4 # default 
            if conf.has_option("AmbiguousSetup", "high_ambiguous"):
                self.high_ambiguous = conf.getboolean("AmbiguousSetup", "high_ambiguous")
            else:
                self.high_ambiguous = True        
            if conf.has_option("AmbiguousSetup", "proba_deviation"):
                self.proba_deviation = conf.getfloat("AmbiguousSetup", "proba_deviation")
            else:
                self.proba_deviation = 0.95 # default 

            if conf.has_option("AmbiguousSetup", "iso_neg_high"):
                self.iso_neg_high = conf.getboolean("AmbiguousSetup", "iso_neg_high")
            else:
                self.iso_neg_high = True    

            # ModelTuningSetup
            if conf.has_option("ModelTuningSetup", "n_trials"):
                self.n_trials = conf.getint("ModelTuningSetup", "n_trials")
            else:
                self.n_trials = 20
                
            if conf.has_option("ModelTuningSetup", "n_estimators"):
                self.n_estimators = conf.getint("ModelTuningSetup", "n_estimators")
            else:
                self.n_estimators = 400

            if conf.has_option("ModelTuningSetup", "test_size"):
                self.test_size = conf.getfloat("ModelTuningSetup", "test_size")
            else:
                self.test_size = 0.1 # default        
            if conf.has_option("ModelTuningSetup", "random_state"):
                self.random_state = conf.getint("ModelTuningSetup", "random_state")
            else:
                self.random_state = 42 # default  

            if conf.has_option("ModelTuningSetup", "tree_method"):
                self.tree_method = conf.get("ModelTuningSetup", "tree_method")
            else:
                self.tree_method = 'auto'            

            if conf.has_option("ModelTuningSetup", "learning_rate_min"):
                self.learning_rate_min = conf.getfloat("ModelTuningSetup", "learning_rate_min")
            else:
                self.learning_rate_min = 1e-3 # default     
            if conf.has_option("ModelTuningSetup", "learning_rate_max"):
                self.learning_rate_max = conf.getfloat("ModelTuningSetup", "learning_rate_max")
            else:
                self.learning_rate_max = 1e-1 # default 

            if conf.has_option("ModelTuningSetup", "max_depth_min"):
                self.max_depth_min = conf.getint("ModelTuningSetup", "max_depth_min")
            else:
                self.max_depth_min = 3 # default     
            if conf.has_option("ModelTuningSetup", "max_depth_max"):
                self.max_depth_max = conf.getint("ModelTuningSetup", "max_depth_max")
            else:
                self.max_depth_max = 15 # default     
            if conf.has_option("ModelTuningSetup", "max_depth_step"):
                self.max_depth_step = conf.getint("ModelTuningSetup", "max_depth_step")
            else:
                self.max_depth_step = 1 # default 

            if conf.has_option("ModelTuningSetup", "num_boost_round_min"):
                self.num_boost_round_min = conf.getint("ModelTuningSetup", "num_boost_round_min")
            else:
                self.num_boost_round_min = 200 # default     
            if conf.has_option("ModelTuningSetup", "num_boost_round_max"):
                self.num_boost_round_max = conf.getint("ModelTuningSetup", "num_boost_round_max")
            else:
                self.num_boost_round_max = 300 # default     
            if conf.has_option("ModelTuningSetup", "num_boost_round_step"):
                self.num_boost_round_step = conf.getint("ModelTuningSetup", "num_boost_round_step")
            else:
                self.num_boost_round_step = 10 # default 

            if conf.has_option("ModelTuningSetup", "subsample_min"):
                self.subsample_min = conf.getfloat("ModelTuningSetup", "subsample_min")
            else:
                self.subsample_min = 0.8 # default     
            if conf.has_option("ModelTuningSetup", "subsample_max"):
                self.subsample_max = conf.getfloat("ModelTuningSetup", "subsample_max")
            else:
                self.subsample_max = 1.0 # default     

            if conf.has_option("ModelTuningSetup", "colsample_bytree_min"):
                self.colsample_bytree_min = conf.getfloat("ModelTuningSetup", "colsample_bytree_min")
            else:
                self.colsample_bytree_min = 0.8 # default     
            if conf.has_option("ModelTuningSetup", "colsample_bytree_max"):
                self.colsample_bytree_max = conf.getfloat("ModelTuningSetup", "colsample_bytree_max")
            else:
                self.colsample_bytree_max = 1.0 # default     
            
            if conf.has_option("ModelTuningSetup", "verbose_eval"):
                self.verbose_eval = conf.getboolean("ModelTuningSetup", "verbose_eval")
            else:
                self.verbose_eval = False
                
            # xgboostclassifier seed
            if conf.has_option("ModelTuningSetup", "xgboost_seed"):
                self.xgboost_seed = conf.getint("ModelTuningSetup", "xgboost_seed")
            else:
                self.xgboost_seed = 42 # default 

            if conf.has_option("ModelTuningSetup", "optuna_seed"):
                self.optuna_seed = conf.getint("ModelTuningSetup", "optuna_seed")
            else:
                self.optuna_seed = 42 # default 

            if conf.has_option("ModelTuningSetup", "best_accuracy"):
                self.best_accuracy = conf.getfloat("ModelTuningSetup", "best_accuracy")
            else:
                self.best_accuracy = 0.85 # default     

            # real umi
            if conf.has_option("RealUMI", "umi_in_read"):
                self.umi_in_read = conf.getboolean("RealUMI", "umi_in_read")
            else:
                self.umi_in_read = False
            if conf.has_option("RealUMI", "umi_start"):
                self.umi_start = conf.getint("RealUMI", "umi_start")
            else:
                self.umi_start = 0 # default 
            if conf.has_option("RealUMI", "umi_end"):
                self.umi_end = conf.getint("RealUMI", "umi_end")
            else:
                self.umi_end = 12 # default 
            if conf.has_option("RealUMI", "non_umi_start"):
                self.non_umi_start = conf.getint("RealUMI", "non_umi_start")
            else:
                self.non_umi_start = 24 # default

            if conf.has_option("RealUMI", "group_read_number"):
                self.group_read_number = conf.getint("RealUMI", "group_read_number")
            else:
                self.group_read_number = 10 # default
            if conf.has_option("RealUMI", "read_edit_dif"):
                self.read_edit_dif = conf.getint("RealUMI", "read_edit_dif")
            else:
                self.read_edit_dif = 2 # default

            ###
            if conf.has_option("RealUMI", "separator1"):
                self.separator1 = conf.get("RealUMI", "separator1")
            else:
                self.separator1 = '_'  
            if conf.has_option("RealUMI", "separator1_idx"):
                self.separator1_idx = conf.getint("RealUMI", "separator1_idx")
            else:
                self.separator1_idx = 2 # default

            if conf.has_option("RealUMI", "separator2"):
                self.separator1 = conf.get("RealUMI", "separator2")
            else:
                self.separator1 = ' '  
            if conf.has_option("RealUMI", "separator2_idx"):
                self.separator2_idx = conf.getint("RealUMI", "separator2_idx")
            else:
                self.separator2_idx = 0 # default

            # amplicon
            if conf.has_option("Amplicon", "amplicon_low_freq"):
                self.amplicon_low_freq = conf.getint("Amplicon", "amplicon_low_freq")
            else:
                self.amplicon_low_freq = 50 # default             
            if conf.has_option("Amplicon", "amplicon_high_freq"):
                self.amplicon_high_freq = conf.getint("Amplicon", "amplicon_high_freq")
            else:
                self.amplicon_high_freq = 1500 # default   
            if conf.has_option("Amplicon", "amplicon_threshold_proba"):
                self.amplicon_threshold_proba = conf.getfloat("Amplicon", "amplicon_threshold_proba")
            else:
                self.amplicon_threshold_proba = 0.85 # default   
            # if conf.has_option("Amplicon", "amplicon_error_node_degree"):
            #     self.amplicon_error_node_degree = conf.getint("Amplicon", "amplicon_error_node_degree")
            # else:
            #     self.amplicon_error_node_degree = 4 # default   

            # Simulation
            if conf.has_option("Simulation", "min_freq"):
                self.min_freq = conf.getint("Simulation", "min_freq")
            else:
                self.min_freq = 4 # default
            if conf.has_option("Simulation", "min_read_count"):
                self.min_read_count = conf.getint("Simulation", "min_read_count")
            else:
                self.min_read_count = 30 # default
                
            # if conf.has_option("Simulation", "substations"):
            #     self.substations = conf.getboolean("Simulation", "substations")
            # else:
            #     self.substations = True
            # if conf.has_option("Simulation", "indels"):
            #     self.indels = conf.getboolean("Simulation", "indels")
            # else:
            #     self.indels = False

            if conf.has_option("Simulation", "error_rate1"):
                self.error_rate1 = conf.getfloat("Simulation", "error_rate1")
            else:
                self.error_rate1 = 0.03 # default

            if conf.has_option("Simulation", "error_rate2"):
                self.error_rate2 = conf.getfloat("Simulation", "error_rate2")
            else:
                self.error_rate2 = 0.005 # default


            # # Evaluation
            # if conf.has_option("Evaluation", "delta"):
            #     self.delta = conf.getint("Evaluation", "delta")
            # else:
            #     self.delta = 10 # default

        if config_file == None:
            # path
            self.result_dir = os.path.join(base_path, 'result/')
            # input
            self.ground_truth_data = None  

            # general
            self.num_workers = -1 
            self.chunks_num = 100
            self.min_iters = 1000
            self.verbose = True 
            self.iso_change_detail = False     
            self.top_n = 10      
            # self.over_sampling = True 
            self.negative_sample_num = 300000
            self.min_read_len = 30

            # GraphSetup
            self.high_freq_thre = 4
            self.max_error_freq = 4
            self.save_graph = False
            self.graph_visualization = False
            self.drawing_graph_num = 50   

            # EmbeddingSetup
            self.entropy_kmer = 3
            self.entropy_q = 2
            self.kmer_freq = 3
            self.read_type = 'DNA'

            # AmbiguousSetup
            self.high_ambiguous = True 
            # high ambiguous predict probability difference
            self.proba_deviation = 0.95      # for base editing and lncRNA quant datasets 0.75 others 0.6 
            # self.ambiguous_error_node_degree = 4   # for base editing and lncRNA quant datasets 3 others 4 
            self.iso_neg_high = True # if True, the high frequency isolated nodes aslso included as negative samples for high ambiguous prediction. This will use quite a lot computational resources (memory) for embeeding and model training

            # ModelTuningSetup
            self.n_trials = 20
            self.n_estimators = 400
            self.test_size = 0.1 # default 
            # random state for SMOTE and train_test_split       
            self.random_state = 42 # default  
            self.tree_method = 'auto'
            self.learning_rate_min = 1e-3 # default     
            self.learning_rate_max = 1e-1 # default 
            self.max_depth_min = 3 # default     
            self.max_depth_max = 15 # default     
            self.max_depth_step = 1 # default 
            self.num_boost_round_min = 200 # default     
            self.num_boost_round_max = 300 # default     
            self.num_boost_round_step = 10 # default 
            self.subsample_min = 0.8 # default     
            self.subsample_max = 1 # default     
            self.colsample_bytree_min = 0.8 # default     
            self.colsample_bytree_max = 1 # default     
            self.verbose_eval = False
            # xgboostclassifier seed
            self.xgboost_seed = 42 # default 
            self.optuna_seed = 42
            # optuna best trial accuracy
            self.best_accuracy = 0.85

            # real umi
            self.umi_in_read = False
            # umi in read
            self.umi_start = 0
            self.umi_end = 12
            self.non_umi_start = 24
            self.group_read_number = 10
            self.read_edit_dif = 2
            # umi in name
            self.separator1 = '_'
            self.separator1_idx = 2
            self.separator2 = ' '
            self.separator2_idx = 0

            #amplicon
            self.amplicon_low_freq = 50
            self.amplicon_high_freq = 1500
            self.amplicon_threshold_proba = 0.85
            # self.amplicon_error_node_degree = 4

            # simulation
            self.min_freq = 4
            self.min_read_count = 30
            # self.substations = True
            # self.indels = False
            self.error_rate1 = 0.03
            self.error_rate2 = 0.005
  
            # # Evaluation
            # self.delta = 1

            # # coverage
            # self.library_layout = 'PE'
            # self.Alignment = '--local' # '--end-to-end'
            
            
