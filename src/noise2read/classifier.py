# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-05-16 18:27:41

import optuna
from optuna.samplers import TPESampler
from sklearn.model_selection import train_test_split
import os
import xgboost as xgb
import numpy as np
# from imblearn.combine import SMOTEENN
from imblearn.over_sampling import SMOTE 
import collections
# from imblearn.under_sampling import RandomUnderSampler
from matplotlib import pyplot
from sklearn.metrics import f1_score
import optuna
import numpy as np
# from matplotlib import pyplot as plt
# from sklearn.metrics import classification_report, confusion_matrix
from noise2read.utils import MemoryMonitor

class MLClassifier:
    """
    A classifier based on xgboost and optuna to predict whether a high-/low-frequency read is an error read mutated from a read when sequencing.
    """
    def __init__(self, logger, config, study_name, data, label, ambi_data):
        """
        initialize the MLClassifier class

        Args:
            logger (class): customized logging
            config (class): parameters setting using configparser
            study_name (str): study name of optuna
            data (arrays): numpy arrays with same length / shape[0]
            label (arrays): numpy arrays, sample labels (n_samples, n_classes(0/1))
            ambi_data (arrays): numpy arrays with same length / shape[0]
        """
        self.logger = logger
        self.config = config
        self.study_name = study_name
        self.data = data
        self.label = label
        self.ambi_data = ambi_data
        self.MM = MemoryMonitor(self.logger)

        # Shuffle the data in a deterministic way
        np.random.seed(self.config.random_state)
        indices = np.random.permutation(len(self.data))
        shuffled_data = self.data[indices]
        shuffled_labels = self.label[indices]
        self.x_train, self.x_test, self.y_train, self.y_test = train_test_split(shuffled_data, shuffled_labels, test_size=self.config.test_size, shuffle=False, random_state=self.config.random_state)

        # self.x_train, self.x_test, self.y_train, self.y_test = train_test_split(self.data, self.label, test_size=self.config.test_size, shuffle=True, random_state=self.config.random_state)

        # self.x_test, self.x_val, self.y_test, self.y_val, self.x_test_weight, self.x_val_weight = train_test_split(test_val_x, test_val_y, test_val_weight, test_size=0.50, random_state=self.config.random_state, shuffle=True)
        # self.x_test_weight, self.x_val_weight = train_test_split(test_val_weight, test_size=0.50, random_state=self.config.random_state, shuffle=False)
        # self.logger.debug(type(self.x_val_weight))
        # self.logger.debug(f'{self.x_val}, {self.y_val.shape}, {len(self.x_val_weight)}')
        self.logger.info(f'The number of negative and positive samples: {sorted(collections.Counter(self.y_train).items())}')

        sm = SMOTE(random_state=self.config.random_state)
        self.X_resampled, self.y_resampled = sm.fit_resample(self.x_train, self.y_train)
        self.logger.info(f'After over-sampling: {sorted(collections.Counter(self.y_resampled).items())}')
        del sm
        self.MM.measure()
        # rus = RandomUnderSampler(random_state=42)
        # self.X_resampled, self.y_resampled = rus.fit_resample(self.x_train, self.y_train)
        # self.logger.info(f'After under-sampling {sorted(collections.Counter(self.y_resampled).items())}')

        # smote_enn = SMOTEENN(random_state=0)
        # self.X_resampled, self.y_resampled = smote_enn.fit_resample(self.x_train, self.y_train)
        # self.logger.info(f'Combination of over- and under-sampling {sorted(collections.Counter(self.y_resampled).items())}')
 
    def objective(self, trial):
        """
        Define an objective function to be minimized in the framewokr of Optuna

        Args:
            trial (class): The trial module of Optuna contains Trial related classes and functions. A Trial instance represents a process of evaluating an objective function. More details please see Optuna API reference: https://optuna.readthedocs.io/en/stable/reference/trial.html
            
        Raises:
            optuna.TrialPruned: Exception for pruned trials of Optuna, more details see https://optuna.readthedocs.io/en/stable/reference/generated/optuna.TrialPruned.html#optuna.TrialPruned

        Returns:
            float: test data accuracy
        """
        # X, y = shuffle(self.data, self.label)
        param = {
            "verbosity": 0,
            # 'gpu_id': 0,
            # "max_delta_step": 1,
            "tree_method": self.config.tree_method,
            "objective": "binary:logistic",
            'n_estimators': self.config.n_estimators,
            'early_stopping_rounds': 100,
            # "eval_metric": ['auc', 'error'],
            # "eval_metric": ['logloss'],
            "eval_metric": ['auc', 'logloss'],
            # "booster": trial.suggest_categorical("booster", ["gbtree", "dart"]),
            "booster": trial.suggest_categorical("booster", ["gbtree", "gblinear", "dart"]),
            "lambda": trial.suggest_float("lambda", 1e-8, 1.0, log=False),
            "alpha": trial.suggest_float("alpha", 1e-8, 1.0, log=False),
            "learning_rate": trial.suggest_float("learning_rate", self.config.learning_rate_min, self.config.learning_rate_max, log=False),
            "max_depth": trial.suggest_int("max_depth", self.config.max_depth_min, self.config.max_depth_max, step = self.config.max_depth_step, log = False),
            "subsample": trial.suggest_float('subsample', self.config.subsample_min, self.config.subsample_max),
            "colsample_bytree": trial.suggest_float('colsample_bytree', self.config.colsample_bytree_min, self.config.colsample_bytree_max),
            # the following setting only for gbtree and dart
            # "eta": trial.suggest_float("eta", 1e-8, 1.0, log=False),
            # "gamma": trial.suggest_float("gamma", 1e-8, 1.0, log=False),
            # "grow_policy": trial.suggest_categorical("grow_policy", ["depthwise", "lossguide"]),
            # "num_boost_round":  trial.suggest_int("num_boost_round", self.config.num_boost_round_min, self.config.num_boost_round_max, step = self.config.num_boost_round_step, log=False)
            # "num_boost_round": 50
        }

        if param["booster"] == "gbtree" or param["booster"] == "dart":
            param["max_depth"] = trial.suggest_int("max_depth", self.config.max_depth_min, self.config.max_depth_max, step = self.config.max_depth_step, log = False)
            param["eta"] = trial.suggest_float("eta", 1e-8, 1.0, log=False)
            param["gamma"] = trial.suggest_float("gamma", 1e-8, 1.0, log=False)
            param["grow_policy"] = trial.suggest_categorical("grow_policy", ["depthwise", "lossguide"])
        if param["booster"] == "dart":
            param["sample_type"] = trial.suggest_categorical("sample_type", ["uniform", "weighted"])
            param["normalize_type"] = trial.suggest_categorical("normalize_type", ["tree", "forest"])
            param["rate_drop"] = trial.suggest_float("rate_drop", 1e-8, 1.0, log=False)
            param["skip_drop"] = trial.suggest_float("skip_drop", 1e-8, 1.0, log=False)
        
        pruning_callback = optuna.integration.XGBoostPruningCallback(trial, "validation_0-auc")
        # pruning_callback = optuna.integration.XGBoostPruningCallback(trial, "validation_1-logloss")
        # pruning_callback = optuna.integration.XGBoostPruningCallback(trial, "validation_1-aucpr")

        # y_count = sorted(collections.Counter(self.y_train).items())
        # _, neg_num = y_count[0]
        # _, pos_num = y_count[1]
        # xgbc = xgb.XGBClassifier(**param, pruning_callback=[pruning_callback], probability=True, scale_pos_weight= neg_num/pos_num)
        # xgbc.fit(self.x_train, self.y_train, eval_set=[(self.x_train, self.y_train),(self.x_test, self.y_test)],  verbose=True)
        # train_accuracy = xgbc.score(self.x_train, self.y_train)

        xgbc = xgb.XGBClassifier(**param, pruning_callback=[pruning_callback], probability=True, seed=self.config.xgboost_seed)
        xgbc.fit(self.X_resampled, self.y_resampled, eval_set=[(self.X_resampled, self.y_resampled), (self.x_test, self.y_test)],  verbose=self.config.verbose_eval)
        
        # xgbc.fit(self.x_train, self.y_train, eval_set=[(self.x_train, self.y_train),(self.x_test, self.y_test)],  verbose=True)
        # xgbc.fit(self.x_train, self.y_train, eval_set=[(self.x_val, self.y_val)], verbose=False)
        # xgbc.fit(self.x_train, self.y_train, sample_weight = self.x_train_weight, eval_set=[(self.x_val, self.y_val)], sample_weight_eval_set = [self.x_val_weight],  verbose=False)

        # accuracy = xgbc.score(self.x_test, self.y_test, sample_weight=self.x_test_weight)
        # train_accuracy = xgbc.score(self.X_resampled, self.y_resampled)
        train_accuracy = xgbc.score(self.x_train, self.y_train)
        test_accuracy = xgbc.score(self.x_test, self.y_test)
        # preds = xgbc.predict(self.x_test)
        # pred_labels = np.rint(preds)
        # accuracy = accuracy_score(self.y_test, pred_labels)

        # train_f1 = f1_score(self.y_resampled, xgbc.predict(self.X_resampled))
        train_f1 = f1_score(self.y_train, xgbc.predict(self.x_train))
        test_f1 = f1_score(self.y_test, xgbc.predict(self.x_test))
        
        trial.report(test_accuracy, trial.number)
        if trial.should_prune():
            raise optuna.TrialPruned()
        
        self.logger.info( " Trial " + str(trial.number))
        self.logger.info("      Train Accuracy: {}, Test Accuracy: {}".format(train_accuracy, test_accuracy))
        self.logger.info("      Train F1: {}, Test F1: {}".format(train_f1, test_f1))
        trial.set_user_attr(key="best_model", value=xgbc) # save model

        # results = xgbc.evals_result()
        # test_logloss = results['validation_1']['logloss']
        # test_auc = results['validation_1']['auc']
        # return test_logloss
        return test_accuracy
        # return test_f1
        
    def tunning(self, n_trials = 150):
        """
        function to create a Optuna study for model training and prediction.

        Args:
            n_trials (int, optional): Number of Optuan trials. Defaults to 150.

        Returns:
            array: numpy ndarray of shape (n_samples, n_classes)
            Estimated probabilities.
        """
        # study = optuna.create_study(study_name = self.study_name,
        #     pruner=optuna.pruners.MedianPruner(n_warmup_steps=5), direction="maximize" 
        # ) #, interval_steps=10 n_startup_trials=5, 
        self.logger.info("-------------------------------------------------------------")
        sampler = TPESampler(seed=self.config.optuna_seed)  # Make the sampler behave in a deterministic way.
        study = optuna.create_study(study_name = self.study_name, direction="maximize", sampler=sampler)
        study.optimize(self.objective, n_trials, show_progress_bar=False, gc_after_trial=True)
        # print(study.best_trial)
        
        self.logger.info(f'Study Name: {self.study_name}')
        self.logger.info('Number of finished trials: {}'.format(len(study.trials)))
        self.logger.info('Best trial:')
        best_trial = study.best_trial
        self.logger.info('  Test Accuracy: {}'.format(best_trial.value))
        self.logger.info('  Params: ')
        for key, value in best_trial.params.items():
            self.logger.info('    {}: {}'.format(key, value))
        self.MM.measure()
        # best_trial_copy = copy.deepcopy(best_trial)

        # re-evaluate
        # self.objective(best_trial)
        # # # the user attribute is overwritten by re-evaluation
        # assert best_trial.user_attrs != best_trial_copy.user_attrs
        # if best_trial.value < self.config.best_accuracy:
        #     self.logger.warning(f"The stduy {self.study_name} output lower best trial accuracy of {best_trial.value} than the pre-set self.config.best_accuracy, which may result in bad error correction performance.")
        # fig1 = optuna.visualization.plot_optimization_history(study)
        # fig1.write_image(os.path.join(self.config.result_dir, "optimization_history.png"))

        # fig2 = optuna.visualization.plot_intermediate_values(study)
        # fig2.write_image(os.path.join(self.config.result_dir, "intermediate_values.png"))

        # fig3 = optuna.visualization.plot_parallel_coordinate(study)
        # fig3.write_image(os.path.join(self.config.result_dir, "parallel_coordinate.png"))

        # fig4 = optuna.visualization.plot_contour(study)
        # fig4.write_image(os.path.join(self.config.result_dir, "contour.png"))
        best_model = best_trial.user_attrs["best_model"]
        #prediction
        results = best_model.evals_result()
        
        # plot learning curves
        fig, ax = pyplot.subplots()
        
        ax.set_xlabel('estimators')
        ax.set_ylabel('logloss')
        
        ax.plot(results['validation_0']['logloss'], label='train')
        ax.plot(results['validation_1']['logloss'], label='test')
        # show the legend
        ax.legend()
        # show the plot
        fig.savefig(os.path.join(self.config.result_dir, self.study_name + '_train-test-logloss.png'))
        predictions = best_model.predict_proba(self.ambi_data)[:, 1]
        pyplot.clsoe()
        self.logger.info("-------------------------------------------------------------")
        del study, sampler, results, best_model
        self.MM.measure()
        self.MM.stop()
        return predictions