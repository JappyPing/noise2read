# -*- coding: utf-8 -*-
# @Author: Pengyao Ping
# @Date:   2023-02-16 11:01:06
# @Last Modified by:   Pengyao Ping
# @Last Modified time: 2023-02-16 11:10:34
import os
from noise2read.utils import *

class SequenceCoverage():
    """
    A class to calculate sequence coverage using read datasets and genome.
    """
    def __init__(self, logger, config):
        """
        initialize the SequenceCoverage class

        Args:
            logger (class): customized logging
            config (class): parameters setting using configparser
        """
        self.logger = logger
        self.config = config
        if config.library_layout == "single" or "SR":
            self.logger.info(f'Read Set: {config.read}')
            bases = config.read.split('/')[-1]
            self.prefix = bases.split('.' + parse_file_type(config.read))[0]
        elif config.library_layout == "paired" or "PE":
            self.logger.info(f'Read Set 1: {config.read1}')
            self.logger.info(f'Read Set 2: {config.read2}')
            bases = config.read1.split('/')[-1]
            self.prefix = bases.split('.' + parse_file_type(config.read1))[0]
        if os.path.exists(self.config.result_dir):
            self.logger.info("Directory '% s' already exists, noise2read will use it" % self.config.result_dir)
        else:
            os.makedirs(self.config.result_dir)
            self.logger.info("Directory '% s' created" % self.config.result_dir)

    def build_index(self, reference):
        self.logger.info("Building genome index using bowtie2")
        ref_index_dir = os.path.join(self.config.result_dir, '/reference_index')
        if os.path.exists(ref_index_dir):
            self.logger.info("Directory '% s' already exists, noise2read will use it" % ref_index_dir)
        else:
            os.makedirs(ref_index_dir)
            self.logger.info("Directory '% s' created" % ref_index_dir)
        os.system("bowtie2-build --large-index %s %s" % (reference, ref_index_dir))
        return ref_index_dir

    def alignment(self, ref_index_dir):
        self.logger.info("Alignment using bowtie2")
        align_result = os.path.join(self.config.result_dir, self.prefix + '_alginment.sam')
        if self.config.library_layout == "single" or "SR":
            os.system("bowtie2 %s -p %s -x %s -0 %s > %s" % (self.config.Alignment, self.config.num_workers, ref_index_dir, self.config.read, align_result))
        elif self.config.library_layout == "paired" or "PE":
            os.system("bowtie2 %s -p %s -x %s -1 %s -2 %s > %s" % (self.config.Alignment, self.config.num_workers, ref_index_dir, self.config.read1, self.config.read2, align_result))
        return align_result

    