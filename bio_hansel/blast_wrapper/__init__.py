# -*- coding: utf-8 -*-

import attr
from datetime import datetime
import logging
import shutil
import os
import pandas as pd
import numpy as np
from pandas.io.common import EmptyDataError
import re
from .const import BLAST_TABLE_COLS
from ..utils import exc_exists, run_command


@attr.s
class BlastRunner:
    fasta_path = attr.ib()
    tmp_work_dir = attr.ib(default='/tmp', validator=attr.validators.instance_of(str))
    blast_db_created = attr.ib(default=False, validator=attr.validators.instance_of(bool))
    makeblastdb_exc = attr.ib(default='makeblastdb', validator=attr.validators.instance_of(str))
    blastn_exc = attr.ib(default='blastn', validator=attr.validators.instance_of(str))

    @makeblastdb_exc.validator
    def _check_makeblastdb_exists(self, attribute, value):
        if not exc_exists(value):
            raise OSError('makeblastdb executable "{}" does not exist in the user $PATH'.format(value))

    @blastn_exc.validator
    def _check_blastn_exists(self, attribute, value):
        if not exc_exists(value):
            raise OSError('blast executable "{}" does not exist in the user $PATH'.format(value))

    @fasta_path.validator
    def _fasta_path_exists(self, attribute, value):
        if not os.path.exists(value):
            raise OSError('FASTA file does not exist at {}'.format(value))

    def _create_tmp_folder(self):
        count = 1
        tmp_dir = self.tmp_work_dir
        while True:
            try:
                logging.info('Trying to create analysis directory at: %s', tmp_dir)
                os.makedirs(tmp_dir)
                break
            except OSError as e:
                logging.warning('Error on creation of tmp analysis directory "{}"! {}'.format(
                    tmp_dir,
                    e
                ))
                tmp_dir = '{}_{}'.format(self.tmp_work_dir, count)
                count += 1
        self.tmp_work_dir = tmp_dir
        return self.tmp_work_dir

    def _copy_fasta_to_work_dir(self):
        filename = os.path.basename(self.fasta_path)
        filename, ext = os.path.splitext(filename)
        filename_no_spaces = re.sub(r'\W', '_', filename)
        dest_path = os.path.join(self.tmp_work_dir, filename_no_spaces + ext)
        if self.fasta_path == dest_path:
            self.tmp_fasta_path = dest_path
            return dest_path
        shutil.copyfile(self.fasta_path, dest_path)
        self.tmp_fasta_path = dest_path
        return dest_path

    def _run_makeblastdb(self):
        work_dir = os.path.dirname(self.tmp_fasta_path)
        filename = os.path.basename(self.tmp_fasta_path)
        nin_filepath = os.path.join(work_dir, filename + '.nin')
        if os.path.exists(nin_filepath):
            self.blast_db_created = True
            return self.tmp_fasta_path
        cmdlist = [self.makeblastdb_exc,
                   '-in', '{}'.format(self.tmp_fasta_path),
                   '-dbtype', 'nucl']
        exit_code, stdout, stderr = run_command(cmdlist)
        if exit_code != 0:
            raise Exception('Error {}: makeblastdb could not create a BLAST DB for {}. stderr: {}'.format(exit_code,
                                                                                                          self.tmp_fasta_path,
                                                                                                          stderr))
        if stdout is not None and stdout != '':
            logging.debug('makeblastdb on {0} STDOUT: {1}'.format(self.tmp_fasta_path, stdout))
        if stderr is not None and stderr != '':
            logging.debug('makeblastdb on {0} STDERR: {1}'.format(self.tmp_fasta_path, stderr))
        if os.path.exists(nin_filepath):
            self.blast_db_created = True
            return self.tmp_fasta_path
        else:
            ex_msg = 'makeblastdb was not able to create a BLAST DB for {0}. STDERR: {1}'.format(filename, stderr)
            logging.error(ex_msg)
            raise Exception(ex_msg)

    def blast_against_query(self, query_fasta_path, blast_task='megablast', evalue=1e-4, min_pid=85, word_size=22):
        if not self.blast_db_created:
            self.prep_blast()

        query_filename = os.path.basename(query_fasta_path)
        db_filename = os.path.basename(self.tmp_fasta_path)
        timestamp = '{:%Y%b%d_%H_%M_%S}'.format(datetime.now())
        outfile = os.path.join(self.tmp_work_dir, '{}-{}-{}.blast'.format(query_filename,
                                                                          db_filename,
                                                                          timestamp))
        cmd_list = [self.blastn_exc,
                    '-task', blast_task,
                    '-query', query_fasta_path,
                    '-db', '{}'.format(self.tmp_fasta_path),
                    '-word_size', '{}'.format(word_size),
                    '-evalue', '{}'.format(evalue),
                    '-dust', 'no',
                    '-perc_identity', '{}'.format(min_pid),
                    '-out', outfile,
                    '-outfmt', '6 {}'.format(' '.join(BLAST_TABLE_COLS))]
        logging.info('Running commandline "{}"'.format(' '.join(cmd_list)))
        exit_code, stdout, stderr = run_command(cmd_list)
        if os.path.exists(outfile):
            return outfile
        else:
            err_msg_fmt = 'error code {}: blastn on db {} and query {} did not produce expected output file at {}. stderr: {}'
            ex_msg = err_msg_fmt.format(exit_code, db_filename, query_filename, outfile, stderr)
            raise Exception(ex_msg)

    def cleanup(self):
        self.blast_db_created = False
        shutil.rmtree(self.tmp_work_dir)

    def prep_blast(self):
        self._create_tmp_folder()
        self._copy_fasta_to_work_dir()
        self._run_makeblastdb()

    def __enter__(self):
        if self.blast_db_created:
            return self
        self.prep_blast()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()


@attr.s
class BlastReader:
    blast_outfile = attr.ib(validator=attr.validators.instance_of(str))
    df = attr.ib(default=None, validator=attr.validators.optional(attr.validators.instance_of(pd.DataFrame)))

    def parse(self):
        """Parse tabular blastn output file into a pandas DataFrame

        Sort the DataFrame by BLAST bitscore, compute query coverage and if result is truncated by end of subject
        sequence/contig.

        Returns:
            pandas.DataFrame: dataframe of tabular BLAST results
            None: if no results could be parsed from BLAST output file

        Exceptions:
            EmptyDataError: No data could be parsed from the `blastn` output file
        """
        try:
            self.df = pd.read_table(self.blast_outfile, header=None)
            self.df.columns = BLAST_TABLE_COLS
            # calculate the coverage for when results need to be validated
            self.df.loc[:, 'coverage'] = self.df.length / self.df.qlen
            self.df.sort_values(by='bitscore', ascending=False, inplace=True)
            self.df.loc[:, 'is_trunc'] = BlastReader.trunc(qstart=self.df.qstart,
                                                           qend=self.df.qend,
                                                           qlen=self.df.qlen,
                                                           sstart=self.df.sstart,
                                                           send=self.df.send,
                                                           slen=self.df.slen)

            return self.df

        except EmptyDataError as exc:
            logging.warning('No BLASTN results to parse from file %s', self.blast_outfile)
            return None

    def to_dict(self):
        if self.df is not None:
            return self.df.to_dict()
        else:
            return None

    @staticmethod
    def trunc(qstart, qend, sstart, send, qlen, slen):
        """Check if a query sequence is truncated by the end of a subject sequence

        Args:
            qstart (int pandas.Series): Query sequence start index
            qend (int pandas.Series): Query sequence end index
            sstart (int pandas.Series): Subject sequence start index
            send (int pandas.Series): Subject sequence end index
            qlen (int pandas.Series): Query sequence length
            slen (int pandas.Series): Subject sequence length

        Returns:
            Boolean pandas.Series: Result truncated by subject sequence end?
        """
        ssum2 = (send + sstart) / 2.0
        sabs2 = np.abs(send - sstart) / 2.0
        smax = ssum2 + sabs2
        smin = ssum2 - sabs2
        q_match_len = np.abs(qstart - qend) + 1
        return (q_match_len < qlen) & ((smax >= slen) | (smin <= 1))

    def perfect_matches(self):
        """
        Return pandas DataFrame with perfect BLAST matches (100% identity and coverage)

        Returns:
            pandas.DataFrame or None: DataFrame of perfect BLAST matches or None if no perfect matches exist
        """
        if self.df is None:
            return None

        df_perfect_matches = self.df[(self.df['coverage'] == 1.0) & (self.df['pident'] == 100.0)]
        if df_perfect_matches.shape[0] == 0:
            return None
        return df_perfect_matches

    def __enter__(self):
        if self.blast_outfile is None:
            return self
        self.parse()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.df = None
        self.blast_outfile = None
