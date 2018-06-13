#!/usr/bin/env python

import os, sys
import click
import random
import string

from pyflux import FluxWorkflowRunner
from subprocess import call


class Runner(FluxWorkflowRunner):
    def __init__(self, query, output_dp, kmer_size, env_fp, max_seqs, hash_size, max_ppn, max_mem):
        self.base_fp = os.path.dirname(os.path.abspath(__file__))
        if not os.path.isdir(output_dp):
            os.makedirs(output_dp)

        self.input_files = self.check_input_files(query)
        if len(self.input_files) == 0:
            sys.exit('No usable files from --query')

        self.max_seqs = max_seqs
        if self.max_seqs == 0:
            sys.exit('Need max seq count')

        self.output_dp = output_dp
        self.kmer_size = kmer_size
        self.env_fp = env_fp
        self.hash_size = hash_size
        self.max_ppn = max_ppn
        self.max_mem = max_mem

    def check_input_files(self, query):
        files = []
        for qry in query:
            if os.path.isdir(qry):
                for filename in os.scandir(qry):
                    if filename.is_file():
                        files.append(filename.path)
            if os.path.isfile(qry):
                files.append(qry)
            else:
                warn('--query "{}" neither file nor directory'.format(qry))
        return files

    def workflow(self):
        """ method invoked on class instance run call """
        job_names = []

        # 1. Subset input files
        subset_dp = os.path.join(self.output_dp, 'subset')
        subset_cmd_fp = os.path.join(self.base_fp, 'autofizkin', 'fa_subset.py')
        subset_files = []
        for i, input_file in enumerate(self.input_files):
            outfile = os.path.join(subset_dp, os.path.basename('subset_{}'.format(input_file)))
            cmd = 'source {} && source activate py3 && '.format(self.env_fp)
            cmd += 'python {} -o {} -n {} -i fastq -t fastq {}'.format(subset_cmd_fp, outfile, self.max_seqs, input_file)
            self.addTask("fa_subset_{}".format(i), nCores=1, memMb=768, command=cmd)
            job_names.append("fa_subset_{}".format(i))


@click.command()
@click.argument('query', nargs=-1, type=click.Path(exists=True))
@click.option('--output_dp', '-o', default='autofizkin_output')
@click.option('--kmer_size', '-k', help='Kmer size', default=20)
@click.option('--max_seqs', '-x', help='Max num of seqs per input file', default=500000)
@click.option('--hash_size', '-s', help='Jellyfish hash size', default='100M')
@click.option('--flux/--no-flux', default=False)
@click.option('--account', '-a')
@click.option('--ppn', '-p', default=4)
@click.option('--mem', '-m', default='20000') # current limitation, only handles mb
@click.option('--walltime', '-w', default='2:00:00')
def main(query, output_dp, kmer_size, max_seqs, hash_size, flux, account, ppn, mem, walltime):
    print query
    if not os.path.exists(output_dp): os.makedirs(output_dp)
    fp = os.path.dirname(os.path.abspath(__file__)).replace('bin','')
    environment = os.path.join(fp, 'dependencies/miniconda/bin/activate')

    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(output_dp, 'logs/runner_{}'.format(r))
    runner = Runner(query=list(query), output_dp=output_dp, kmer_size=kmer_size, env_fp=environment,
                    max_seqs=max_seqs, hash_size=hash_size, max_ppn=ppn, max_mem=mem)
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ppn, memMb=mem) 

if __name__ == "__main__":
    main()
