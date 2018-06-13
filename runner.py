#!/usr/bin/env python

import os, sys
import click
import random
import string

from metago.pyflux import FluxWorkflowRunner
from metago.quality_control import RunQualityControl
from metago.quality_control import SampleQualityControl
from subprocess import call


class Runner(FluxWorkflowRunner):
    def __init__(query, output_dp, kmer_size, max_seqs, hash_size, max_ppn, max_mem):
        if not os.path.isdir(output_dp):
            os.makedirs(output_dp)

        self.input_files = self.check_input_files(self.query)
        if len(self.input_files) == 0:
            sys.exit('No usable files from --query')

        self.max_seqs = max_seqs
        if len(self.max_seqs == 0:
            sys.exit('Need max seq count')

        self.output_dp = output_dp
        self.kmer_size = kmer_size
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
            if os.path.isfile(query):
                files.append(qry)
            else:
                warn('--query "{}" neither file nor directory'.format(qry))
        return files

    def workflow(self):
        """ method invoked on class instance run call """
        job_names = []

        # 1. Subset input files
        subset_dp = os.path.join(self.output_dp, 'subset')
        subset_files = []
        for i, input_file in enumerate(self.input_files):
            outfile = os.path.join(subset_dp, os.path.basename('subset_{}'.format(input_file)))
            cmd = 'python fa_subset -o {} -n {} -i fastq -t fastq {}'.format(outfile, self.max_seqs, input_file)
            self.addTask("fa_subset_{}".format(i), nCores=1, memMb=768, command=self.cmd)
            job_names.append("fa_subset_{}".format(i))


@click.group()
@click.option('--query', '-q', help='Input fastq/fasta files', multiple=True)
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
    if not os.path.exists(output_dp): os.makedirs(output_dp)
    fp = os.path.dirname(os.path.abspath(__file__)).replace('bin','')
    environment = os.path.join(fp, 'dependencies/miniconda/bin/activate')

    if flux:
        if not account: sys.exit('To attempt a submission to the flux cluster you need to supply an --account/-a')
        full_dp = os.path.dirname(os.path.abspath(__file__))
        metago_fp = os.path.join(full_dp, 'metago')
        metago_index = -1
        for i, s in enumerate(sys.argv):
            if 'metago' in s:
                metago_index = i
                break
        cmd = sys.argv
        cmd[metago_index] = metago_fp
        cmd.remove('--flux')
        cmd = ' '.join(cmd)
        qsub = 'qsub -N metaGO -A {} -q fluxod -l nodes=1:ppn={}:largemem,mem={}mb,walltime={}'.format(
                                                                   account, ppn, mem, walltime)
        call('echo "source {} && python {}" | {}'.format(environment, cmd, qsub), shell=True)
        sys.exit('Launched command via Flux')

    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(output_dp, 'logs/runner_{}'.format(r))
    runner = Runner(query=query, output_dp=output_dp, kmer_size=kmer_size,
                    max_seqs=max_seqs, hash_size=hash_size, max_ppn=ppn, max_mem=mem)
    workflow_runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ppn, memMb=mem) 

if __name__ == "__main__":
    main()
