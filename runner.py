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
        if not os.path.exists(subset_dp): os.makedirs(subset_dp)

        subset_cmd_fp = os.path.join(self.base_fp, 'autofizkin', 'fa_subset.py')
        subset_files = []
        for i, input_file in enumerate(self.input_files):
            outdir = os.path.join(subset_dp, os.path.basename(input_file).split('.')[0])
            outfile = os.path.join(outdir, input_file)
            subset_files.append(outfile)
            if os.path.exists(outfile): continue
            cmd = 'source {} && source activate py3 && '.format(self.env_fp)
            cmd += 'python {} -o {} -n {} -i fastq -t fastq {}'.format(subset_cmd_fp, outdir, self.max_seqs, input_file)
            self.addTask("fa_subset_{}".format(i), nCores=1, memMb=768, command=cmd)
            job_names.append("fa_subset_{}".format(i))

        # 2. Count kmers
        kmer_dp = os.path.join(self.output_dp, 'kmer_counts')
        if not os.path.exists(kmer_dp): os.makedirs(kmer_dp)

        for i, subset_fp in enumerate(subset_files):
            kmer_fp = os.path.join(kmer_dp, os.path.basename(subset_fp).split('.')[0] + '.jf')
            if os.path.exists(kmer_fp): continue
            cmd = 'source {} && source activate py3 && '.format(self.env_fp)
            cmd += 'jellyfish count -m {} -t {} -s {} '.format(self.kmer_size, self.max_ppn, self.hash_size)
            cmd += '-o {} {}'.format(kmer_fp, subset_fp)
            self.addTask('kmer_count_{}'.format(i), nCores=self.max_ppn, memMb=self.max_mem,
                                                          command=cmd, dependencies=job_names)
            job_names.append('kmer_count_{}'.format(i))

        # 3. kmer pairwise compare
        keep_dp = os.path.join(self.output_dp, 'reads_kept')
        if not os.path.exists(keep_dp): os.makedirs(keep_dp)
        reject_dp = os.path.join(self.output_dp, 'reads_rejected')
        if not os.path.exists(reject_dp): os.makedirs(reject_dp)

        # kmer_files check needs to occur in scheduled process as they won't be there in the beginning
        kmer_files = [os.path.join(kmer_dp, f) for f in os.listdir(kmer_dp) if os.path.isfile(os.path.join(kmer_dp, f))]
        if not kmer_files: sys.exit('No kmer files to compare!')

        """ TODO - rewrite query_per_sequence or only run this in a container """
        for i, kmer_fp in enumerate(kmer_files):
            index_name = os.path.basename(kmer_fp)
            keep = os.path.join(keep_dp, index_name)
            if not os.path.exists(keep): os.makedirs(keep)
            reject = os.path.join(reject_dp, index_name)
            if not os.path.exists(reject): os.makedirs(reject)

            for k, subset_fp in enumerate(subset_files):
                query_name = os.path.basename(subset_fp)
                keep_fp = os.path.join(keep, query_name)
                reject_fp = os.path.join(reject, query_name)
                if os.path.exists(keep_fp) and os.path.exists(reject_fp): continue
                cmd = 'source {} && source activate py3 && '.format(self.env_fp)
                cmd += 'query_per_sequence 1 10 {} {} 1>{} 2>{}'.format(kmer_fp, subset_fp, keep_fp, reject_fp)
                self.addTask('kmdr_count_{}'.format(i*2+k), nCores=self.max_ppn, memMb=self.max_mem,
                                                                  command=cmd, dependencies=job_names)
                

            
        


@click.command()
@click.argument('query', nargs=-1, type=click.Path(exists=True))
@click.option('--reference', '-r', default=None, help='reference to use for removing unwanted reads')
@click.option('--output_dp', '-o', default='autofizkin_output', help='directory path to save output to')
@click.option('--kmer_size', '-k', default=20, help='size of kmer to cluster by')
@click.option('--max_seqs', '-x', help='Max num of seqs per input file', default=500000)
@click.option('--hash_size', '-s', help='Jellyfish hash size', default='100M')
@click.option('--flux/--no-flux', default=False, help='whether to attempt to run on UofM flux cluster')
@click.option('--account', '-a', help='name of UofM Flux account to use')
@click.option('--ppn', '-p', default=4, help='number of processors to request from cluster')
@click.option('--mem', '-m', default='20000', help='amount of memory (in MB) to request from cluster')
@click.option('--walltime', '-w', default='2:00:00', help='amount of time to request for job to run (HH:MM:SS)')
def main(query, reference, output_dp, kmer_size, max_seqs, hash_size, flux, account, ppn, mem, walltime):
    if not os.path.exists(output_dp): os.makedirs(output_dp)
    fp = os.path.dirname(os.path.abspath(__file__)).replace('bin','')
    environment = os.path.join(fp, 'dependencies/miniconda/bin/activate')

    if flux:
        if not account: sys.exit('To attempt a submission to the flux cluster you need to supply an --account/-a')
        full_dp = os.path.dirname(os.path.abspath(__file__))
        cmd_fp = os.path.join(full_dp, sys.argv[0])

        args = sys.argv
        args.remove('--flux')
        args.remove(sys.argv[0])
        args = ' '.join(args)
        qsub = 'qsub -N {} -A {} -q fluxod -l nodes=1:ppn={}:largemem,mem={}mb,walltime={}'.format(
                                    cmd_fp.split('/')[-1].split('.')[0], account, ppn, mem, walltime)
        call('echo "source {} && python {} {}" | {}'.format(environment, cmd_fp, args, qsub), shell=True)
        sys.exit('Launched via Flux')

    r = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
    log_output_dp = os.path.join(output_dp, 'logs/runner_{}'.format(r))
    runner = Runner(query=list(query), output_dp=output_dp, kmer_size=kmer_size, env_fp=environment,
                    max_seqs=max_seqs, hash_size=hash_size, max_ppn=ppn, max_mem=mem)
    runner.run(mode='local', dataDirRoot=log_output_dp, nCores=ppn, memMb=mem) 

if __name__ == "__main__":
    main()
