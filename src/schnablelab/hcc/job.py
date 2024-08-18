"""
create, submit, canceal jobs
"""

import sys
import subprocess
from pathlib import Path
from subprocess import call
from subprocess import Popen
from schnablelab.apps.base import ActionDispatcher, OptionParser


def main():
    actions = (
        ('submit', 'submit a batch of jobs or all of them'),
        ('quickjob', 'create a quick slurm job'),
        ('cancel', 'canceal running, pending or all jobs'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def submit(args):
    """
    %prog submit dir

    Submit jobs under the dir
    """
    p = OptionParser(submit.__doc__)
    p.add_option("--pattern", default="*.slurm",
                 help="specify the filename pattern of slurm files, "
                 "remember to add quotes [default: %default]")
    p.add_option("--partition", default='jclarke',
                 choices=('batch', 'jclarke', 'gpu', 'schnablelab'),
                 help="choose which partition you are going to submit")
    p.add_option("--range", default='all',
                 help="how many jobs you gonna submit. "
                 "e.g., '1-10', '11-20', 'all'. use 1-based coordinate")
    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())

    folder, = args
    partition = '-p %s' % opts.partition
    alljobs = ['sbatch %s %s' % (partition, i)
               for i in Path(folder).glob(opts.pattern)]
    print("Total %s jobs under '%s'" % (len(alljobs), folder))

    if opts.range == 'all':
        for i in alljobs:
            print(i)
            call(i, shell=True)
    else:
        start = int(opts.range.split('-')[0])
        end = int(opts.range.split('-')[1])
        if end <= len(alljobs):
            for i in alljobs[start - 1: end]:
                print(i)
                call(i, shell=True)
            print(f'{len(alljobs[start - 1: end])} of total {len(alljobs)}'
                  f' were submitted. [{start} to {end}] this time.')
        else:
            print('jobs exceed the limit')


def cancel(args):
    """
    %prog cancel

    cancel running jobs on HCC
    """
    p = OptionParser(cancel.__doc__)
    p.add_option("--status", default='running', choices=('running', 'pending'),
                 help="specify the status of the jobs you want to cancel")
    p.add_option("--partition", default='jclarke',
                 choices=('gpu', 'batch', 'jclarke'),
                 help="specify the partition where jobs are runnig")
    opts, args = p.parse_args(args)
    if len(args) != 0:
        sys.exit(not p.print_help())
    myjobs = Popen('squeue -u cmiao', shell=True,
                   stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    running_jobs, pending_jobs, others = [], [], []
    for i in myjobs.split('\n'):
        j = i.strip().split()
        if j[4] == 'R':
            running_jobs.append(j[0])
        elif j[4] == 'PD':
            pending_jobs.append(j[0])
        else:
            others.append(j[0])
    cmd = 'scancel %s' % (' '.join(running_jobs)) \
        if opts.status == 'running' \
        else 'scancel %s' % (' '.join(pending_jobs))

    print(cmd)


quick_job_header = '''#!/bin/bash
#SBATCH --time=120:00:00
#SBATCH --mem-per-cpu=10000
#SBATCH --job-name={jn}
#SBATCH --error={jn}.err
#SBATCH --output={jn}.out

'''


def quickjob(args):
    """
    %prog quickjob job_name cmd

    generate a qucik slurm job for the cmd
    """
    if len(args) == 0:
        sys.exit('specify job_name and command')
    jobname, *cmd = args
    cmd = ' '.join(cmd)
    header = quick_job_header.format(jn=jobname)
    header += cmd
    with open(f'{jobname}.slurm', 'w') as f:
        f.write(header)
    print('slurm file {jobname}.slurm has been created...')


if __name__ == "__main__":
    main()
