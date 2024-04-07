import os.path as op
import sys
from schnablelab.apps.base import SLURM_header, ActionDispatcher, OptionParser
from .base import FarmCPU_header, GAPIT_header, MVP_Data_header, MVP_Run_header

gemma = op.abspath(op.dirname(__file__)) + '/../apps/gemma'
tassel = op.abspath(op.dirname(__file__)) +\
      '/../apps/tassel-5-standalone/run_pipeline.pl'


def main():
    actions = (
        ('farmcpu', 'run GWAS using FarmCPU (muti-loci mixed model)'),
        ('gapit_cMLM', 'run GWAS using GAPIT compressed Mixed Linear Model'),
        ('gemma_GLM', 'run GWAS using GEMMA general linear model'),
        ('gemma_MLM', 'run GWAS using GEMMA mixed linear model '),
        ('mvp_data', 'prepare files for running MVP'),
        ('mvp', 'Run both MLM and FarmCPU models in MVP'),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def farmcpu(args):
    """
    %prog farmcpu phenotype.csv genotype_data_prefix PCA_file

    phenotype.csv with header, tab delimited
    genotype_data_prefix for GM and GD files
    Chr must be numbers in the GM file
    """
    p = OptionParser(farmcpu.__doc__)
    p.add_slurm_opts(job_prefix='farmcpu')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())

    pheno, geno_prefix, PCA = args
    mem = '.'.join(pheno.split('/')[-1].split('.')[0:-1])
    with open(f'{mem}.FarmCPU.R', 'w') as f:
        farmcpu_cmd = FarmCPU_header.format(pheno, geno_prefix, geno_prefix,
                                            PCA, mem)
        f.write(farmcpu_cmd)

    with open('{mem}.FarmCPU.slurm', 'w') as f:
        h = SLURM_header
        h += 'module load R/3.3\n'
        header = h.format(opts.time, opts.memory, opts.prefix, opts.prefix,
                          opts.prefix)
        f.write(header)
        cmd = 'R CMD BATCH %s.FarmCPU.R' % mem
        f.write(cmd)
    print(f'R script {mem}.FarmCPU.R and slurm file {mem}.FarmCPU.slurm has '
          'been created. You can sbatch your job file now.')


def gapit_cMLM(args):
    """
    %prog gapit_cMLM phenotype.csv genotype_data_prefix PCA_file Kinship_file

    Run GAPIT compressed Mixed Linear Model
    phenotype.csv is tab delimited with header
    genotype_data_prefix for GM and GD files
    """
    p = OptionParser(gapit_cMLM.__doc__)
    p.add_slurm_opts(job_prefix='gapit_cMLM')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    pheno, geno_prefix, PCA, Kinship = args
    mem = '.'.join(pheno.split('.')[0:-1])
    with open('%s.cMLM.R' % mem, 'w') as f:
        gapit_cmd = GAPIT_header.format(pheno, geno_prefix, geno_prefix, PCA,
                                        Kinship, mem)
        f.write(gapit_cmd)

    with open('%s.cMLM.slurm' % mem, 'w') as f:
        h = SLURM_header
        h += 'module load R/3.3\n'
        header = h.foramt(opts.time, opts.memory, opts.prefix, opts.prefix,
                          opts.prefix)
        f.write(header)
        cmd = 'R CMD BATCH %s.cMLM.R\n' % mem
        f.write(cmd)
    print(f'R script {mem}.cMLM.R and slurm file {mem}.cMLM.slurm has been '
          'created. You can sbatch your job file now.')


def gemma_GLM(args):
    """
    %prog gemma_GLM phenotype.csv genotype_data_prefix output_dir

    RUN GEMMA General Linear Model
    genotype_data_prefix of GEMMA mean and annotation files
    """
    p = OptionParser(gemma_GLM.__doc__)
    p.add_slurm_opts(job_prefix='gemma_GLM')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    pheno, geno_prefix, output_dir = args
    meanG, annoG = geno_prefix+'.mean', geno_prefix+'.annotation'
    output_prefix = pheno.split('.')[0]
    cmd = '%s -g %s -p %s -a %s -lm 4 -outdir %s -o %s'\
        % (gemma, meanG, pheno, annoG, output_dir, output_prefix)
    print('Command running on the local node:\n', cmd)

    h = SLURM_header
    header = h.format(opts.time, opts.memory, opts.prefix, opts.prefix,
                      opts.prefix)
    header += cmd
    with open('%s.glm.slurm' % output_prefix, 'w') as f:
        f.write(header)
    print(f'slurm file {output_prefix}.glm.slurm has been created.'
          'you can sbatch your job now.')


def gemma_MLM(args):
    """
    %prog gemma_MLM phenotype.csv genotype_data_prefix output_dir

    RUN GEMMA Mixed Linear Model
    genotype_data_prefix of GEMMA mean and annotation files
    """
    p = OptionParser(gemma_MLM.__doc__)
    p.add_option('--kinship', default=False,
                 help='specify the relatedness matrix file name')
    p.add_option('--pca', default=False,
                 help='specify the principle components file name')
    p.add_slurm_opts(job_prefix='gemma_MLM')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    pheno, geno_prefix, output_dir = args
    meanG, annoG = geno_prefix+'.mean', geno_prefix+'.annotation'
    output_prefix = '.'.join(pheno.split('/')[-1].split('.')[0:-1])
    cmd = '%s -g %s -p %s -a %s -lmm 4 -outdir %s -o %s'\
        % (gemma, meanG, pheno, annoG, output_dir, output_prefix)
    if opts.kinship:
        cmd += ' -k %s' % opts.kinship
    if opts.pca:
        cmd += ' -c %s' % opts.pca
    print('The command running on the local node:\n', cmd)

    h = SLURM_header
    header = h.format(opts.time, opts.memory, opts.prefix, opts.prefix,
                      opts.prefix)
    header += cmd
    with open('%s.mlm.slurm' % output_prefix, 'w') as f:
        f.write(header)
    print(f'slurm file {output_prefix}.mlm.slurm has been created.'
          'you can sbatch your job file.')


def mvp_data(args):
    """
    %prog mvp_data hmp out_prefix

    Prepare the genotype, map, kinship, and covariants for running MVP
    """
    p = OptionParser(mvp_data.__doc__)
    p.add_slurm_opts(job_prefix='mvp_data')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    hmp, output_prefix, = args
    with open(f'{output_prefix}.mvp.data.R', 'w') as f:
        f.write(MVP_Data_header % (hmp, op))

    with open(f'{output_prefix}.mvp.data.slurm', 'w') as f:
        header = SLURM_header.format(opts.time, opts.memory, output_prefix,
                                     output_prefix, output_prefix)
        header += 'module load R\n'
        header += f'R CMD BATCH {output_prefix}.mvp.data.R\n'
        f.write(header)
    print(f'{output_prefix}.mvp.data.R and {output_prefix}.mvp.data.slurm '
          'have been created. You can submit slurm job now.')


def mvp(args):
    """
    %prog mvp phenotype.csv GMKC_prefix

    GMKC_prefix: the prefix of Genotype_Map_Kinship_Covariates files.
    """
    p = OptionParser(mvp.__doc__)
    p.add_slurm_opts(job_prefix='mvp')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    pheno, output_prefix, = args
    with open(f'{output_prefix}.mlm.farmcpu.R', 'w') as f:
        f.write(MVP_Run_header.format(pheno, output_prefix, output_prefix,
                                      output_prefix, output_prefix))
    with open(f'{output_prefix}.mlm.farmcpu.slurm', 'w') as f:
        header = SLURM_header.format(opts.time, opts.memory, output_prefix,
                                     output_prefix, output_prefix)
        header += 'module load R\n'
        header += f'R CMD BATCH {output_prefix}.mlm.farmcpu.R\n' % opts.prefix
        f.write(header)

    print(f'{output_prefix}.mlm.farmcpu.R and '
          f'{output_prefix}.mlm.farmcpu.slurm have been created.'
          'You can submit slurm job now.')


if __name__ == "__main__":
    main()
