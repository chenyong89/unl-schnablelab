import os.path as op
import sys
from schnablelab.apps.base import ActionDispatcher, OptionParser, create_slurm
from .base import FarmCPU_header, GAPIT_header, MVP_Data_header, MVP_Run_header

gemma = op.abspath(op.dirname(__file__)) + '/../apps/gemma'
tassel = op.abspath(op.dirname(__file__)) +\
      '/../apps/tassel-5-standalone/run_pipeline.pl'


def main():
    actions = (
        ('farmcpu', 'run GWAS using FarmCPU (muti-loci mixed model)'),
        ('gapit_cMLM', 'run GWAS using GAPIT (compressed Mixed Linear Model)'),
        ('gemma_GLM', 'run GWAS using GEMMA (general linear model)'),
        ('gemma_MLM', 'run GWAS using GEMMA (mixed linear model)'),
        ('mvp_data', 'prepare files for running MVP'),
        ('mvp', 'Run both MLM and FarmCPU models in MVP'),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def farmcpu(args):
    """
    %prog farmcpu phenotype_file genotype_data_prefix PCA_file

    Args:
        phenotype_file:
            phenotype file is tab delimited with header
        genotype_data_prefix:
            file prefix for GM and GD genotype files
            Chr must be numbers in the GM file
        PCA_file:
            PCA file
    """
    p = OptionParser(farmcpu.__doc__)
    p.add_slurm_opts(job_prefix='farmcpu')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())

    pheno, geno_prefix, PCA = args
    memo = '.'.join(pheno.split('/')[-1].split('.')[0:-1])
    r_fn = f'{memo}.FarmCPU.R'
    with open(r_fn, 'w') as f:
        farmcpu_cmd = FarmCPU_header.format(pheno=pheno,
                                            geno_prefix=geno_prefix,
                                            pca=PCA,
                                            memo=memo)
        f.write(farmcpu_cmd)
    print(f'R file to run on premises: {r_fn}')

    cmd_header = 'module load R/3.3\n'
    cmd = f'R CMD BATCH {r_fn}'
    slurm_dict = vars(opts)
    slurm_dict['cmd_header'] = cmd_header
    create_slurm([cmd], slurm_dict)


def gapit_cMLM(args):
    """
    %prog gapit_cMLM phenotype_file genotype_data_prefix PCA_file Kinship_file

    Run GAPIT compressed Mixed Linear Model

    Args:
        phenotype_csv: phenotype file is tab delimited with header
        genotype_data_prefix: prefix of GM and GD files
        PCA_file: input PCA file
        Kinship_file: input Kinship matrix file
    """
    p = OptionParser(gapit_cMLM.__doc__)
    p.add_slurm_opts(job_prefix='gapit_cMLM')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    pheno, geno_prefix, PCA, Kinship = args
    memo = '.'.join(pheno.split('.')[0:-1])
    r_fn = f'{memo}.cMLM.R'
    with open(r_fn, 'w') as f:
        gapit_cmd = GAPIT_header.format(pheno=pheno,
                                        geno_prefix=geno_prefix,
                                        pca=PCA,
                                        kinship=Kinship,
                                        memo=memo)
        f.write(gapit_cmd)
    print(f'R file to run on premises: {r_fn}')

    cmd_header = 'module load R/3.3\n'
    cmd = f'R CMD BATCH {r_fn}\n'
    slurm_dict = vars(opts)
    slurm_dict['cmd_header'] = cmd_header
    create_slurm([cmd], slurm_dict)


def gemma_GLM(args):
    """
    %prog gemma_GLM phenotype.csv genotype_data_prefix output_dir

    RUN GEMMA General Linear Model
    Args:
        genotype_data_prefix:
            filename prefix of GEMMA mean and annotation files
    """
    p = OptionParser(gemma_GLM.__doc__)
    p.add_slurm_opts(job_prefix='gemma_GLM')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    pheno, geno_prefix, output_dir = args
    meanG, annoG = geno_prefix+'.mean', geno_prefix+'.annotation'
    output_prefix = pheno.split('.')[0]
    cmd = f'{gemma} -g {meanG} -p {pheno} -a {annoG} -lm 4 -outdir '\
          f'{output_dir} -o {output_prefix}'
    print(f'command to run on premises:\n\t{cmd}')

    slurm_dict = vars(opts)
    create_slurm([cmd], slurm_dict)


def gemma_MLM(args):
    """
    %prog gemma_MLM phenotype.csv genotype_data_prefix output_dir

    RUN GEMMA Mixed Linear Model
    Args:
        genotype_data_prefix:
            filename prefix of GEMMA mean and annotation files
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
    print('command to run on premises:\n', cmd)

    slurm_dict = vars(opts)
    create_slurm([cmd], slurm_dict)


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
    r_fn = f'{output_prefix}.mvp.data.R'
    with open(r_fn, 'w') as f:
        f.write(MVP_Data_header.format(hapmap=hmp,
                                       output=output_prefix))
    print(f'R file to run on premises: {r_fn}')

    cmd_header = 'module load R\n'
    cmd = f'R CMD BATCH {r_fn}\n'
    slurm_dict = vars(opts)
    slurm_dict['cmd_header'] = cmd_header
    create_slurm([cmd], slurm_dict)


def mvp(args):
    """
    %prog mvp phenotype.csv GMKC_prefix

    Args:
        GMKC_prefix: the prefix of Genotype_Map_Kinship_Covariates files.
    """
    p = OptionParser(mvp.__doc__)
    p.add_slurm_opts(job_prefix='mvp')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    pheno, output_prefix, = args
    r_fn = f'{output_prefix}.mlm.farmcpu.R'
    with open(r_fn, 'w') as f:
        f.write(MVP_Run_header.format(pheno=pheno,
                                      output_prefix=output_prefix))
    print(f'R file to run on premises: {r_fn}')

    cmd_header = 'module load R\n'
    cmd = f'R CMD BATCH {r_fn}\n'
    slurm_dict = vars(opts)
    slurm_dict['cmd_header'] = cmd_header
    create_slurm([cmd], slurm_dict)


if __name__ == "__main__":
    main()
