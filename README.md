# Nextflow RNA-seq Pipeline (rna_nf)

## Conda environment
To create the 'rna' conda environment for Nextflow
```
conda env create -f environment.yml -p [path_and_environment name]
```
In the environment.yml, the environment name, conda channels and packages with or without restriction of versions.

Include this environment path in ~/.condarc, that is
```
envs_dirs:
  - /hpc/users/liuy22/LOAD/rna/py2_conda_env
```
Make sure conda is loaded via module load ananconda or miniconda in HPC, or locally available.

## Conda integration in Nextflow
In the Nextflow configuration, include 'conda' in the 'process' scope as the following. This ablates specification in 'processes' of the Nextflow script.
```
process {
	conda = '/home/yiyuan/Data/Nextflow/py2_conda_env/rna'
}
```

## Nextflow configuration file (*.config)

## Parameter file (*.yml)

## Local usage
As in the call_nf.sh, to run the RNA-seq pipeline locally,
```
nextflow rna.nf -c rna_local.config -params-file project_local.yml
```
or with more complexity,

```
nextflow rna.nf -c rna_local.config -params-file project_local.yml -with-trace -with-timeline -with-dag flowchart.html
```

## HPC usage
To be continued.
