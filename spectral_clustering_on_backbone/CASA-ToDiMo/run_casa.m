delete(gcp('nocreate'))

c = parcluster('local');
c.NumWorkers = str2double(getenv('SLURM_CPUS_PER_TASK'));

parpool(c)

classAverageDisordered

