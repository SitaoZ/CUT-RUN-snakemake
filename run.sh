# snakemake --dag  | dot -Tsvg > dag.svg
nohup snakemake -s Snakefile --cores 50 --rerun-incomplete --latency-wait 100 &
