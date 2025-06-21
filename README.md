# JetTagger

Quick guide to preparing data for training:

```bash
source /cvmfs/fcc.cern.ch/sw/latest/setup.sh
git clone https://github.com/Apranikstar/JetTagger.git
cd 2.weaver
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ --outdir ./output/  --sample mgp8_pp_tt_HT_2000_100000_5f_84TeV --ncpus 16
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ --outdir ./output/  --sample mgp8_pp_jj_HT_2000_100000_5f_84TeV --ncpus 16
python stage_plots.py --indir ./output --outdir ./plots/
```


Submitting job to HTCondor

run_tagger.sh (replace <initials>/<user>/... with your PATH
```bash
#!/bin/bash
# File: run_tagger.sh

# Environment setup
source /cvmfs/fcc.cern.ch/sw/latest/setup.sh

# Clone your repository
git clone https://github.com/Apranikstar/JetTagger.git
cd JetTagger/2.weaver

# Run your stages
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ --outdir /afs/cern.ch/work/<initials>/<user>/.../output/ --sample mgp8_pp_tt_HT_2000_100000_5f_84TeV --ncpus 16
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ --outdir /afs/cern.ch/work/<initials>/<user>/.../output/ --sample mgp8_pp_jj_HT_2000_100000_5f_84TeV --ncpus 16

# Generate plots
python stage_plots.py --indir /afs/cern.ch/work/h/hfatehi/JET/output --outdir /afs/cern.ch/work/h/hfatehi/JET/plots/

```
```bash
chmod +x run_tagger.sh
```

Next step:
```vim run_tagger.sub```
With this content:
```ini
executable = run_tagger.sh
output     = logs/run_tagger.out
error      = logs/run_tagger.err
log        = logs/run_tagger.log
request_cpus = 16
request_memory = 4GB
+JobFlavour = "longlunch"
queue
```

Create logs directory:
```mkdir -p logs```

Submit the job:
```condor_submit run_tagger.sub```



