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
