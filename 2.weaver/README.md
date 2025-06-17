### Running Weaver on the Samples

#### Stage 1
```bash
fccanalysis run stage1.py --output output.root --files-list input.root --ncpus 16
```
#### Stage 2
```bash
python stage2.py output.root out.root 0 100
```

### stage_all
```bash
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ --outdir ./output/  --sample mgp8_pp_tt_HT_2000_100000_5f_84TeV --ncpus 16 
```
