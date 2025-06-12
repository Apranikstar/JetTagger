### Running Weaver on the Samples

#### Stage 1
```bash
fccanalysis run stage1.py --output output.root --files-list input.root --ncpus 16
```
#### Stage 2
```bash
python stage2.py output.root out.root 0 100
```
