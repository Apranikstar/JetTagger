### How to run weaver on the samples
##stage1:
```python
fccanalysis run stage1.py --output output.root --files-list input.root --ncpus 16
```
##stage2:
```python
python stage2.py  output.root out.root 0 100
```
