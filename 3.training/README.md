#### How to run the model
```bash
 torchrun --standalone --nnodes=1 \
 --nproc_per_node=1 -m weaver.train \
 --data-train /afs/cern.ch/work/h/hfatehi/weaverTraining/database/stage2_*.root \
 --data-config /afs/cern.ch/work/h/hfatehi/weaverTraining/topTagger.yaml \
 --network-config /afs/cern.ch/work/h/hfatehi/weaverTraining/example_ParticleTransformer.py \
 --model-prefix /afs/cern.ch/work/h/hfatehi/weaverTraining/trainings \
 --num-workers 0 --gpus 0 --batch-size 16 \
 --start-lr 1e-3 --num-epochs 20 \
 --optimizer ranger --fetch-step 0.01 \
 --backend nccl
```


| Flavor | Branch Name     | Class Index |
| ------ | --------------- | ----------- |
| g      | `recojet_isG`   | 0           |
| q      | `recojet_isQ`   | 1           |
| s      | `recojet_isS`   | 2           |
| c      | `recojet_isC`   | 3           |
| b      | `recojet_isB`   | 4           |
| t      | `recojet_isT`   | 5           |
| j      | `recojet_isJ`   | 6           |
| tau    | `recojet_isTAU` | 7           |
