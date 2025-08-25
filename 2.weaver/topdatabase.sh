#!/bin/bash

# Create logs directory
mkdir -p logs

#tops
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_tt_HT_2000_100000_5f_84TeV \
 --sample mgp8_pp_tt_HT_2000_100000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_tt_HT_2000_100000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_tt_Q_3000_10000_5f_84TeV \
 --sample mgp8_pp_tt_Q_3000_10000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_tt_Q_3000_10000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_tt_Q_10000_30000_5f_84TeV \
 --sample mgp8_pp_tt_Q_10000_30000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_tt_Q_10000_30000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_tt_Q_30000_84000_5f_84TeV \
 --sample mgp8_pp_tt_Q_30000_84000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_tt_Q_30000_84000_5f_84TeV.log 2>&1

#bb
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_bb_Q_3000_10000_5f_84TeV \
 --sample mgp8_pp_bb_Q_3000_10000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_bb_Q_3000_10000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_bb_Q_10000_30000_5f_84TeV \
 --sample mgp8_pp_bb_Q_10000_30000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_bb_Q_10000_30000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_bb_Q_30000_84000_5f_84TeV \
 --sample mgp8_pp_bb_Q_30000_84000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_bb_Q_30000_84000_5f_84TeV.log 2>&1

#cc
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_cc_Q_3000_10000_5f_84TeV \
 --sample mgp8_pp_cc_Q_3000_10000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_cc_Q_3000_10000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_cc_Q_10000_30000_5f_84TeV \
 --sample mgp8_pp_cc_Q_10000_30000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_cc_Q_10000_30000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_cc_Q_30000_84000_5f_84TeV \
 --sample mgp8_pp_cc_Q_30000_84000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_cc_Q_30000_84000_5f_84TeV.log 2>&1

#uuddggss
python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_uuddssgg_Q_3000_10000_5f_84TeV \
 --sample mgp8_pp_uuddssgg_Q_3000_10000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_uuddssgg_Q_3000_10000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_uuddssgg_Q_10000_30000_5f_84TeV \
 --sample mgp8_pp_uuddssgg_Q_10000_30000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_uuddssgg_Q_10000_30000_5f_84TeV.log 2>&1

python stage_all.py --indir /eos/experiment/fcc/hh/generation/DelphesEvents/fcc_v07/II/ \
 --outdir /eos/user/h/hfatehi/topTaggingGeneral/mgp8_pp_uuddssgg_Q_30000_84000_5f_84TeV \
 --sample mgp8_pp_uuddssgg_Q_30000_84000_5f_84TeV --ncpus 16 \
 > logs/mgp8_pp_uuddssgg_Q_30000_84000_5f_84TeV.log 2>&1
