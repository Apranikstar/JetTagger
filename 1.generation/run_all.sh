#!/bin/bash

base_config=$(cat << EOF
! 1) Settings that will be used in a main program.
Main:numberOfEvents = 1            ! number of events to generate
Main:timesAllowErrors = 100        ! abort run after this many flawed events
#Random:seed = 1234                ! initialize random generator with a seed

! 2) Settings related to output in init(), next() and stat() functions.
Init:showChangedSettings = on
Init:showAllSettings = off
Init:showChangedParticleData = on
Init:showAllParticleData = off
Next:numberCount = 10
Next:numberShowLHA = 1
Next:numberShowInfo = 1
Next:numberShowProcess = 1
Next:numberShowEvent = 1
Stat:showPartonLevel = off

! 3) Tell Pythia that LHEF input is used
Beams:frameType             = 4
Beams:setProductionScalesFromLHEF = off
EOF
)

mkdir -p cmdFiles
export EOS_MGM_URL="root://eospublic.cern.ch"
EOSDIR="/eos/user/h/hfatehi/generatedFiles/tt"

for lhefile in mg_pp_tt_HT_2000_100000_5f_84TeV/*.lhe; do
    base=$(basename "$lhefile" .lhe)
    cardfile="cmdFiles/card_${base}.cmd"

    echo "$base_config" > "$cardfile"
    echo "Beams:LHEF = $lhefile" >> "$cardfile"

    seed=$(( 100000 + RANDOM % 900000 ))
    echo "Random:seed = $seed" >> "$cardfile"
    echo "Main:numberOfEvents = 20000" >> "$cardfile"

    outputroot="events_${seed}.root"

    DelphesPythia8_EDM4HEP card.tcl edm4hep_output_config.tcl "$cardfile" "$outputroot"

    cp "$outputroot" "$EOSDIR"
    rm "$outputroot"
done



