#!/bin/bash

cutadapt -j 0 -g XGATAACTGATCATAATCAGCCATACC --overlap 10 --discard-untrimmed -m 50 -e 0.2 -o $1 $2