#!/bin/bash
# --------------------------------------------------------------------
# Single-local Chapel
# --------------------------------------------------------------------
export CHPL_HOME=~/chapel-1.24.0
CHPL_BIN_SUBDIR=`"$CHPL_HOME"/util/chplenv/chpl_bin_subdir.py`
export PATH="$PATH":"$CHPL_HOME/bin/$CHPL_BIN_SUBDIR"
export MANPATH="$MANPATH":"$CHPL_HOME"/man
# --------------------------------------------------------------------
# path for modules as libraries
# --------------------------------------------------------------------
export CHPL_MODULE_PATH=~/modules
# --------------------------------------------------------------------
# use all cores available
# --------------------------------------------------------------------
export CHPL_RT_NUM_THREADS_PER_LOCALE=MAX_LOGICAL
