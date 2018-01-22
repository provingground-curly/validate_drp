#!/bin/bash

PRODUCT_DIR="$VALIDATE_DRP_DIR"
VALIDATION_DATA_DIR="$VALIDATION_DATA_HSC_DIR/raw"
CALIB_DATA="$VALIDATION_DATA_HSC_DIR/CALIB"

CAMERA=Hsc
CONFIG_FILE="${PRODUCT_DIR}/config/hscConfig.py"
MAPPER=lsst.obs.hsc.HscMapper
ASTROMDIR=sdss-dr9-fink-v5b

"${PRODUCT_DIR}/examples/processData.sh" \
    -c "$CAMERA" \
    -m "$MAPPER" \
    -v "$VALIDATION_DATA_DIR" \
    -f "$CONFIG_FILE" \
    -e "fits" \
    -a "$ASTROMDIR" \
    -d "$CALIB_DATA" \
    -r \
    -- "$@"
"${PRODUCT_DIR}/examples/validateRepo.sh" \
    -c "$CAMERA" \
    -- "$@"
