#!/bin/bash

PRODUCT_DIR="$VALIDATE_DRP_DIR"
VALIDATION_DATA_DIR="$VALIDATION_DATA_DECAM_DIR/instcal"

CAMERA=Decam
CONFIG_FILE="${PRODUCT_DIR}/config/decamConfig.py"
MAPPER=lsst.obs.decam.DecamMapper

"${PRODUCT_DIR}/examples/processData.sh" \
    -c "$CAMERA" \
    -m "$MAPPER" \
    -v "$VALIDATION_DATA_DIR" \
    -f "$CONFIG_FILE" \
    -i ingestImagesDecam.py \
    -- "$@"
"${PRODUCT_DIR}/examples/processData.sh" \
    -c "$CAMERA" \
    -- "$@"
