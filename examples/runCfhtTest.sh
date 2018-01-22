#!/bin/bash

PRODUCT_DIR="$VALIDATE_DRP_DIR"
VALIDATION_DATA_DIR="$VALIDATION_DATA_CFHT_DIR/raw"

CAMERA=Cfht
CONFIG_FILE="${PRODUCT_DIR}/config/cfhtConfig.py"
MAPPER=lsst.obs.cfht.MegacamMapper

"${PRODUCT_DIR}/examples/processData.sh" \
    -c "$CAMERA" \
    -m "$MAPPER" \
    -v "$VALIDATION_DATA_DIR" \
    -f "$CONFIG_FILE" \
    -- "$@"
"${PRODUCT_DIR}/examples/validateRepo.sh" \
    -c "$CAMERA" \
    -- "$@"
