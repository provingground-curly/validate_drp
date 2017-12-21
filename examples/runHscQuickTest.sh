#!/bin/bash

PRODUCT_DIR="$VALIDATE_DRP_DIR"
VALIDATION_DATA_DIR="$CI_HSC_DIR/raw"
CALIB_DATA="$CI_HSC_DIR/CALIB"

CAMERA=HscQuick
CONFIG_FILE="${PRODUCT_DIR}/config/hscConfig.py"
MAPPER=lsst.obs.hsc.HscMapper

echo "This will not work because ci_hsc has gone over to using new style reference catalogs"
echo "Please remove this when DM-13298 is addressed"
exit 1

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
