#!/bin/bash

PRODUCT_DIR="$VALIDATE_DRP_DIR"
VALIDATION_DATA_DIR="$VALIDATION_DATA_HSC_DIR/raw"
CALIB_DATA="$VALIDATION_DATA_HSC_DIR/CALIB"

CAMERA=Hsc
CONFIG_FILE="${PRODUCT_DIR}/config/hscConfig.py"
MAPPER=lsst.obs.hsc.HscMapper
ASTROMDIR=sdss-dr9-fink-v5b

print_error() {
    >&2 echo "$@"
}

DOPROCESS=true
DOVERIFY=true

usage() {
    print_error
    print_error "Usage: $0 [-pv] [-h] [-- <extra options to validateDrp.py>]"
    print_error
    print_error "Specifc options:"
    print_error "   -p          Skip processing?"
    print_error "   -v          Skip verification?"
    print_error "   -h          show this message"
    exit 1
}

# thank OSX for not including getopt
while getopts "vp" option; do
    case "$option" in
        p)  DOPROCESS=false;;
        v)  DOVERIFY=false;;
        h)  usage;;
        *)  usage;;
    esac
done
shift $((OPTIND-1))

if [ "$DOPROCESS" = true ] ; then
    "${PRODUCT_DIR}/examples/processData.sh" \
        -c "$CAMERA" \
        -m "$MAPPER" \
        -v "$VALIDATION_DATA_DIR" \
        -f "$CONFIG_FILE" \
        -e "fits" \
        -a "$ASTROMDIR" \
        -d "$CALIB_DATA" \
        -r
fi

if [ "$DOVERIFY" = true ] ; then
    "${PRODUCT_DIR}/examples/validateRepo.sh" \
        -c "$CAMERA" \
        -- "$@"
fi
