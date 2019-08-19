#!/bin/bash

ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

source $ROOT_DIR/lib/common.sh

VERSION=""

# Arguments =============================================================================================

function print_help {
cat <<EOF

-v VERSION
--version VERSION           The package version. (required)

--help                      Show this help and exit.

EOF
}

ARGS=""
while [ $# -gt 0 ]; do
	case $1 in
		-v | --version ) shift; VERSION="$1" ;;

		-h | --help ) print_help; exit 0 ;;

		* )
			if [ -n "$ARGS" ]; then
				ARGS="$ARGS $1"
			else
				ARGS="$1"
			fi
			;;
	esac

	shift
done

if [ ! -n "$VERSION" ]; then
	log error "Version is required. Use --help for further information."
	exit -1
fi

# Package version ==================================================================================================

cat >$ROOT_DIR/lib/python/intogensm/version.py <<EOF
# This file is generated automatically by /distrib scripts
VERSION="$VERSION"
RELEASE_DATE="$(date +%d/%m/%Y)"
AUTHORS="Biomedical Genomics"
AUTHORS_EMAIL=""
EOF
