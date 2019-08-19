#!/bin/bash

ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

source $ROOT_DIR/lib/common.sh

DISTRIB_DIR="$ROOT_DIR/distrib"

OVERWRITE="no"

VERSION=""

# Arguments =============================================================================================

function print_help {
cat <<EOF

-v VERSION
--version VERSION           The package version. (required)

--overwrite                 If package file already exists then overwrite otherwise do not regenerate.

--help                      Show this help and exit.

EOF
}

ARGS=""
while [ $# -gt 0 ]; do
	case $1 in
		-v | --version ) shift; VERSION="$1" ;;

		--overwrite ) OVERWRITE="yes" ;;
		
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

# Package version ========================================================================================

$DISTRIB_DIR/update-version.sh -v $VERSION

# Create package =========================================================================================

PACKAGE_NAME="intogen-sm-$VERSION"
PACKAGE_PATH="$DISTRIB_DIR/$PACKAGE_NAME.tar.gz"

if [ ! -e $PACKAGE_PATH -o "$OVERWRITE" == "yes" ]; then
	log info "Creating package $PACKAGE_NAME ..."

	TMP_DIR=$(mktemp -d)
	
	ln -s $ROOT_DIR $TMP_DIR/$PACKAGE_NAME
	tar -cvhzf $PACKAGE_PATH -X $DISTRIB_DIR/exclude-tar.txt --exclude="$PACKAGE_NAME/runtime" -C $TMP_DIR $PACKAGE_NAME
	rm -rf $TMP_DIR
else
	log warn "Package already exists at $PACKAGE_PATH. Not overwritting unless --overwrite specified"
fi
