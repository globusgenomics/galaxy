#!/bin/bash

ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

source $ROOT_DIR/lib/common.sh

DISTRIB_DIR="$ROOT_DIR/distrib"

TUNNEL="no"
TUNNEL_USER=$USER
TUNNEL_HOST=tunnel-host
TUNNEL_LOCAL_PORT=2022
TUNNEL_REMOTE_PORT=22
TUNNEL_REMOTE_HOST=remote-host

UPLOAD_PORT=2022
UPLOAD_HOST=localhost
UPLOAD_PATH=~
UPLOAD_USER=$USER

[ -f $DISTRIB_DIR/upload-defaults.sh ] && source $DISTRIB_DIR/upload-defaults.sh

VERSION=""

# Arguments =============================================================================================

function print_help {
cat <<EOF

-v VERSION
--version VERSION           The package version. (required)

-t
--tunnel                    Create a ssh tunnel for upload

--tunnel-local-port PORT    Tunnel local port. Default $TUNNEL_LOCAL_PORT
--tunnel-remote-port PORT   Tunnel remote port. Default $TUNNEL_REMOTE_PORT
--tunnel-remote-host HOST   Tunnel remote host. Default $TUNNEL_REMOTE_HOST
--tunnel-user USER          Tunnel user. Default $TUNNEL_USER
--tunnel-host HOST          Tunnel host. Default $TUNNEL_HOST

--upload-port PORT          Upload port. Default $UPLOAD_PORT
--upload-host HOST          Upload host. Default $UPLOAD_HOST
--upload-path PATH          Upload path. Default $UPLOAD_PATH
--upload-user USER          Upload user. Default $UPLOAD_USER

--help                      Show this help and exit.

EOF
}

ARGS=""
while [ $# -gt 0 ]; do
	case $1 in
		-v | --version ) shift; VERSION="$1" ;;

		-t | --tunnel ) TUNNEL="yes" ;;

		--tunnel-local-port ) shift; TUNNEL_LOCAL_PORT="$1" ;;
		--tunnel-remote-port ) shift; TUNNEL_REMOTE_PORT="$1" ;;
		--tunnel-remote-host ) shift; TUNNEL_REMOTE_HOST="$1" ;;
		--tunnel-user ) shift; TUNNEL_USER="$1" ;;
		--tunnel-host ) shift; TUNNEL_HOST="$1" ;;

		--upload-port ) shift; UPLOAD_PORT="$1" ;;
		--upload-host ) shift; UPLOAD_HOST="$1" ;;
		--upload-path ) shift; UPLOAD_PATH="$1" ;;
		--upload-user ) shift; UPLOAD_USER="$1" ;;

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

# Package =========================================================================================

PACKAGE_NAME="intogen-sm-$VERSION"
PACKAGE_PATH="$DISTRIB_DIR/$PACKAGE_NAME.tar.gz"

# Upload to server =====================================================================================

if [ "$TUNNEL" == "yes" ]; then
    log info "Creating a tunnel through $TUNNEL_USER@$TUNNEL_HOST to $TUNNEL_REMOTE_HOST:$TUNNEL_REMOTE_PORT ..."
    ssh -N -f -L $TUNNEL_LOCAL_PORT:$TUNNEL_REMOTE_HOST:$TUNNEL_REMOTE_PORT $TUNNEL_USER@$TUNNEL_HOST
fi

log info "Uploading package into $UPLOAD_USER@$UPLOAD_HOST:$UPLOAD_PATH ..."

scp -P$UPLOAD_PORT $PACKAGE_PATH $UPLOAD_USER@$UPLOAD_HOST:$UPLOAD_PATH
