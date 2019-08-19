
function log {
	local COLOR=""
	case $1 in
		debug ) local COLOR="\033[1;34m";;
		info ) local COLOR="\033[1;32m";;
		warn ) local COLOR="\033[1;33m";;
		error ) local COLOR="\033[1;31m";;
	esac
	local TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
	echo -e "$COLOR$TIMESTAMP $2\033[0m"
}

function cmd_exists {
	if command -v $1 >/dev/null 2>&1; then
		return 0
	else
		return 1
	fi
}

function pylib_exists {
	local freeze=$(pip freeze | grep -E "^$1$")
	if [ -n "$freeze" ]; then
		return 0
	else
		return 1
	fi
}

function pylib_install {
	LIB=$1
	if ! pylib_exists $LIB; then
		log info "Installing required Python library '$LIB' ..."
		if ! pip install --upgrade $LIB; then
			log error "ERROR: Installation of Python library '$LIB' failed."
			log error "Please check that your internet connection is working."
			exit -1
		fi
	else
		log info "Required Python library '$LIB' already installed"
	fi
}

function pylib_install_url {
	LIB=$1
	URL=$2
	if ! pylib_exists $LIB; then
		log info "Installing required Python library '$LIB' ..."
		pip uninstall --yes $LIB 2>/dev/null
		if ! pip install --upgrade $URL#egg=$LIB; then
			log error "ERROR: Installation of Python library '$LIB' failed."
			log error "Please check that your internet connection is working."
			exit -1
		fi
	else
    		log info "Required Python library '$LIB' already installed"
	fi
}

function pylib_install_dev {
	LIB=$1
	URL=$2
	REF=$3
	[ -z "$REF" ] && REF="develop"
	log info "Installing required Python library '$LIB' ..."
	pip uninstall --yes $LIB 2>/dev/null
	if ! pip install --upgrade -e git+$URL@${REF}#egg=$LIB==dev; then
		log error "ERROR: Installation of Python library '$LIB' failed."
		log error "Please check that your internet connection is working."
		exit -1
	fi
}

function perlib_exists {
	if $PERL_BIN -M$1 -e "" >/dev/null 2>&1; then
		return 0
	else
		return 1
	fi
}

function add_python_path {
	if [ -n "$PYTHONPATH" ]; then
		PYTHONPATH=$PYTHONPATH:$1
	else
		PYTHONPATH=$1
	fi
}

function add_perl_path {
	if [ -n "$PERL5LIB" ]; then
		PERL5LIB=$PERL5LIB:$1
	else
		PERL5LIB=$1
	fi
}

function lowercase {
# echo "$1" | sed "y/ABCDEFGHIJKLMNOPQRSTUVWXYZ/abcdefghijklmnopqrstuvwxyz/"
# echo "$1" | perl -e 'print lc <>;'
python <<EOF
print "$1".lower()
EOF
}

function browse {
	if [ -n "$BROWSER" ]; then
		eval $BROWSER "$1"
	elif cmd_exists xdg-open; then
		xdg-open "$1"
	elif cmd_exists gnome-open; then
		gnome-open "$1"
	elif cmd_exists firefox; then
		firefox "$1"
	elif cmd_exists google-chrome; then
		google-chrome "$1"
	else
		log warn "WARNING: Could not find any web browser to open the url $1"
	fi
}

function pid_alive {
	if [ -n "$(ps -eo pid | grep $1)" ]; then
		return 0
	else
		return 1
	fi
}

function stop {
	local PID_FILE=$1
	local SIGNAL=$2
	[ "$SIGNAL" == "" ] && SIGNAL="-INT"
	
	if [ -f $PID_FILE ]; then
		PID=$(cat $PID_FILE)
		log info "Stopping the process with pid $PID ..."
		kill $SIGNAL $PID >/dev/null 2>&1
	
		count=60
		while true; do
			if ! pid_alive $PID; then
				break
			fi
		
			sleep 1
		
			let count=$count-1
			if [ $count -eq 0 ]; then
				break
			fi
		done
	
		if pid_alive $PID; then
			log warn "WARNING: The process with pid $PID couldn't be stopped. Killing it !"
			kill -9 $PID >/dev/null 2>&1
			exit -1
		else
			rm -f $PID_FILE
			exit 0
		fi
	else
		log warn "WARNING: No process was found to be running."
	fi
}

function path2json {
python <<EOF
path = "$1".split(":")
import json
print json.dumps(path)
EOF
}


