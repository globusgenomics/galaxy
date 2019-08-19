#!/bin/sh
#
# Galaxy init script written by Michael Rusch and Peter Li
#

### BEGIN INIT INFO
# Provides:                     galaxy
# Required-Start:               $local_fs $remote_fs $network
# X-UnitedLinux-Should-Start:
# Required-Stop:                $local_fs $remote_fs $network
# X-UnitedLinux-Should-Stop:
# Default-Start:                3 5
# Default-Stop:                 0 1 2 6
# Short-Description:            Galaxy daemon
# Description:                  Start the Galaxy daemon
### END INIT INFO

#For log information, see ./paster.log

GALAXY_RUN="./run.sh"
GALAXY_USER="huayan"

case "$1" in
    start)
        echo "Starting galaxy..."
        sudo -u $GALAXY_USER $GALAXY_RUN --daemon
        rc_done=$rc_done_up
        #Check status
        sudo -u $GALAXY_USER $GALAXY_RUN --status
        ;;
    stop)
        echo "Stopping galaxy..."
        sudo -u $GALAXY_USER $GALAXY_RUN --stop-daemon
        ;;
    restart)
        $0 stop
        $0 start
        ;;
    *)
        echo "Usage: $0 start|stop|restart"
        exit 1
esac
