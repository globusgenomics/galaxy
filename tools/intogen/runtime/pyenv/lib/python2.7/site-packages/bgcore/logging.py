logging = __import__("logging")

_LOG_LEVELS = {
	"debug" : logging.DEBUG,
	"info" : logging.INFO,
	"warn" : logging.WARN,
	"error" : logging.ERROR,
	"critical" : logging.CRITICAL,
	"notset" : logging.NOTSET }

def add_logging_arguments(parser):
	parser.add_argument("-L", "--log-level", dest="log_level", metavar="LEVEL", default="info",
						choices=["debug", "info", "warn", "error", "critical", "notset"],
						help="Define log level: debug, info, warn, error, critical, notset")

def selected_level(args):
	log_level = args.log_level.lower() if args.log_level is not None else "notset"
	return log_level

def initialize(args):
	args.log_level = selected_level(args)

	logging.basicConfig(
		format = "%(asctime)s %(name)s %(levelname)-5s : %(message)s",
		datefmt = "%Y-%m-%d %H:%M:%S",
		level = _LOG_LEVELS[args.log_level])

def get_logger(name, level=None):
	log = logging.getLogger(name)
	if level is not None:
		level = level.lower()
		log.setLevel(_LOG_LEVELS[level])

	return log