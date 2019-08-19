import os
import logging
import importlib

def create_email_sender(name, conf):
	"""
		Creates an email sender dynamically from a module name.
		conf is a dictionary that will be passed to the class constructor as **conf.
	"""
	
	log = logging.getLogger("intogensm.smtp")
	
	name = name.lower()
	mfile = os.path.join(os.path.dirname(__file__), "{0}.py".format(name))
	if not os.path.exists(mfile):
		log.error("SMTP module not found: {0}".format(name))
		return None
	
	try:
		m = importlib.import_module("intogensm.smtp.{0}".format(name))

		cls_name = name[0].upper() + name[1:]
		if not hasattr(m, cls_name):
			log.error("SMTP class not found: {0}".format(cls_name))
			return None

		cls = getattr(m, cls_name)
		instance = cls(**conf)
	except Exception as ex:
		log.error("Error creating email sender: {0}".format(name))
		log.error(" ".join(ex.args))
		return None
	
	return instance
