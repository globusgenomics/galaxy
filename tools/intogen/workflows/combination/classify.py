import os

from wok.task import task

from intogensm.projdb import ProjectDb
from intogensm.paths import create_combination_folders

CONF_KEY="combination.classifiers"

@task.begin()
def begin():
	log = task.logger
	conf = task.conf

	log.info("Creating combination folders ...")

	create_combination_folders(conf)

	log.info("Checking classifiers ...")

	if CONF_KEY in conf:
		classifiers = conf[CONF_KEY].to_native()
	else:
		classifiers = []

	if not isinstance(classifiers, list):
		log.error("'{0}' should be a list".format(CONF_KEY))
		return -1

	results = []

	if len(classifiers) == 0:
		log.warn("No classifiers have been defined !!!")
		log.warn("Specify them in the configuration with '{0}'".format(CONF_KEY))

	for index, classifier in enumerate(classifiers):
		if not isinstance(classifier, dict):
			raise Exception("Classifier at index {0} should be a dictionary".format(index))

		if "id" not in classifier:
			classifier["id"] = str(index)

		classifier_id = classifier["id"]

		if "name" not in classifier:
			classifier["name"] = "Classifier {0}".format(index)

		if "keys" not in classifier:
			classifier["keys"] = []
		elif not isinstance(classifier["keys"], list):
			raise Exception("'keys' for classifier at index {0} should be a list".format(classifier_id))

		keys_len = len(classifier["keys"])

		if "default_key_values" not in classifier:
			classifier["default_key_values"] = [""] * keys_len
		elif not isinstance(classifier["default_key_values"], list):
			raise Exception("'default_key_values' for classifier {0} should be a list".format(classifier_id))
		elif len(classifier["default_key_values"]) != keys_len:
			raise Exception("Number of values for 'default_key_values' should be the same of 'keys' in classifier '{0}'".format(classifier_id))
		
		if "short_names" not in classifier:
			classifier["short_names"] = classifier["keys"]
		elif not isinstance(classifier["short_names"], list):
			raise Exception("'short_names' for classifier {0} should be a list".format(classifier_id))
		elif len(classifier["short_names"]) != keys_len:
			raise Exception("Number of keys for 'short_names' should be the same of 'keys' in classifier '{0}'".format(classifier_id))

		if "default_short_values" not in classifier:
			classifier["default_short_values"] = classifier["default_key_values"]
		elif not isinstance(classifier["default_short_values"], list):
			raise Exception("'default_short_values' for classifier {0} should be a list".format(classifier_id))
		elif len(classifier["default_short_values"]) != keys_len:
			raise Exception("Number of values for 'default_short_values' should be the same of 'keys' in classifier '{0}'".format(classifier_id))
		
		if "long_names" not in classifier:
			classifier["long_names"] = classifier["short_names"]
		elif not isinstance(classifier["long_names"], list):
			raise Exception("'long_names' for classifier {0} should be a list".format(classifier_id))
		elif len(classifier["long_names"]) != keys_len:
			raise Exception("Number of keys for 'long_names' should be the same of 'keys' in classifier '{0}'".format(classifier["id"]))

		if "default_long_values" not in classifier:
			classifier["default_long_values"] = classifier["default_key_values"]
		elif not isinstance(classifier["default_long_values"], list):
			raise Exception("'default_long_values' for classifier {0} should be a list".format(classifier_id))
		elif len(classifier["default_long_values"]) != keys_len:
			raise Exception("Number of values for 'default_long_values' should be the same of 'keys' in classifier '{0}'".format(classifier_id))

		results += [{}]

	task.context["classifiers"] = classifiers
	task.context["results"] = results

def str_join(delimiter, iterable):
	return delimiter.join([str(x) for x in iterable])

@task.foreach()
def classify(project):
	log = task.logger
	conf = task.conf

	project_id = project["id"]

	log.info("--- [{0}] --------------------------------------------".format(project_id))

	classifiers = task.context["classifiers"]
	groups = []
	results = task.context["results"]

	annotations = project["annotations"]

	for index, classifier in enumerate(classifiers):
		name = classifier["name"]

		keys = classifier["keys"]
		default_key_values = classifier["default_key_values"]

		short_names = classifier["short_names"]
		default_short_values = classifier["default_short_values"]

		long_names = classifier["long_names"]
		default_long_values = classifier["default_long_values"]

		group_values = []
		short_values = []
		long_values = []
		for i, key in enumerate(keys):
			group_values += [annotations.get(key, default_key_values[i])]
			short_values += [annotations.get(short_names[i], default_short_values[i])]
			long_values += [annotations.get(long_names[i], default_long_values[i])]

		group_values = tuple(group_values)
		short_values = tuple(short_values)
		long_values = tuple(long_values)

		log.info("{0}: key=({1}), short=({2}), long=({3})".format(name,
						str_join("; ", group_values),
						str_join("; ", short_values),
						str_join("; ", long_values)))

		cresults = results[index]

		if group_values in cresults:
			cresults[group_values][2].append(project)
		else:
			cresults[group_values] = (short_values, long_values, [project])

@task.end()
def end():
	log = task.logger
	conf = task.conf

	project_sets_port = task.ports("project_sets")

	classifiers = task.context["classifiers"]
	results = task.context["results"]

	log.info("Classification results ...")

	for index, classifier in enumerate(classifiers):
		log.info("  {0}:".format(classifier["name"]))
		cresults = results[index]
		groups = sorted(cresults.keys())
		for group_values in groups:
			group_short_values, group_long_values, projects = cresults[group_values]

			group_name = str_join("-", group_values)
			group_short_name = str_join("-", group_short_values)
			group_long_name = str_join("; ", group_long_values)

			projects_set = (dict(classifier,
								 group_name=group_name,
								 group_values=group_values,
								 group_short_name=group_short_name,
								 group_short_values=group_short_values,
								 group_long_name=group_long_name,
								 group_long_values=group_long_values),
							projects)

			log.info("    ({0}) -> {1} projects".format(group_name, len(projects)))
			project_sets_port.send(projects_set)

task.run()