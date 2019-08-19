import os
import json

class ExamplesManager(object):
	def __init__(self, path):
		self.path = path
		self.index_path = os.path.join(path, "index.json")

		self.examples = []
		self.examples_by_id = {}
		
		self._load_examples()
		
	def _check_required(self, d, keys):
		ok = True
		for k in keys:
			ok = ok and (k in d)
		return ok
		
	def _load_examples(self):
		if self.index_path is None:
			return
		
		if not os.path.exists(self.index_path):
			return
		
		f = open(self.index_path)
		examples = json.load(f)
		f.close()
		
		req_fields = ["id", "title", "file"]
		
		for example in examples:
			if self._check_required(example, req_fields):
				example["format"] = example.get("format", "tab")
				example["assembly"] = example.get("assembly", "hg19")
				example["path"] = os.path.join(self.path, example["file"])
				example["mime"] = example.get("mime", "text/plain")
				desc = example.get("desc", "")
				if isinstance(desc, list):
					desc = "\n".join(desc)
				example["desc"] = desc
				self.examples += [example]
				self.examples_by_id[example["id"]] = example

	def reload(self):
		self.examples = []
		self.examples_by_id = {}
		self._load_examples()

	def get_example(self, eid):
		if eid not in self.examples_by_id:
			return None
		
		return self.examples_by_id[eid]
	
	def get_examples(self):
		return self.examples

	def close(self):
		pass

