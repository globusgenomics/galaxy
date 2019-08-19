from stream import *

class OboParserException(Exception):
	def __init__(self, message, line = -1):
		Exception.__init__(self, message, line)
		self.message = message
		self.line = line

class Tag(object):
	def __init__(self, name, content, line = -1, comment = None):
		self.name = name
		self.content = content
		self.line = line
		self.comment = comment

	def __str__(self):
		return self.content

	def __repr__(self):
		sb = []
		sb.append("%s = %s" % (self.name, self.content))
		if self.comment is not None:
			sb.append(" ! %s" % self.comment)
		return "".join(sb)

class TagCollection(object):
	def __init__(self):
		self.tag_map = {}

	def get_tags(self):
		tags = []
		for e in self.tag_map.values():
			if type(e) == list:
				tags += e
			else:
				tags += [e]
		return tags

	def get_tag(self, name):
		return self.tag_map[name]
		
	def contains_tag(self, name):
		return name in self.tag_map
		
	def is_tag_list(self, name):
		return self.contains_tag(name) and type(self.get_tag(name)) == list
	
	def add_tag(self, tag):
		if tag.name in self.tag_map:
			prev_tag = self.tag_map[tag.name]
			if type(prev_tag) == list:
				prev_tag.append(tag)
			else:
				self.tag_map[tag.name] = [prev_tag, tag]
		else:
			self.tag_map[tag.name] = tag
	
	def add_tag_list(self, tag):
		if tag.name in self.tag_map:
			self.add_tag(tag)
		else:
			self.tag_map[tag.name] = [tag]

	def __repr__(self):
		sb = []
		if self.contains_tag("id"):
			sb.append("%s\n" % repr(self.get_tag("id")))

		for name,tag in self.tag_map.items():
			if name != "id":
				if type(tag) != list:
					sb.append("%s\n" % repr(tag))
				else:
					for tag2 in tag:
						sb.append("%s\n" % repr(tag2))

		return "".join(sb)
		
class Stanza(TagCollection):
	def __init__(self, name, line = -1, comment = None):
		TagCollection.__init__(self)
		self.name = name
		self.line = line
		self.comment = comment

	def has_id(self):
		return self.contains_tag("id")

	def get_id(self):
		return str(self.get_tag("id"))
	
	def __str__(self):
		return self.get_id()
			
	def __repr__(self):
		sb = ["[%s]" % self.name]
		if self.comment is not None:
			sb.append(" ! %s" % self.comment)
		sb.append("\n")
		sb.append(TagCollection.__repr__(self))

		return "".join(sb)

class Ontology(TagCollection):
	def __init__(self):
		TagCollection.__init__(self)
		self.files = []
		self.stanza_map = {}

	def get_stanzas(self, stype):
		return self.stanza_map[stype].values()
	
	def contains_stanza(self, stanza):
		return stanza.name in self.stanza_map and stanza.get_id() in self.stanza_map[stanza.name]

	def get_stanza(self, stype, stanza_id):
		return self.stanza_map[stype][stanza_id]
	
	def add_stanza(self, stanza):
		if stanza.name not in self.stanza_map:
			self.stanza_map[stanza.name] = {}

		if self.contains_stanza(stanza):
			prev_stanza = self.get_stanza(stanza.name, stanza.get_id())
			for tag in stanza.get_tags():
				prev_stanza.add_tag(tag)
		else:
			self.stanza_map[stanza.name][stanza.get_id()] = stanza
		
	def __str__(self):
		return ", ".join(self.files)

	def __repr__(self):
		sb = ["Files: %s\n\n" % ", ".join(self.files)]
		TagCollection.__init__(self)

		sb.append("\n")
		for stanza_type in self.stanza_map.values():
			for stanza in stanza_type.values():
				sb.append("%s\n" % repr(stanza))

		return "".join(sb)

# is_a

class DefTag(Tag):
	def __init__(self, name, content, line = -1, comment = None):
		Tag.__init__(self, name, content, line, comment)
		m = re.match(r'^"(.*)"\s+\[(.*)\]$', content)
		if m is None:
			raise OboParserException("'def' syntax error at line {}".format(line), line)
		self.definition = m.group(1)
		self.dbxref = m.group(2)

	def __repr__(self):
		sb = []
		sb.append('%s = "%s" [%s]' % (self.name, self.definition, self.dbxref))
		if self.comment is not None:
			sb.append(" ! %s" % self.comment)
		return "".join(sb)

tag_info = {
	"header" : {
		"format-version" : {"required" : True, "unique" : True},
		"" : {"required" : False, "unique" : True},
	},
	"term" : {
		"id" : {"required" : True, "unique" : True},
		"name" : {"required" : True, "unique" : True},
		"is_anonymous" : {"required" : False, "unique" : True},
		"alt_id" : {"required" : False, "unique" : False},
		"def" : {"required" : False, "unique" : True, "class" : DefTag},
		"comment" : {"required" : False, "unique" : True},
		"subset" : {"required" : False, "unique" : False},
		"synonym" : {"required" : False, "unique" : False},
		"xref" : {"required" : False, "unique" : False},
		"is_a" : {"required" : False, "unique" : False},
		"intersection_of" : {"required" : False, "unique" : False},
		"union_of" : {"required" : False, "unique" : False},
		"disjoint_from" : {"required" : False, "unique" : False},
		"relationship" : {"required" : False, "unique" : False},
		"is_obsolete" : {"required" : False, "unique" : True},
		"replaced_by" : {"required" : False, "unique" : False},
		"consider" : {"required" : False, "unique" : False}
	},
	"typedef" : {
		"id" : {"required" : True, "unique" : True},
		"name" : {"required" : True, "unique" : True},
		"is_anonymous" : {"required" : False, "unique" : True},
		"alt_id" : {"required" : False, "unique" : False},
		"def" : {"required" : False, "unique" : True},
		"comment" : {"required" : False, "unique" : True},
		"subset" : {"required" : False, "unique" : False},
		"synonym" : {"required" : False, "unique" : False},
		"xref" : {"required" : False, "unique" : False},
		"is_a" : {"required" : False, "unique" : False},
		"relationship" : {"required" : False, "unique" : False},
		"is_obsolete" : {"required" : False, "unique" : True},
		"replaced_by" : {"required" : False, "unique" : False},
		"consider" : {"required" : False, "unique" : False},
		"domain" : {"required" : False, "unique" : True},
		"range" : {"required" : False, "unique" : True},
		"inverse_of" : {"required" : False, "unique" : False},
		"transitive_over" : {"required" : False, "unique" : False},
		"is_cyclic" : {"required" : False, "unique" : True},
		"is_reflexive" : {"required" : False, "unique" : True},
		"is_symmetric" : {"required" : False, "unique" : True},
		"is_anti_symmetric" : {"required" : False, "unique" : True},
		"is_transitive" : {"required" : False, "unique" : True},
		"is_metadata_tag" : {"required" : False, "unique" : True}
	},
	"instance" : {
		"id" : {"required" : True, "unique" : True},
		"name" : {"required" : True, "unique" : True},
		"instance_of" : {"required" : True, "unique" : True},
		"property_value" : {"required" : False, "unique" : False},
		"is_anonymous" : {"required" : False, "unique" : True},
		"alt_id" : {"required" : False, "unique" : False},
		"comment" : {"required" : False, "unique" : True},
		"synonym" : {"required" : False, "unique" : False},
		"xref" : {"required" : False, "unique" : False},
		"is_obsolete" : {"required" : False, "unique" : True},
		"replaced_by" : {"required" : False, "unique" : False},
		"consider" : {"required" : False, "unique" : False}
	}
}

class OboSimpleParser(object):
	"""
	VERY SIMPLE OBO parser
	"""

	def __init__(self):
		self.ontology = Ontology()

	def _parse_tag(self, evt, tag_collection, stanza_name, stanza_id = None):
		valid_tags = tag_info[stanza_name]

		if evt.tag_name not in valid_tags:
			tag_collection.add_tag(Tag(evt.tag_name, evt.content, evt.line, evt.comment))
		else:
			tag_desc = valid_tags[evt.tag_name]
			tag_class = tag_desc.get("class", Tag)
			unique = tag_desc.get("unique", False)

			tag = tag_class(evt.tag_name, evt.content, evt.line_pos, evt.comment)
			if tag_collection.contains_tag(evt.tag_name):
				prev_tag = tag_collection.get_tag(evt.tag_name)
				if unique:
					if stanza_id is not None:
						msg = "Tag '%s' specified more than one time for %s %s" % (evt.tag_name, stanza_name, stanza_id)
					else:
						msg = "Tag '%s' specified more than one time for a %s" % (evt.tag_name, stanza_name)
					raise OboParserException(msg, evt.line_pos)

			if not unique:
				tag_collection.add_tag_list(tag)
			else:
				tag_collection.add_tag(tag)

	def parse_header_tag(self, evt):
		self._parse_tag(evt, self.ontology, "header")

	def parse_stanza_tag(self, evt, stanza):
		if stanza.has_id():
			stanza_id = stanza.get_id()
		else:
			stanza_id = None

		self._parse_tag(evt, stanza, stanza.name, stanza_id)
		
	def parse(self, filename):
		reader = OboStreamReader(filename)
		self.ontology.files.append(filename)

		evt = reader.next_event()

		# Parse header
		
		while evt is not None and evt.type != OboEventTypes.STANZA_START:
			if evt.type == OboEventTypes.TAG_START:
				self.parse_header_tag(evt)

			evt = reader.next_event()

		# Parse stanzas

		while evt is not None:
			stype = evt.name.lower()
			stanza = Stanza(stype, evt.line)
			evt = reader.next_event()
			while evt is not None and evt.type != OboEventTypes.STANZA_START:
				if evt.type == OboEventTypes.TAG_START:
					self.parse_stanza_tag(evt, stanza)

				evt = reader.next_event()

			if not stanza.contains_tag("id"):
				raise OboParserException("Stanza without id", stanza.line)

			self.ontology.add_stanza(stanza)
		
		reader.close()
		
		return self.ontology

	#TODO
	def validate(self):
		return false

