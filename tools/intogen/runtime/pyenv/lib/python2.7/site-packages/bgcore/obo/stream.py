import re

from collections import deque

class OboEventTypes(object):
	UNKNOWN = -1
	COMMENT = 1
	DOCUMENT_START = 10
	DOCUMENT_END = 19
	HEADER_START = 20
	HEADER_END = 29
	STANZA_START = 30
	STANZA_END = 39
	TAG_START = 40
	TAG_END = 49

class OboEvent(object):
	def __init__(self, type, line_pos = -1, line = "", comment = None):
		self.type = type
		self.line_pos = line_pos
		self.line = line
		self.comment = comment

	def __str__(self):
		name = {
			OboEventTypes.UNKNOWN: "UNKNOWN",
			OboEventTypes.COMMENT: "COMMENT",
			OboEventTypes.DOCUMENT_START: "DOCUMENT_START",
			OboEventTypes.DOCUMENT_END: "DOCUMENT_END",
			OboEventTypes.HEADER_START: "HEADER_START",
			OboEventTypes.HEADER_END: "HEADER_END",
			OboEventTypes.STANZA_START: "STANZA_START",
			OboEventTypes.STANZA_END: "STANZA_END",
			OboEventTypes.TAG_START: "TAG_START",
			OboEventTypes.TAG_END: "TAG_END",
			}.get(self.type, "UNKNOWN TYPE")

		return "%04i : %s" % (self.line_pos, name)

class OboDocumentStartEvent(OboEvent):
	def __init__(self):
		OboEvent.__init__(self, OboEventTypes.DOCUMENT_START)

class OboDocumentEndEvent(OboEvent):
	def __init__(self, line_pos):
		OboEvent.__init__(self, OboEventTypes.DOCUMENT_END, line_pos=line_pos)

class OboStanzaStartEvent(OboEvent):
	def __init__(self, name, **kv):
		OboEvent.__init__(self, OboEventTypes.STANZA_START, **kv)
		self.name = name

	def __str__(self):
		return OboEvent.__str__(self) + " " + self.name

class OboStanzaEndEvent(OboEvent):
	def __init__(self, name, **kv):
		OboEvent.__init__(self, OboEventTypes.STANZA_END, **kv)
		self.name = name

	def __str__(self):
		return OboEvent.__str__(self) + " " + self.name

class OboTagStartEvent(OboEvent):
	def __init__(self, stanza_name, tag_name, content, **kv):
		OboEvent.__init__(self, OboEventTypes.TAG_START, **kv)
		self.stanza_name = stanza_name
		self.tag_name = tag_name
		self.content = content

	def __str__(self):
		return OboEvent.__str__(self) + " " + self.tag_name + " = " + self.content

class OboTagEndEvent(OboEvent):
	def __init__(self, stanza_name, tag_name, content, **kv):
		OboEvent.__init__(self, OboEventTypes.TAG_END, **kv)
		self.stanza_name = stanza_name
		self.tag_name = tag_name
		self.content = content

	def __str__(self):
		return OboEvent.__str__(self) + " " + self.tag_name #+ " = " + self.content

class OboStream(object):

	def __init__(self, filename):
		self.filename = filename
		self.file = open(filename, "r")
		self.line_pos = 0

	# returns a non empty line or null if EOF
	def next_line(self):
		line = self.read_line()
		while line is not None and self.is_empty_line(line):
			line = self.read_line()

		return line;

	def internal_read_line(self):
		line = self.file.readline()
		if not line:
			line = None
		else:
			line = line.rstrip()
		return line

	# returns the next line or null if EOF
	def read_line(self):
		completeLine = []

		line = self.internal_read_line()

		self.line_pos += 1

		if line is not None:
			line = line.strip()

		while line is not None and line.endswith("\\"):
			line = line[0:len(line) - 1]
			completeLine.append(line);
			line = self.internal_read_line()
			linePos += 1
			if line is not None:
				line = line.strip()

		if line is not None:
			completeLine.append(line)

		s = ''.join(completeLine)
		if line is None and len(s) == 0:
			return None

		return s

	def is_empty_line(self, line):
		return len(re.sub("\s", "", line)) == 0
		
	def close(self):
		self.file.close()

class OboStreamReader(object):

	STANZA_NAME_PATTERN = "^\[(.*)\][ \t]*(?:!(.*))?$"
	LINE_COMMENT_PATTERN = "^\s*!(.*)$"
	TAG_NAME_PATTERN = "^[0-9a-zA-Z_]$"

	ESCAPE_MAP = {
		'n': '\n',
		'W': ' ',
		't': '\t',
		':': ':',
		',': ',',
		'"': '"',
		'\\': '\\',
		'(': '(',
		')': ')',
		'[': '[',
		']': ']',
		'{': '{',
		'}': '}' }

	def __init__(self, filename):
		self.stream = OboStream(filename)
		self.tokens = deque([OboDocumentStartEvent()])

		self.header_started = False
		self.header_ended = False
		self.document_ended = False
		self.stanza_name = None

	def next_event(self):
		if len(self.tokens) > 0:
			return self.tokens.pop()

		while len(self.tokens) == 0:
			line = self.stream.next_line()
			pos = self.stream.line_pos;

			if line is None:
				if not self.document_ended:
					self.document_ended = True
					self.tokens.appendleft(OboDocumentEndEvent(line_pos = pos))
					break
				else:
					return None

			stz_match = re.match(OboStreamReader.STANZA_NAME_PATTERN, line)
			comment_match = re.match(OboStreamReader.LINE_COMMENT_PATTERN, line)

			if stz_match is not None:
				if self.stanza_name is None:
					if not self.header_ended:
						self.tokens.appendleft(OboEvent(OboEventTypes.HEADER_END, line_pos=pos))
						self.header_ended = True

				stz_name = stz_match.group(1)
				self.tokens.appendleft(OboStanzaStartEvent(stz_name, line_pos=pos))
				self.stanza_name = stz_name

			elif comment_match is not None:
				self.tokens.appendleft(OboEvent(OboEventTypes.COMMENT, line_pos=pos))
			else:
				if self.stanza_name is None and not self.header_started:
					self.tokens.appendleft(OboEvent(OboEventTypes.HEADER_START, line_pos=pos))
					self.header_started = True

				self.next_tag(line, pos)

		return self.tokens.pop()

	def next_tag(self, line, linepos):
		pos = line.find(':')
		if pos < 0:
			self.tokens.appendleft(OboEvent(OboEventTypes.UNKNOWN, line_pos=linepos))
			return

		self.tag_name = line[0:pos]
		content = line[pos + 1:]
		sb = []

		self.escape_chars_and_remove_comments(content, sb)
		content = ''.join(sb).strip();
		self.tokens.appendleft(OboTagStartEvent(self.stanza_name, self.tag_name, content, line_pos=linepos));
		self.tokens.appendleft(OboTagEndEvent(self.stanza_name, self.tag_name, content, line_pos=linepos));

	# replace escape characters and remove comments
	def escape_chars_and_remove_comments(self, line, sb):
		ln = len(line);
		pos = 0;
		while pos < ln:
			c = line[pos]
			pos += 1
			if c == '!':
				pos = ln
			elif c != '\\' or (c == '\\' and pos == ln - 1):
				sb.append(c);
			else:
				c = line[pos]
				pos += 1

				sb.append(self.ESCAPE_MAP.get(c, c))
				
	def close(self):
		self.stream.close()

