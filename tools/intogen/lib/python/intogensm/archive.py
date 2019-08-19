import os.path
import re
import zipfile
import tarfile
import gzip
import bz2

_TAR_EXT = re.compile(r".*\.tar(?:\.(gz|bz2))?$", re.I)

class AbstractArchive(object):
	def __init__(self, path, mode):
		self.path = path
		self.mode = mode

	def close(self):
		if self._arc is not None:
			self._arc.close()
			self._arc = None

class ArchiveEntry(object):
	def __init__(self, arc, info, name, size):
		self._arc = arc
		self._info = info
		self.name = name
		self.size = size

# Zip ------------------------------------------------------------------------------------------------------------------

class ZipEntry(ArchiveEntry):
	def __init__(self, arc, info):
		ArchiveEntry.__init__(self, arc, info, info.filename, info.file_size)

	def extract(self, path=None):
		self._arc.extract(self._info, path)

class ZipArchive(AbstractArchive):
	def __init__(self, path, mode):
		AbstractArchive.__init__(self, path, mode)

		self._arc = zipfile.ZipFile(self.path, mode)
	
	def add(self, filename, arcname):
		self._arc.write(filename, arcname)

	def list(self):
		return [ZipEntry(self._arc, info) for info in self._arc.infolist()]

# Tar ------------------------------------------------------------------------------------------------------------------

class TarEntry(ArchiveEntry):
	def __init__(self, arc, info):
		ArchiveEntry.__init__(self, arc, info, info.name, info.size)

	def extract(self, path=None):
		self._arc.extractall(path, [self._info])

class TarArchive(AbstractArchive):
	def __init__(self, path, mode):
		AbstractArchive.__init__(self, path, mode)

		m = _TAR_EXT.match(path)
		if m:
			self.compression = ":{0}".format(m.group(1))
		else:
			self.compression = ""

		mode += self.compression
		self._arc = tarfile.open(self.path, mode)
	
	def add(self, filename, arcname):
		self._arc.add(filename, arcname)

	def list(self):
		return [TarEntry(self._arc, info) for info in self._arc.getmembers() if info.isfile()]

# Gzip -----------------------------------------------------------------------------------------------------------------

class GzipEntry(ArchiveEntry):
	def __init__(self, arc, name, size):
		ArchiveEntry.__init__(self, arc, None, name, size)

	def extract(self, path=None):
		if path is None:
			path = os.getcwd()
		path = os.path.join(path, self.name)
		f = open(path, "wb")
		buf_size = 1024 * 1024
		buf = self._arc.read(buf_size)
		while buf:
			f.write(buf)
			buf = self._arc.read(buf_size)
		f.close()

class GzipArchive(AbstractArchive):
	def __init__(self, path, mode):
		AbstractArchive.__init__(self, path, mode)

		self._arc = gzip.GzipFile(path, mode=mode)

		size = os.path.getsize(path)

		name = os.path.basename(path)
		if name.lower().endswith(".gz"):
			name = name[:-3]

		self._entries = [GzipEntry(self._arc, name, size)]

	def add(self, filename, arcname):
		raise Exception("Unsupported operation")

	def list(self):
		return self._entries

# Bzip2 ----------------------------------------------------------------------------------------------------------------

class Bzip2Entry(ArchiveEntry):
	def __init__(self, arc, name, size):
		ArchiveEntry.__init__(self, arc, None, name, size)

	def extract(self, path=None):
		if path is None:
			path = os.getcwd()
		path = os.path.join(path, self.name)
		f = open(path, "wb")
		buf_size = 1024 * 1024
		buf = self._arc.read(buf_size)
		while buf:
			f.write(buf)
			buf = self._arc.read(buf_size)
		f.close()

class Bzip2Archive(AbstractArchive):
	def __init__(self, path, mode):
		AbstractArchive.__init__(self, path, mode)

		self._arc = bz2.BZ2File(path, mode=mode)

		size = os.path.getsize(path)

		name = os.path.basename(path)
		if name.lower().endswith(".bz2"):
			name = name[:-4]

		self._entries = [Bzip2Entry(self._arc, name, size)]

	def add(self, filename, arcname):
		raise Exception("Unsupported operation")

	def list(self):
		return self._entries

__FORMAT = {
	"zip" : ZipArchive,
	"tar" : TarArchive,
	"gz" : GzipArchive,
	"bz2" : Bzip2Archive
}

def Archive(path, mode="r", fmt=None):
	if fmt is None:
		name, ext = os.path.splitext(path)
		if len(ext) > 0 and ext[0] == ".":
			ext = ext[1:]

		if _TAR_EXT.match(path):
			fmt = "tar"
		elif ext in __FORMAT:
			fmt = ext
		else:
			raise Exception("Unrecognised compression format for file {0}".format(path))
	
	if fmt not in __FORMAT:
		raise Exception("Unsupported format: {0}".format(fmt))

	return __FORMAT[fmt](path, mode)

