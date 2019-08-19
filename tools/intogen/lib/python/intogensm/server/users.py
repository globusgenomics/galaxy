import re
import hashlib
import uuid

from datetime import datetime

from flask.ext.login import UserMixin, AnonymousUser

# Constants

_DT_FORMAT = "%Y-%m-%d %H:%M:%S"

_TOKEN_RE = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}")

_EMAIL_RE = re.compile("([^@]+)@[^@]+\.[^@]+")

# Roles

ADMIN_ROLE = "admin"
USER_ROLE = "user"

# Groups

ADMIN_GROUP = "admin"
USERS_GROUP = "users"

# Defaults

_MAX_ANALYSIS = 2
_SIZE_LIMIT = 1000

# Helper functions

def limit_max(a, b):
	if a == -1 or b == -1:
		return -1
	elif a is None:
		return b
	elif b is None:
		return a
	else:
		return max(a, b)

# Exceptions

class SignupException(Exception):
	pass

# Model
								
class User(UserMixin):
	def __init__(self, id, name, roles=None, groups=None, created_time=None, validated=False, max_analysis=None, size_limit=None, active=True):
		self.id = id
		self.name = name
		self.roles = roles or set()
		self.groups = groups or set()
		self.created_time = created_time
		self.validated = validated
		self.max_analysis = max_analysis
		self.size_limit=size_limit
		self.active = active

	def is_active(self):
		return self.active
		
	def __repr__(self):
		return "{1} <{0}>".format(self.id, self.name)

class Anonymous(AnonymousUser):
	name = u"Anonymous"
	
	def __init__(self):
		self.roles = set(["anonymous"])
		self.groups = set(["anonymous"])
		self.created_time = datetime.now()
		self.validated = True
		self.max_analysis = _MAX_ANALYSIS
		self.size_limit = _SIZE_LIMIT

# Manager

class UsersManager(object):
	def __init__(self, email_sender_factory=None):
		self.email_sender_factory = email_sender_factory
		
	def init_db(self, db, conn):
		c = conn.cursor()
		c.execute("SELECT count(*) FROM sqlite_master WHERE type='table' AND name='users'") # FIXME there are more tables
		res = c.fetchone()
		c.close()
		count, = res
		if count == 0:
			self.__create_tables(conn)

	def __create_tables(self, conn):
		c = conn.cursor()
		c.execute("""
		CREATE TABLE roles (
			role_id          TEXT PRIMARY KEY,
			name             TEXT,
			desc             TEXT
		)""")
		
		c.execute("INSERT INTO roles VALUES ('{0}', 'Administrator', 'Can perform administrative tasks')".format(ADMIN_ROLE))
		c.execute("INSERT INTO roles VALUES ('{0}', 'User', 'Can perform analysis')".format(USER_ROLE))
		
		c.execute("""
		CREATE TABLE groups (
			group_id         TEXT PRIMARY KEY,
			name             TEXT,
			
			max_analysis     INTEGER,
			size_limit       INTEGER
		)""")
		
		c.execute("INSERT INTO groups (group_id, name) VALUES ('{0}', 'Administrators')".format(ADMIN_GROUP))
		c.execute("INSERT INTO groups (group_id, name) VALUES ('{0}', 'Users')".format(USERS_GROUP))
		
		c.execute("""
		CREATE TABLE users (
			user_id          TEXT PRIMARY KEY,
			pass_hash        TEXT,
			name             TEXT,
			created_time     TEXT,
			validated        INTEGER DEFAULT 0,
			validation_token TEXT UNIQUE,
			
			max_analysis     INTEGER,
			size_limit       INTEGER
		)""")
		
		c.execute("""
		CREATE TABLE user_roles (
			user_id          TEXT,
			role_id          TEXT,
			
			PRIMARY KEY (user_id, role_id)
		)""")
		
		c.execute("""CREATE TABLE group_roles (
			group_id         TEXT,
			role_id          TEXT,

			PRIMARY KEY (group_id, role_id)
		)""")

		c.execute("INSERT INTO group_roles (group_id, role_id) VALUES ('{0}', '{1}')".format(ADMIN_GROUP, ADMIN_ROLE))
		c.execute("INSERT INTO group_roles (group_id, role_id) VALUES ('{0}', '{1}')".format(USERS_GROUP, USER_ROLE))
		
		c.execute("""
		CREATE TABLE user_groups (
			user_id          TEXT,
			group_id         TEXT,
			
			PRIMARY KEY (user_id, group_id)
		)""")
		
		c.close()

	def __pass_hash(self, password):
		return hashlib.md5(password).hexdigest()
	
	def __send_validation(self, user, validation_url, validation_token):
		if self.email_sender_factory is None or validation_url is None:
			return

		sender = self.email_sender_factory()
		if sender is None:
			return
		
		url = validation_url + "?token={0}".format(validation_token)
		
		subject = "[intogensm] Welcome to IntOGen SM, please validate your email address"
		msg = """
		Hello {0},
		
		welcome to IntOGen SM pipeline demo. Please follow the next link to validate your account:
		
		{1}
		
		If you didn't signed up forget this email and the account will be deleted in some days.
		
		Sincerely yours,
		
		The IntOGen SM team
		""".format(user.name, url)
		
		sender.send(to=user.id, subject=subject, msg=msg)
		
		
	def authenticate(self, conn, user_id, password):
		c = conn.cursor()
		c.execute("SELECT pass_hash FROM users WHERE user_id=?", (user_id,))
		res = c.fetchone()
		c.close()
		
		if res is None:
			return False
		
		pass_hash, = res
		if pass_hash is None:
			return False
		
		return pass_hash == self.__pass_hash(password)
	
	def register(self, conn, user_id, password, name=None, validation_url=None):
	
		m = _EMAIL_RE.match(user_id)
		if not m:
			raise SignupException("Invalid user email")
		
		if password is None or len(password) < 6:
			raise SignupException("Password must contain at least 6 characters")
			
		if name is None or len(name) == 0:
			name = m.group(1)
		
		if self.exists_user(conn, user_id):
			raise SignupException("User {0} already exists".format(user_id))
		
		pass_hash = self.__pass_hash(password)
		
		created_time = datetime.now()
		
		user = User(user_id, name, created_time=created_time)
		
		validation_token = str(uuid.uuid4())
		
		if validation_url is not None:
			self.__send_validation(user, validation_url, validation_token)
			validated = 0
		else:
			validated = 1
		
		c = conn.cursor()
		c.execute("""INSERT INTO users(user_id, pass_hash, name, created_time, validated, validation_token)
						VALUES (?, ?, ?, ?, ?, ?)""", (
						user_id, pass_hash, name, created_time.strftime(_DT_FORMAT),
						validated, validation_token))
		c.execute("INSERT INTO user_groups VALUES (?, ?)", (user_id, USERS_GROUP))
		c.close()
		
		return user
	
	def validate(self, conn, token):
		if not _TOKEN_RE.match(token):
			return None

		c = conn.cursor()
		c.execute("SELECT user_id FROM users WHERE validation_token=?", (token,))
		res = c.fetchone()
		c.close()
		if res is None:
			return None
		user_id, = res
		c = conn.cursor()
		c.execute("UPDATE users SET validated=1 WHERE user_id=?", (user_id,))
		c.close()
		return self.get_user(conn, user_id)
	
	def exists_user(self, conn, user_id):
		c = conn.cursor()
		c.execute("SELECT COUNT(*) FROM users WHERE user_id=?", (user_id,))
		res = c.fetchone()
		c.close()
		if res is None:
			return False
		return res[0] > 0
		
	def get_user(self, conn, user_id):
		c = conn.cursor()
		
		# user
		c.execute("""SELECT user_id, name, created_time, validated, max_analysis, size_limit
						FROM users WHERE user_id=?""", (user_id,))
		res = c.fetchone()
		if res is None:
			c.close()
			return None
		user_id, name, created_time, validated, max_analysis, size_limit = res
		
		# roles
		roles = set()
		
		c.execute("""SELECT role_id FROM user_roles WHERE user_id = ?
						UNION SELECT role_id FROM user_groups JOIN group_roles USING (group_id) WHERE user_id = ?""", (user_id, user_id))
		res = c.fetchone()
		while res:
			role_id, = res
			roles.add(role_id)
			res = c.fetchone()
	
		# groups
		groups = set()
		max_analysis = max_analysis or _MAX_ANALYSIS
		size_limit = size_limit or _SIZE_LIMIT
		
		c.execute("SELECT g.group_id, max_analysis, size_limit FROM groups g JOIN user_groups WHERE user_id = ?", (user_id,))
		res = c.fetchone()
		while res:
			group_id, group_max_analysis, group_size_limit = res
			groups.add(group_id)
			max_analysis = limit_max(max_analysis, group_max_analysis)
			size_limit = limit_max(size_limit, group_size_limit)
			res = c.fetchone()
			
		c.close()
					
		return User(user_id, name, roles=roles, groups=groups,
					created_time=datetime.strptime(created_time, _DT_FORMAT),
					validated=(validated != 0), max_analysis=max_analysis, size_limit=size_limit)

	def remove_user(self, conn, user_id):
		if not self.exists_user(conn,user_id):
			return

		c = conn.cursor()
		c.execute("DELETE FROM user_groups WHERE user_id=?", (user_id,))
		c.execute("DELETE FROM user_roles WHERE user_id=?", (user_id,))
		c.execute("DELETE FROM users WHERE user_id=?", (user_id,))
		c.close()

	def get_user_ids(self, conn):
		c = conn.cursor()
		c.execute("SELECT user_id FROM users ORDER BY created_time")
		user_ids = []
		res = c.fetchone()
		while res is not None:
			user_ids += [res[0]]
			res = c.fetchone()
		c.close()
		return user_ids

	def close(self):
		pass

