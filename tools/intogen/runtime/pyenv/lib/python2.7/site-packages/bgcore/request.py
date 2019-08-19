import time
import urllib
import urllib2

class Request(object):
	def __init__(self, max_retries=1, max_freq=None):
		"""
		:param max_retries: Maximum number of retries per request
		:param max_freq: Maximum number of requests per second
		:return:
		"""
		self.max_retries = max_retries
		self.max_freq = max_freq

		if max_freq is not None:
			self.__wait_period = 1.0 / max_freq
			self.__last_time = None
		else:
			self.__wait_period = None
			self.__last_time = None

	def request(self, url, params=None, data=None, headers=None):
		retries = self.max_retries

		if self.__wait_period is not None and self.__last_time is not None:
			elapsed_time = time.time() - self.__last_time
			if elapsed_time < self.__wait_period:
				time.sleep(self.__wait_period - elapsed_time)

		if params is not None:
			if isinstance(params, dict):
				params = urllib.urlencode(params)
			url += "?" + params
			
		if data is not None and isinstance(data, dict):
			data = urllib.urlencode(data)

		if headers is None:
			headers = {}

		response = None

		while response is None and retries > 0:
			self.__last_time = time.time()

			try:
				req = urllib2.Request(url, data, headers)
				response = urllib2.urlopen(req)
			except urllib2.URLError:
				response = None

			retries -= 1

		return response

	def get(self, url, params=None, headers=None):
		return self.request(url, params, None, headers)

	def post(self, url, data=None, headers=None):
		return self.request(url, None, data, headers)
