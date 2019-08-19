from empirical import MeanEmpiricalTest, MedianEmpiricalTest
from zscore import MeanZscoreTest, MedianZscoreTest

__METHODS = [
	MeanEmpiricalTest, MedianEmpiricalTest,
	MeanZscoreTest, MedianZscoreTest
]

__METHODS_BY_NAME = {}
for m in __METHODS:
	__METHODS_BY_NAME[m.NAME] = m

def method_names():
	return sorted(__METHODS_BY_NAME.keys())

def create_method(method):
	method = method.lower()
	if method not in __METHODS_BY_NAME:
		return None
	return __METHODS_BY_NAME[method]()