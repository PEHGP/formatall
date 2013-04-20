#!/usr/bin/python
import re
class rfasta:
	def __init__(self,f,ze):
		self.f=f
		self.ze=ze
	def rfa(self):
		self.d={}
		for x in open(self.f):
			x=x.rstrip()
			m=re.search("%s"%self.ze,x)
			if m:
				p=m.group(1)
				self.d[p]=""
				continue
			self.d[p]+=x
