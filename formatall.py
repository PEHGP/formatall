import re
class format:
	def readfasta(self,f):
		self.d={}
		self.name=f
		for x in open(f):
			x=x.rstrip()
			m=re.search(">(.*)",x)
			if m:
				p=m.group(1)
				self.d[p]=""
				continue
			self.d[p]+=x
	def formataa(self):
		self.r={}
		for x in self.d:
			da={}
			for i in self.d[x]:
				if not i in da.keys():
					da[i]=1
				else:
					da[i]+=1
			for i in da:
				da[i]=da[i]/float(len(self.d[x]))
			self.r[x]=da
	def hodc(self,seq,n):
		l=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
		self.dc={}
		for x in l:
			for y in l:
				self.dc[x+"x"*(n-1)+y]=0.
		i=0
		while i<=len(seq)-n-1:
			s=seq[i]+"x"*(n-1)+seq[i+n]
			if s in self.dc.keys():
				self.dc[s]+=1
			i+=1
		for x in self.dc.keys():
			self.dc[x]=self.dc[x]/(len(seq)-n)
	def formatcksaap(self,n):
		self.r={}
		for x in self.d:
			dll={}
			for i in range(1,n):
				self.hodc(self.d[x],i)
				dll.update(self.dc)
			self.r[x]=dll
	def formatcsv(self,f):
		aa=self.r[self.r.keys()[0]].keys()
		t="protein"+"\t"
		for x in aa:
			t+=x+"\t"
		f.write(t[:-1]+"\n")
		for x in self.r:
			fm=x+"\t"
			for i in self.r[x]:
				fm+=str(self.r[x][i])+"\t"
			f.write(fm[:-1]+"\n")
	def formatarff(self,f,cl):
		f.write("@relation "+self.name+"\n\n@attribute protein string\n")
		aa=self.r[self.r.keys()[0]].keys()
		for x in aa:
			f.write("@attribute "+x+" numeric\n")
		f.write("@attribute class {"+cl+"}\n\n@data\n\n")
		for x in self.r:
			fm=x+","
			for i in self.r[x]:
				fm+=str(self.r[x][i])+","
			f.write(fm[:-1]+"\n")
	def formatlibsvm(self,f,cl):
		if cl=="yes":
			flag="1"
		elif cl=="no":
			flag="-1"
		else:
			print "class is not right"
			exit(0)
		for x in self.r:
			fm=flag+" "
			t=1
			for i in self.r[x]:
				fm+=str(t)+":"+str(self.r[x][i])+" "
				t+=1
			f.write(fm[:-1]+"\n")
if __name__=="__main__":
	p=format()
	p.readfasta("eff104.fasta")
	p.formataa()
	f=open("aa.csv","w")
	p.formatcsv(f)
	f.close()
	#f=open("t.arff","w")
	#p.formatarff(f,"yes,no")
	#f.close()
	#f=open("t.libsvm","w")
	#p.formatlibsvm(f,"yes")
	
