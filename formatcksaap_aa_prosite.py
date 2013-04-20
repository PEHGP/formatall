#!/usr/bin/python
import re,sys,os
def prosite(f):
	d={}
	for x in open(f):
		x=x.rstrip()
		m=re.search(">(.*)",x)
		if m:
			p=m.group(1)
			d[p]={}
	l=[]
	for x in open("motifs.list"):
		x=x.rstrip()
		l.append(x)
	for p in d.keys():
		for x in l:
			d[p][x]=0
	#os.system("./ps_scan.pl %s >%s"%(f,"temp_"+f))#careful
	lr=[]
	for x in open("temp_"+f):#careful
		x=x.rstrip()
		m=re.search(">(.*?)\s:\s(.*?)\s",x)
		if m:
			p=m.group(1)
			if not p in lr:
				lr.append(p)
			s=m.group(2)
			continue
		d[p][s]+=1
	#print len(lr)
	#os.system("rm temp")
	return d
def aafrequency(seq):
	#ACDEFGHIKLMNPQRSTVWY
	d={}
	for x in "ACDEFGHIKLMNPQRSTVWY":
		d[x]=0.
	for x in seq:
		if x in d.keys():
			d[x]+=1
	for x in d.keys():
		#d[x]=(d[x]-dmin)/float(dmax) #may need change
		d[x]=d[x]/len(seq)
	#print len(seq)
	return d
def hodc(seq,n):
	l=["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"]
	d={}
	for x in l:
		for y in l:
			d[x+"x"*(n-1)+y]=0.
	i=0
	while i<=len(seq)-n-1:
		s=seq[i]+"x"*(n-1)+seq[i+n]
		if s in d.keys():
			d[s]+=1
		i+=1
	for x in d.keys():
		#d[x]=d[x]/len(seq) #may need change
		d[x]=d[x]/(len(seq)-n) #formatcksaap_aa2
	return d
def cksaap_aa(f):
	d={}
	r={}
	for x in open(f):
		x=x.rstrip()
		m=re.search(">(.*)",x)
		if m:
			p=m.group(1)
			d[p]=""
			continue
		d[p]+=x
	for x in d.keys():
		dll={}
		dll.update(aafrequency(d[x]))
		for i in range(1,7):#need change
			dll.update(hodc(d[x],i))
		r[x]=dll
	return r
def formatcsv(f,flag,fr):
	dp=prosite(f)
	dc=cksaap_aa(f)
	t="protein"+"\t"
	p=dp.keys()[0]
	l=dp[p].keys()+dc[p].keys()
	for x in l:
		t+=x+"\t"
	t+="class"
	if flag=="yes":
		#print t
		fr.write(t+"\n")
	l=dp.keys()
	for p in l:
		s=p+"\t"
		lv=[]
		#te="protein"+"\t"
		for x in dp[p]:
			#te+=x+"\t"
			#lv.append(dp[p][x])
			s+=str(dp[p][x]/2342.0)+"\t"
		for x in dc[p]:
			#te+=x+"\t"
			#lv.append(dc[p][x])
			s+=str(dc[p][x])+"\t"
		#te+="class"
		"""
		ma=max(lv)
		mi=min(lv)
		for v in lv:
			v=(v-mi)/float(ma-mi)
			s+=str(v)+"\t"
		"""
		s+=flag
		#print te
		print s
		fr.write(s+"\n")
if __name__=="__main__":
	for i in range(1,11):
		fr=open("pe_"+str(i)+"_cksaap_aa_prosite.csv","w")
		formatcsv("eff104.fasta","yes",fr)
		formatcsv(str(i),"no",fr)
		fr.close()

