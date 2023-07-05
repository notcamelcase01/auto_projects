
from auto import *
def myRun(demoname):
	r = run(e = demoname, c = demoname)
	branchpoints = r("BP")
	for solution in branchpoints:
		bp = load(solution, ISW=-1)
		# Compute forwards
		print "Solution label", bp["LAB"], "forwards"
		fw = run(bp)
		# Compute backwards
		print "Solution label", bp["LAB"], "backwards"
		bw = run(bp,DS='-')
		both = fw + bw
		merged = merge(both)
		r = r + merged
	r=relabel(r)
	save(r, demoname)
	plot(r)
	wait()
myRun("fol")

'''
from auto import *
def myRun(demoname):
	r = run(e = demoname, c = demoname)
	branchpoints = r("BP")
	for solution in branchpoints:
		bp = load(solution)
		# Compute forwards
		print "Solution label", bp["LAB"], "forwards"
		fw = run(bp)
		# Compute backwards
		print "Solution label", bp["LAB"], "backwards"
		bw = run(bp)
		both = fw + bw
		merged = merge(both)
		r = r + merged
	r=relabel(r)
	save(r, demoname)
	plot(r)
	wait()
myRun("ext")
'''
