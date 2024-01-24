from statistics import mean
import itertools


class Tricluster:

    def __init__(self, nTimes, nSamples, nPatients):
        self.nTimes = int(nTimes)
        self.nSamples = int(nSamples)
        self.nPatients = int(nPatients)
        self.times = []
        self.samples = []
        self.patients = []

        self.cluster = {}  # format : {(t,s,g): val ...}

    def addTime(self, time):
        if time not in self.times and len(self.times) < self.nTimes:
            self.times.append(time)

    def addSample(self, sample):
        if sample not in self.samples and len(self.samples) < self.nSamples:
            self.samples.append(sample)

    def addPatient(self, patient):
        if len(self.patients) < self.nPatients:
            self.patients.append(patient)

    def addValue(self, time, sample, patient, value):
        if time not in self.times:
            self.times.append(time)

        if sample not in self.samples:
            self.samples.append(sample)

        if patient not in self.patients:
            self.patients.append(patient)

        if (time, sample, patient) not in self.cluster:
            self.cluster[(time, sample, patient)] = value

    def hasPatient(self, patient):
        return patient in self.patients

    def hasTimePoint(self, time):
        return time in self.times

    def getTricluster(self):
        return self.cluster

    def getTriclusterCoord(self):
        return self.cluster.keys()

    def getPatients(self):
        return self.patients

    def getGPatients(self):
        return list(map(lambda x: int(float(x[2:])), self.patients))

    def getSamples(self):
        return self.samples

    def getTimes(self):
        return self.times

    def __repr__(self) -> str:
        return str(self.cluster)

    # def compute_MSR(self):
    # 	#genes-patients #contexts-features(samples) #times
    # 	r_s=[]
    # 	for t in self.cluster.keys():
    # 		#print("Processing Time: ", str(t))
    # 		for s in self.cluster[t].keys():
    # 			#print("Processing Feature: ", str(s))
    # 			for g in self.cluster[t][s].keys():
    # 				r_s.append(self.compute_R(g,s,t)**2)

        # return sum(r_s) / (int(self.nPatients) * int(self.nSamples) * int(self.nTimes))

    def getFeatValues(self, tp, feat):
        vals = filter(lambda x: x[0] == tp and x[1]
                      == feat, self.cluster.keys())

        return [self.cluster[v] for v in vals]

    def getPatientsVals(self, t):
        vals = filter(lambda x: x[0] == t, self.cluster.keys())
        return [self.cluster[v] for v in vals]

    def getSlice(self, g=None, c=None, t=None):
        coords = self.cluster
        if g != None and c != None and t != None:
            return coords[(t, c, g)]
        elif g != None and c != None:
            vals = filter(lambda x: x[2] == g and x[1] == c, coords.keys())
        elif g != None and t != None:
            vals = filter(lambda x: x[2] == g and x[2] == t, coords.keys())
        elif c != None and t != None:
            vals = filter(lambda x: x[1] == c and x[0] == t, coords.keys())
        elif g != None:
            vals = filter(lambda x: x[2] == g, coords.keys())
        elif c != None:
            vals = filter(lambda x: x[1] == c, coords.keys())
        elif t != None:
            vals = filter(lambda x: x[0] == t, coords.keys())
        else:
            return None
        return mygrouper(self.nSamples, [coords[v] for v in vals])

# 	def compute_R(self, g,c,t):
# 		TC_v = self.getSlice(g,c,t)
# 		M_ct = mean(self.getSlice(g=g))
# 		M_gt = mean(self.getSlice(c=c))
# 		M_gc = mean(self.getSlice(t=t))
# 		M_g = mean(self.getSlice(c=c, t=t))
# 		M_c = mean(self.getSlice(g=g, t=t))
# 		M_t = mean(self.getSlice(c=c, g=g))
# 		M_gct = mean(list(self.getTriclusterCoord().values()))

# 		return TC_v + M_ct + M_gt + M_gc - M_g - M_c - M_t - M_gct


def mygrouper(n, iterable):
    args = [iter(iterable)] * n
    return (list(([e for e in t if e != None] for t in itertools.zip_longest(*args))))


# if __name__ == "__main__":
# 	t = Tricluster(3,2,2)
# 	t.addPatient("P1")
# 	t.addPatient("P2")
# 	t.addSample("S1")
# 	t.addSample("S2")
# 	t.addTime("T1")
# 	t.addTime("T2")
# 	t.addTime("T3")
# 	t.addValue("T1", "S1", "P1", 3)
# 	t.addValue("T2", "S1", "P1", 3)
# 	t.addValue("T3", "S1", "P1", 3)
# 	t.addValue("T1", "S2", "P1", 4)
# 	t.addValue("T2", "S2", "P1", 5)
# 	t.addValue("T3", "S2", "P1", 6)
# 	t.addValue("T1", "S1", "P2", 3)
# 	t.addValue("T2", "S1", "P2", 3)
# 	t.addValue("T3", "S1", "P2", 3)
# 	t.addValue("T1", "S2", "P2", 4)
# 	t.addValue("T2", "S2", "P2", 5)
# 	t.addValue("T3", "S2", "P2", 6)

# 	import identify_triclusters as idt

# 	print(idt.compute_representative_pattern(t.getSlice(t="T3")))
