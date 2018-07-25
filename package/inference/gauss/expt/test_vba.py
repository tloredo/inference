from numpy import array, ones, transpose, dot, sqrt
from inference.gauss import vba


nsamp = 10
v1 = ones(nsamp, float)
v2 = ones(nsamp,float)
data = 0.5*ones(nsamp,float)
for i in range(nsamp):
    if (i % 2 == 0):
        v2[i] = -v2[i]
v2[4] = 2.
v2[6] = 2.
dmat = array([v1,v2])
data = ones(nsamp, float)
data[4] = -1.
data[6] = -1.
print dmat
print data

def results(dmat,data):
    (metric, L, J, P, A, S) = vba.metricAnalysis(dmat,data)

    print 'Metric: ', metric
    print 'L: ', L
    print 'J: ', J
    print sqrt(metric[0,0]*metric[1,1] - metric[1,0]*metric[0,1])
    U = transpose(L)
    print 'LU: ', dot(L,U)
    print 'P: ', P
    print 'A: ', A
    print 'S: ', S
    C = vba.covar(L)
    print 'C: ', C
    print 'C*metric: ', dot(C,metric)

results(dmat,data)
print '----------------------'
v3 = ones(nsamp,float)
for i in range(nsamp):
    if (i >=5):
        v3[i] = -v3[i]
dmat = array([v1,v2,v3])
print dmat
print data
results(dmat,data)
