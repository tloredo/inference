from numpy import array, load

ifile = open('lmc.dat','r')

lines = ifile.readlines()[1:]

vals = []
for line in lines:
    vals.append(float(line.strip()))

vals = array(vals)

print 'Found', len(vals), 'arrival times'

ifile.close()

vals.dump(open('lmc.pkl','w'))

v = load('lmc.pkl')
