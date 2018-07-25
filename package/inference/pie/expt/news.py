
class HoldString:
	x = 'ab cd'

hs1 = HoldString()
hs2 = HoldString()

print hs1.x
print hs1.x.split()
hs2.x = 'ef gh'
print hs1.x





class Parameter(object):
    """Insert named descriptor into calling instance."""
    def __init__(self, name):
        self.name = name + 'val'
    def __get__(self, inst, owner):
        # print 'getting val'
        return getattr(inst, self.name)
    def __set__(self, inst, val):
        # print 'setting val'
        setattr(inst, self.name, val)
    def __del__(self,inst):
        del inst.value


class test(object):
    x = Parameter('x')
    y = Parameter('y')

class mystring(str):
    pass

a=test()
a.x = 'ab cd'
print a.x
print a.x.split()
