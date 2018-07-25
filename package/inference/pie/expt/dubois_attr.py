"""
[Python-Dev] Python interface to attribute descriptors 
Paul F Dubois paul@pfdubois.com 
Wed, 13 Nov 2002 14:08:57 -0800 

Question: how could the descriptor "know" the name "x" if it is created
by a descriptor-creating statement such as x = descriptor_creator(...).
I guess one could do this by making a metaclass that would look for the
descriptors in the class and "poke" the name into them but is there
another way? In the example below I evaded this question by making the
name an argument to positive's constructor.
"""

class positive (object):
    "Attribute that can only be positive."
    def __init__(self, name, default, doc=""):
        self.default = default
        self.__name__ = name
        self.__doc__ = doc

    def __get__ (self, obj, metaobj=None):
        if obj is None:
            return metaobj.__dict__[self.__name__]
        else:
            return getattr(obj, "__" + self.__name__, self.default)

    def __set__ (self, obj, value):
        try:
            v = type(self.default)(value)
            if not (v > 0):
                raise ValueError, "Value not positive."
        except:
            raise ValueError, "Cannot convert to positive %s" %
repr(type(self.default))
        setattr(obj, "__" + self.__name__, v)

class Simple(object) :
    """This class has two attributes, x and y, that must be
     positive floats.
    """
    x = positive("x", 1.0, "x-coordinate")
    y = positive("y", 2.0, "y-coordinate")

    def meth (self):
        return self.x * self.y

s = Simple()
print s.x, s.y, s.meth()
s.x, s.y = (2.,3)
print s.x, s.y, s.meth()
print Simple.x.__doc__
s2 = Simple()
print s2.x, s2.y, s2.meth()
try:
    s.x = -1.0
    raise RuntimeError, "Evaded validation"
except ValueError, e:
    print e

# When run this prints:
# 1.0 2.0 2.0
# 2.0 3.0 6.0
# x-coordinate
# 1.0 2.0 2.0
# Cannot convert to positive &lt;type 'float'>
