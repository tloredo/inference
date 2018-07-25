
class test1:
    x = 1
    y = 2
    z = 3

if 0:
    a=test1()
    b=test1()
    print 'a:', a.x, a.y, a.z
    print 'b:', b.x, b.y, b.z
    a.x, a.y, a.z = 10, 20, 30  # inserts new floats into a's __dict__
    print 'a:', a.x, a.y, a.z
    print 'b:', b.x, b.y, b.z   # b keeps initialized values
    print 'a dict:', a.__dict__
    print 'b dict:', b.__dict__ # b's dict is empty; it uses test1.__dict__

# behavior unchanged with new classes, e.g., class test1(object).
print '----------------'

class test2(object):
    x = [1]
    y = [2]
    z = [3]

if 0:
    a=test2()
    b=test2()
    print 'a:', a.x, a.y, a.z
    print 'b:', b.x, b.y, b.z
    a.x[0], a.y[0], a.z[0] = 10, 20, 30  # inserts new floats into a's __dict__
    print 'a:', a.x, a.y, a.z
    print 'b:', b.x, b.y, b.z   # b shows the changes made to a
    print 'a dict:', a.__dict__ # a's dict is empty too this time; we modified
                                # values stored in test2.__dict__
    print 'b dict:', b.__dict__ # b's dict is empty; it uses test1.__dict__

# behavior unchanged with new classes, e.g., class test2(object).
print '----------------'

class test3(object):
    def getx(self):
        return self.xval    # be sure not to use self.x here!
    def setx(self, val):
        self.xval = val
    def delx(self):
        del self.xval
    x = property(getx, setx, delx, 'x attribute')

if 1:
    print 'test3...'
    a=test3()
    b=test3()
    a.x = 1
    b.x = 2
    print 'a.x=',a.x,'  b.x=',b.x
    # Prints "a.x= 1   b.x= 2" -- good, we have dif values across instances.

    print '----------------'

# Now we want to have different instances of something like x *within*
# a single class, each maintaining a different value.  We'd like to
# define x outside the containing class, with a descriptor.

class Parameter(object):
    """A parameter descriptor that stores the param value in the 
    descriptor dict."""
    def __init__(self, default=None):
        self.value = default
    def __get__(self, inst, owner):
        # print 'getting val'
        return self.value
    def __set__(self, inst, val):
        # print 'setting val'
        self.value = val
    def __delete__(self,inst):
        del self.value

class test4(object):
    x = Parameter(1)

if 0:
    print 'test4...'
    a=test4()
    b=test4()
    a.x = 1
    b.x = 2
    print 'a.x=',a.x,'  b.x=',b.x
    # prints "a.x= 2   b.x= 2" since the b assignment changes self.value
    # in the (mutable!) Parameter instance stored in test4.__dict__
    print a.__dict__    # empty!
    print test4.__dict__    # includes x as a Parameter instance

print '----------------'

# We need to do something like test3, storing state for each param
# with the class *instance* dict, not in the class dict.

class Parameter2(object):
    def __get__(self, inst, owner):
        # print 'getting val'
        return inst.value    # note use of *inst* here
    def __set__(self, inst, val):
        # print 'setting val'
        inst.value = val
    def __delete__(self,inst):
        del inst.value

class test5(object):
    x = Parameter2()
    y = Parameter2()

if 0:
    print 'test5...'
    a=test5()
    b=test5()
    a.x = 1
    b.x = 2
    print 'a.x=',a.x,'  b.x=',b.x
    # prints "a.x= 1   b.x= 2" since the assignments store values in a.__dict__
    # and b.__dict__.  But this mechanism can't have multiple Parameter2 instances
    # in a class:
    b.y = 3
    print 'b.x, b.y:', b.x, b.y #   prints "3 3"!!
    print 'a dict:', a.__dict__     # contains value=1
    print 'b dict:', b.__dict__     # contains value=2
    print test4.__dict__    # includes x as a Parameter instance but no "value"

    # We need a way for Parameter2 to use different names for each instance.

    print '----------------'

class Parameter3(object):
    """Simpler version of Paul Dubois's Python-Dev posting."""
    def __init__(self, name):
        self.name = name + 'val'    # don't use just "name" here---recursion!
    def __get__(self, inst, owner):
        # print 'getting val'
        return getattr(inst, self.name)
    def __set__(self, inst, val):
        # print 'setting val'
        setattr(inst, self.name, val)
    def __delete__(self,inst):
        del inst.value

class test6(object):
    x = Parameter3('x')
    y = Parameter3('y')

print 'test6...'
a=test6()
b=test6()
a.x = 1
b.x = 2
print 'a.x=',a.x,'  b.x=',b.x
# prints "a.x= 1   b.x= 2" since the assignments store values in a.__dict__
# and b.__dict__.  This mechanism *allows* multiple Parameter3 instances
# in a class:
b.y = 3
print 'b.x, b.y:', b.x, b.y #   prints "2 3" -- yeah!!
print 'a dict:', a.__dict__     # contains xval=1
print 'b dict:', b.__dict__     # contains xval=2, yval=3

# But it's awkward to have to pass the names of x, y when they are right there.

print '----------------'

# Based on Python-Dev posting by Phillib Eby.

class metaAutoName(type):
    """Instances of this type will be notified of the name of any attributes
    that are instances of AutoNamed."""

    def __init__(cls, name, bases, dict):
        for key,val in dict.items():
            if isinstance(val, AutoNamed):
                val.notify(cls, key)  # Tell the class "key" is an AutoNamed name.
        super(metaAutoName, cls).__init__(name,bases,dict)

class AutoNamed(object):
    """An abstract base class intended for automatically passing name info
    to descriptors.  Subclasses should override "notify" to do
    something useful with the name of the subclass instance."""

    def notify(self, cls, name):
        """Override this method to do something useful with name."""
        raise NotImplementedError

# Note:  AutoName cannot appear before AutoNamed; the interpreter will
# complain that AutoNamed is not defined when it uses metaAutoName.

class HasAutoNamed(object):
    """Empty base class for classes with AutoNamed descriptors."""
    __metaclass__ = metaAutoName

class Parameter4(AutoNamed):

    def notify(self, cls, name):
        print 'Parameter4 instance in',cls,'is named',name
        self.name = name + 'val'    # don't use just "name" here---recursion!
    def __get__(self, inst, owner):
        # print 'getting val'
        return getattr(inst, self.name)
    def __set__(self, inst, val):
        # print 'setting val'
        setattr(inst, self.name, val)
    def __delete__(self,inst):
        del inst.value

class test6(HasAutoNamed):
    x = Parameter4()
    y = Parameter4()

print 'test6...'
a=test6()
b=test6()
a.x = 1
b.x = 2
print 'a.x=',a.x,'  b.x=',b.x
# prints "a.x= 1   b.x= 2" since the assignments store values in a.__dict__
# and b.__dict__.  This mechanism *allows* multiple Parameter3 instances
# in a class:
b.y = 3
print 'b.x, b.y:', b.x, b.y #   prints "2 3" -- yeah!!
print 'a dict:', a.__dict__     # contains xval=1
print 'b dict:', b.__dict__     # contains xval=2, yval=3
# But now the names "xval" and "yval" are created automatically!!!

#----------------------------------------------------------------
# Can descriptors have "normal" methods?

class Parameter5(AutoNamed):

    def notify(self, cls, name):
        print 'Parameter5 instance in',cls,'is named',name
        self.name = name + 'val'    # don't use just "name" here---recursion!
    def check(self):
        print 'Calling check method'
    def __get__(self, inst, owner):
        # print 'getting val'
        return getattr(inst, self.name)
    def __set__(self, inst, val):
        # print 'setting val'
        setattr(inst, self.name, val)
    def __delete__(self,inst):
        del inst.value

class test7(HasAutoNamed):
    x = Parameter5()
    y = Parameter5()

print 'test7...'
a=test7()
b=test7()
a.x = 1
b.x = 2
print 'a.x=',a.x,'  b.x=',b.x
a.x.check()
# This fails because a.x gets returned by __get__ as just an int;
# no .check method exists!
# To get around this, instead of using an int, we'll derive a class
# from an int and capture its non-int methods to send to a handler
# (or perhaps just add non-int methods directly---which has better
# performance?).
