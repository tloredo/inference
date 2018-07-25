"""
Automatically named descriptors

The following three classes comprise a type (a metaclass) and two base
classes for implementing classes with descriptors that automatically get
notified of the name of the attribute they are assigned to. The intent is
for the descriptor instance to use this name to facilitate storing state in an instance
of the class containing the descriptor.  This way a class may contain multiple
instances of the same descriptor (as distinctly named attributes), each with
its own state and behavior in every instance of the containing class.

The motivation is to allow a class to represent a parameterized model
that may have a number of parameters of the same type (e.g., real scalars).
These parameters could not only have different values, but also different
allowed ranges (i.e., state) and different behavior (i.e., some fixed
and some varying or stepped as the model is fit to data).

Useage in a nutshell:  If a class inherting "HasAutoNamed" has descriptors
inheriting "AutoNamed", the descriptor instances will be notified of their
attribute names via calls to their "_report_name" methods.

This implementation is based on similar material in Python-Dev list posts
ca. the release of Python 2.2.
"""


class metaAutoName(type):
    """
    A metaclass for classes with descriptors automatically notified of thier names.

    Instances of this type can have AutoName'd descriptors that get notified of the
    names of the attributes they get assigned to.
    """

    def __init__(cls, name, bases, dict):
        for key,val in dict.items():
            # Search for AutoNamed descriptors.
            if isinstance(val, AutoNamed):
                # Tell the descriptor "key" is its name.
                val._report_name(cls, key)
        super(metaAutoName, cls).__init__(name,bases,dict)


class AutoNamed(object):
    """
    An abstract base class intended for automatically passing name info
    to descriptors.  Subclasses should override "_report_name" to do
    something useful with the name of the subclass instance.
    """

    def _report_name(self, cls, name):
        """
        Override this method to do something useful with name.
        """
        raise NotImplementedError('Class must implement _reportname!')


# Note:  HasAutoNamed cannot appear before AutoNamed; the interpreter will
# complain that AutoNamed is not defined when it uses metaAutoName to construct
# the HasAutoNamed class.

class HasAutoNamed(object):
    """
    Empty base class that acts as a type for classes with AutoNamed
    descriptors.
    """
    __metaclass__ = metaAutoName

