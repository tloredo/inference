import _obsinfo as f2pymodule

#FortranType = type(cbmlike.gammln)

# We want the doc strings of objects of type 'fortran', but
# this will do.

for name in dir(f2pymodule):
    obj = getattr(f2pymodule, name)
    try:
        doc = obj.__doc__
    except:
        continue
    # Ignore __doc__, __file__, __name__, etc.
    if type(obj) == str:
        continue
    print(':'*72)
    print(f'Name: {name}\n')
    print(obj.__doc__)
    print()
