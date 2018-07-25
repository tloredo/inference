import os
import os.path as path
import pickle as pkl
from numpy import ndarray, savez, load

__all__ = ['backup', 'safe_open', 'ArrayStore', 'AttrStore']

def backup(fname):
    """
    Backup any existing file with name fname, up to two levels.  If a
    file fname exists, the existing copy is renamed with a '%'appended
    to its name.  If such a backup already exists, it is renamed with an
    additional '%' appended, overwriting any previous such file.
    """
    if path.exists(fname):
        bname = fname + '%'
        if path.exists(bname):
            bbname = bname + '%'
            os.rename(bname, bbname)
        os.rename(fname, bname)

def safe_open(fname, binary=False, opener=None, *args, **kwds):
    """
    Open a file by name for writing, backing up an existing version up
    to two levels.  If the file exists, the existing copy is renamed with
    a '%'appended to its name.  If such a backup already exists, it is
    renamed with an additional '%' added, overwriting any previous such
    file.  Return the opened file object.

    By default, the file is opened with the `open` built-in, with 'w'
    mode.  If `opener` is specified, it is used to open the file,
    passing any additional arguments after the file name.
    """
    backup(fname)
    if opener:
        if args:
            fobj = opener(fname, *args, **kwds)
        else:
            fobj = opener(fname)
    else:
        if binary is True or binary == 'b':
            fobj = open(fname, 'wb')
        else:
            fobj = open(fname, 'w')
    return fobj

def read_column_vectors(dfile, types, as_arrays=True):
    """
    Read vectors of scalar-valued data stored as whitespace-delimited strings
    in columns in an open file.  Return a list of vectors storing the values
    read from each column.
    
    The file is presumed open so that header lines can be discarded as needed
    before calling read_vectors().
    
    `types` is a list of callables, one for each column, that convert a
    string representation of a scalar into the actual value.
    
    If `as_arrays` is True, the returned values are Numpy arrays.
    
    If all columns are of the same type, similar capability is available from
    numpy.loadtxt with unpack=True, or by reading rows as records and
    extracting vectors from the resulting array of records.
    """

    # TODO:  Replace with numpy.loadtxt?
    # TODO:  Add skiprows, usecols as in loadtxt?
    # TODO:  Use csv; add delimiter as in loadtxt?

    nc = len(types)
    cols = [[] for i in range(nc)]  # list of nc empty lists
    for line in dfile.readlines():
        words = line.strip().split()
        for i in range(nc):
            cols[i].append(types[i](words[i]))
    if as_arrays:
        for i in range(nc):
            cols[i] = array(cols[i])
    return [cols[i] for i in range(nc)]


class ArrayStore(object):
    """
    A container class that provides lazy persistance of array instance
    attributes via NumPy's .npz archives.

    This is doubly lazy:  arrays are not read into the namespace until
    actually requested, and they are not saved to the archive until the
    ArrayStore is explicitly saved.
    """

    # *** Is loading really lazy, i.e., does NumPy's load() read in array data
    # before actual access via dict lookup?

    # TODO:  Perhaps have a __del__ method that saves on deletion?
    # It would have to handle fname=None gracefully.
    # Perhaps an init arg should specify whether to save on deletion.

    # TODO:  Support creating with a name to a non-existing file;
    # subsequent "save()" will use the stored name.

    def __init__(self, fname=None):
        """
        Prepare to load arrays from storage if a file name is provided;
        otherwise support saving of arrays assigned as attributes.
        """
        # All internal attributes start with '_' to avoid __setattr__
        # array filtering.
        if fname:
            if not fname.endswith('.npz'):
                fname = fname + '.npz'
                self._npz = load(fname)
        else:
            self._npz = None
        self._fname = fname
        self._archived = {}  # arrays pulled from the archive
        self._new = {}  # arrays to be archived

    def __getattr__(self, name):
        """
        Catch references to array attributes that have not yet been
        loaded from the archive.
        """
        # This is called only if name is not already in the instance dict.
        if self._npz is None:  # no archive to grab attribute from
            raise AttributeError(name)
        else:  # get value from archive and keep a reference to it
            try:
                value = self._npz[name]
                # Set the attribute directly so __setattr__ won't add it
                # to self.new.
                object.__setattr__(self, name, value)
                self._archived[name] = value
                return value
            except:
                raise AttributeError(name)

    def __setattr__(self, name, value):
        """
        Catch assignments to new attributes, marking them for saving when
        the store is next saved.

        Names starting with '_' have their corresponding attributes set
        without marking; such names are intended for internal use only,
        to keep track of the state of the store.
        """
        # *** Should this prevent over-writing existing names, either
        # directly assigned or yet to be loaded from the archive?
        # Right now we allow reassignments; saving will save the
        # reassigned value.
        if name.startswith('_'):
            object.__setattr__(self, name, value)
        elif isinstance(value, ndarray):
            self._new[name] = value
            object.__setattr__(self, name, value)
        else:
            raise ValueError('Only ndarray objects may be stored!')

    def save(self, fname=None):
        """
        Save array attributes to a NumPy .npz archive.  If a name is
        provided, it is used (adding '.npz' if needed); otherwise it is
        presumed this store was created from an existing archive, and
        that archive's name is used.

        In either case, any existing version is backed up before the
        new one is created, with up to two levels of backups (with '%'
        and '%%' suffixes).
        """
        # *** Perhaps this should be called by __del__; or this possibility
        # could be set by a flag on init???

        # If no name given, use the name provided at creation.
        if fname is None:
            if self._fname is None:
                raise ValueError('Need a file name for saving!')
            fname = self._fname
        else:
            if not fname.endswith('.npz'):
                fname = fname + '.npz'
        # Gather arrays to store; note new versions supersede old ones.
        arrays = {}
        for name, val in list(self._archived.items()):
            arrays[name] = val
        # Don't forget to get any archived values not already accessed.
        if self._npz:
            for name in self._npz.files:
                if name not in self._archived:
                    arrays[name] = self._npz[name]
        for name, val in list(self._new.items()):
            arrays[name] = val
        # Backup any existing file, and save the arrays to a new file.
        backup(fname)
        savez(fname, **arrays)

    def contents(self):
        """
        List the names of arrays accessible by this store, including
        both previously archived arrays and newly defined arrays.
        """
        names = []
        if self._npz:
            names.extend(self._npz.files)  # archive file names will be array names
        names.extend(list(self._new.values()))
        return names


class AttrStore(object):
    """
    A container class that provides lazy persistance of instance attributes
    as a Python pickle.

    This is doubly lazy:  attributes are not placed in the namespace until
    actually requested, and they are not saved to the archive until the
    AttrStore is explicitly saved.
    
    If only NumPy arrays are to be stored, ArrayStore may be more efficient.
    """

    # *** Is there any virtue to laziness in loading the namespace, since
    # we aren't unpacking a NumPy npz file archive?

    # TODO:  Perhaps have a __del__ method that saves on deletion?
    # It would have to handle fname=None gracefully.
    # Perhaps an init arg should specify whether to save on deletion.

    # TODO:  Support creating with a name to a non-existing file;
    # subsequent "save()" will use the stored name.
    
    # TODO:  Add kwd args to init for setting attributes to store.

    def __init__(self, fname=None, **kwds):
        """
        Prepare to load attributes from storage if a file name is provided;
        otherwise support saving of arrays assigned as attributes.
        """
        # All internal attributes start with '_' to avoid __setattr__
        # array filtering.
        if fname:
            if not fname.endswith('.pkl'):
                fname = fname + '.pkl'
                ifile = open(fname, 'rb')
                self._archive = pkl.load(ifile)
                ifile.close()
        else:
            self._archive = None
        self._fname = fname
        self._pulled = {}  # attributes pulled from the archive
        self._new = {}  # attributes to be archived
        if kwds:
            for key, value in list(kwds.items()):
                setattr(self, key, value)

    def __getattr__(self, name):
        """
        Catch references to attributes that have not yet been loaded from
        the store.
        """
        # This is called only if name is not already in the instance dict.
        if self._archive is None:  # no archive to grab attribute from
            raise AttributeError(name)
        else:  # get value from archive and keep a reference to it
            try:
                value = self._archive[name]
                # Set the attribute directly so __setattr__ won't add it
                # to self.new.
                object.__setattr__(self, name, value)
                self._pulled[name] = value
                return value
            except:
                raise AttributeError(name)

    def __setattr__(self, name, value):
        """
        Catch assignments to new attributes, marking them for saving when
        the store is next saved.

        Names starting with '_' have their corresponding attributes set
        without marking; such names are intended for internal use only,
        to keep track of the state of the store.
        """
        # *** Should this prevent over-writing existing names, either
        # directly assigned or yet to be loaded from the archive?
        # Right now we allow reassignments; saving will save the
        # reassigned value.
        if name.startswith('_'):
            object.__setattr__(self, name, value)
        else:
            self._new[name] = value
            object.__setattr__(self, name, value)

    def save(self, fname=None):
        """
        Save attributes to a high-protocol pickle.  If a name is
        provided, it is used (adding '.pkl' if needed); otherwise it is
        presumed this store was created from an existing archive, and
        that archive's name is used.

        In either case, any existing version is backed up before the
        new one is created, with up to two levels of backups (with '%'
        and '%%' suffixes).
        """
        # *** Perhaps this should be called by __del__; or this possibility
        # could be set by a flag on init???

        # If no name given, use the name provided at creation.
        if fname is None:
            if self._fname is None:
                raise ValueError('Need a file name for saving!')
            fname = self._fname
        else:
            if not fname.endswith('.pkl'):
                fname = fname + '.pkl'
        # Gather attributes to store; note new versions supersede old ones.
        attrs = {}
        for name, val in list(self._pulled.items()):
            attrs[name] = val
        # Don't forget to get any archived values not already accessed.
        if self._archive:
            for name in list(self._archive.keys()):
                if name not in self._pulled:
                    attrs[name] = self._archive[name]
        for name, val in list(self._new.items()):
            attrs[name] = val
        # Backup any existing file, and save the attrs to a new file.
        ofile = safe_open(fname, 'b')
        pkl.dump(attrs, ofile, -1)  # use highest protocol
        ofile.close()

    def contents(self):
        """
        List the names of arrays accessible by this store, including
        both previously archived arrays and newly defined arrays.
        """
        names = []
        if self._archive:
            names.extend(list(self._archive.keys()))
        names.extend(list(self._new.values()))
        return names



if __name__ == '__main__':
    from numpy import *

    # Store some arrays.
    print('Testing ArrayStore:')
    store = ArrayStore()
    store.a = array([1,2,3])
    store.b = array([[1., 2.], [3., 4.]])
    store.save('junk')

    # Read the stored arrays.
    store2 = ArrayStore('junk')
    print('a:', store2.a)
    print('b:', store2.b)
    store2.c = eye(6)
    store2.a = ones(5, float)
    store2.save()
    store2 = ArrayStore('junk')
    print('store2 should have altered a and new c.')

    # Store a mix of arrays and other pickleable objects.
    print('\nTesting AttrStore:')
    store3 = AttrStore(a=array([1,2,3]), c=42)
    store3.b = array([[1., 2.], [3., 4.]])
    store3.d = ['a', 'b', 'c']
    store3.s = 'This is a test'
    store3.save('junk2')

    # Read the stored objects.
    store4 = AttrStore('junk2')
    print('a:', store4.a)
    print('b:', store4.b)
    print('c:', store4.c)
    print('d:', store4.d)
    print('s:', store4.s)
    store4.e = eye(6)
    store4.a = ones(5, float)
    store4.save()
    store4 = AttrStore('junk2')
    print('store4 should have altered a and new e.')
