import warnings
from glob import glob
import numpy as np

#from six import  string_types
"""Utilities for writing code that runs on Python 2 and 3"""

# Copyright (c) 2010-2015 Benjamin Peterson
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#from __future__ import absolute_import

#import functools
#import itertools
#import operator
import sys
import types

__author__ = "Benjamin Peterson <benjamin@python.org>"
__version__ = "1.9.0"


# Useful for very coarse version differentiation.
PY2 = sys.version_info[0] == 2
PY3 = sys.version_info[0] == 3

if PY3:
    string_types = str,
    integer_types = int,
    class_types = type,
    text_type = str
    binary_type = bytes

    MAXSIZE = sys.maxsize
else:
    string_types = basestring,
    integer_types = (int, long)
    class_types = (type, types.ClassType)
    text_type = unicode
    binary_type = str

    if sys.platform.startswith("java"):
        # Jython always uses 32 bits.
        MAXSIZE = int((1 << 31) - 1)
    else:
        # It's possible to have sizeof(long) != sizeof(Py_ssize_t).
        class X(object):
            def __len__(self):
                return 1 << 31
        try:
            len(X())
        except OverflowError:
            # 32-bit
            MAXSIZE = int((1 << 31) - 1)
        else:
            # 64-bit
            MAXSIZE = int((1 << 63) - 1)
        del X


def _glob(filenames):
    """Glob a filename or list of filenames but always return the original
    string if the glob didn't match anything so URLs for remote file access
    are not clobbered.
    """
    if isinstance(filenames, string_types):
        filenames = [filenames]
    matches = []
    for name in filenames:
        matched_names = glob(name)
        if not matched_names:
            # use the original string
            matches.append(name)
        else:
            matches.extend(matched_names)
    return matches


def tree2array(tree,
               branches=None,
               selection=None,
               object_selection=None,
               start=None,
               stop=None,
               step=None,
               include_weight=False,
               weight_name='weight',
               cache_size=-1):
    """Convert a tree into a numpy structured array.

    Convert branches of strings and basic types such as bool, int, float,
    double, etc. as well as variable-length and fixed-length multidimensional
    arrays and 1D or 2D vectors of basic types and strings. ``tree2array`` can
    also create columns in the output array that are expressions involving the
    TTree branches (i.e. ``'vect.Pt() / 1000'``) similar to ``TTree::Draw()``.
    See the notes below for important details.

    Parameters
    ----------
    tree : ROOT TTree instance
        The ROOT TTree to convert into an array.
    branches : list of strings and tuples or a string or tuple, optional (default=None)
        List of branches and expressions to include as columns of the array or
        a single branch or expression in which case a nonstructured array is
        returned. If None then include all branches that can be converted.
        Branches or expressions that result in variable-length subarrays can be
        truncated at a fixed length by using the tuple ``(branch_or_expression,
        fill_value, length)`` or converted into a single value with
        ``(branch_or_expression, fill_value)`` where ``length==1`` is implied.
        ``fill_value`` is used when the original array is shorter than
        ``length``. This truncation is after any object selection performed
        with the ``object_selection`` argument.
    selection : str, optional (default=None)
        Only include entries fulfilling this condition. If the condition
        evaluates to multiple values per tree entry (e.g. conditions on array
        branches) then an entry will be included if the condition evaluates to
        true for at least one array element.
    object_selection : dict, optional (default=None)
        A dictionary mapping selection strings to branch names or lists of
        branch names. Only array elements passing the selection strings will be
        included in the output array per entry in the tree. The branches
        specified must be variable-length array-type branches and the length of
        the selection and branches it acts on must match for each tree entry.
        For example ``object_selection={'a > 0': ['a', 'b']}`` will include all
        elements of 'a' and corresponding elements of 'b' where 'a > 0' for
        each tree entry. 'a' and 'b' must have the same length in every tree
        entry.
    start, stop, step: int, optional (default=None)
        The meaning of the ``start``, ``stop`` and ``step`` parameters is the
        same as for Python slices. If a range is supplied (by setting some of
        the ``start``, ``stop`` or ``step`` parameters), only the entries in
        that range and fulfilling the ``selection`` condition (if defined) are
        used.
    include_weight : bool, optional (default=False)
        Include a column containing the tree weight ``TTree::GetWeight()``.
        Note that this will be the same value for all entries unless the tree
        is actually a TChain containing multiple trees with different weights.
    weight_name : str, optional (default='weight')
        The field name for the weight column if ``include_weight=True``.
    cache_size : int, optional (default=-1)
        Set the size (in bytes) of the TTreeCache used while reading a TTree. A
        value of -1 uses ROOT's default cache size. A value of 0 disables the
        cache.

    Notes
    -----
    Types are converted according to the following table:

    .. _conversion_table:

    ========================  ===============================
    ROOT                      NumPy
    ========================  ===============================
    ``Bool_t``                ``np.bool``
    ``Char_t``                ``np.int8``
    ``UChar_t``               ``np.uint8``
    ``Short_t``               ``np.int16``
    ``UShort_t``              ``np.uint16``
    ``Int_t``                 ``np.int32``
    ``UInt_t``                ``np.uint32``
    ``Float_t``               ``np.float32``
    ``Double_t``              ``np.float64``
    ``Long64_t``              ``np.int64``
    ``ULong64_t``             ``np.uint64``
    ``<type>[2][3]...``       ``(<nptype>, (2, 3, ...))``
    ``<type>[nx][2]...``      ``np.object``
    ``string``                ``np.object``
    ``vector<t>``             ``np.object``
    ``vector<vector<t> >``    ``np.object``
    ========================  ===============================

    * Variable-length arrays (such as ``x[nx][2]``) and vectors (such as
      ``vector<int>``) are converted to NumPy arrays of the corresponding
      types.

    * Fixed-length arrays are converted to fixed-length NumPy array fields.

    **Branches with different lengths:**

    Note that when converting trees that have branches of different lengths
    into numpy arrays, the shorter branches will be extended to match the
    length of the longest branch by repeating their last values. If all
    requested branches are shorter than the longest branch in the tree, this
    will result in a "read failure" since beyond the end of the longest
    requested branch no additional bytes will be read from the file and
    root_numpy is unable to distinguish this from other ROOT errors that result
    in no bytes being read. In this case, explicitly set the ``stop`` argument
    to the length of the longest requested branch.


    See Also
    --------
    root2array
    array2root
    array2tree

    """
    import ROOT
    if not isinstance(tree, ROOT.TTree):
        raise TypeError("tree must be a ROOT.TTree")

    if isinstance(branches, string_types):
        # single branch selected
        flatten = branches
        branches = [branches]
    elif isinstance(branches, tuple):
        if len(branches) not in (2, 3):
            raise ValueError(
                "invalid branch tuple: {0}. "
                "A branch tuple must contain two elements "
                "(branch_name, fill_value) or three elements "
                "(branch_name, fill_value, length) "
                "to yield a single value or truncate, respectively".format(branches))
        flatten = branches[0]
        branches = [branches]
    else:
        flatten = False

    copytree = tree.CopyTree(selection)
    arr = copytree.AsMatrix(columns=branches )
    #arr = _librootnumpy.root2array_fromtree(
    #    cobj, branches, selection, object_selection,
    #    start, stop, step,
    #    include_weight,
    #    weight_name,
    #    cache_size)

    if flatten:
        # select single column
        return arr[flatten]
    return arr



def array2tree(arr, name='tree', tree=None):
    """Convert a numpy structured array into a ROOT TTree.

    Fields of basic types, strings, and fixed-size subarrays of basic types are
    supported. ``np.object`` and ``np.float16`` are currently not supported.

    Parameters
    ----------
    arr : array
        A numpy structured array
    name : str (optional, default='tree')
        Name of the created ROOT TTree if ``tree`` is None.
    tree : ROOT TTree (optional, default=None)
        An existing ROOT TTree to be extended by the numpy array. Any branch
        with the same name as a field in the numpy array will be extended as
        long as the types are compatible, otherwise a TypeError is raised. New
        branches will be created and filled for all new fields.

    Returns
    -------
    root_tree : a ROOT TTree

    Notes
    -----
    When using the ``tree`` argument to extend and/or add new branches to an
    existing tree, note that it is possible to create branches of different
    lengths. This will result in a warning from ROOT when root_numpy calls the
    tree's ``SetEntries()`` method. Beyond that, the tree should still be
    usable. While it might not be generally recommended to create branches with
    differing lengths, this behaviour could be required in certain situations.
    root_numpy makes no attempt to prevent such behaviour as this would be more
    strict than ROOT itself. Also see the note about converting trees that have
    branches of different lengths into numpy arrays in the documentation of
    :func:`tree2array`.

    See Also
    --------
    array2root
    root2array
    tree2array

    Examples
    --------

    Convert a numpy array into a tree:

    >>> from root_numpy import array2tree
    >>> import numpy as np
    >>>
    >>> a = np.array([(1, 2.5, 3.4),
    ...               (4, 5, 6.8)],
    ...              dtype=[('a', np.int32),
    ...                     ('b', np.float32),
    ...                     ('c', np.float64)])
    >>> tree = array2tree(a)
    >>> tree.Scan()
    ************************************************
    *    Row   *         a *         b *         c *
    ************************************************
    *        0 *         1 *       2.5 *       3.4 *
    *        1 *         4 *         5 *       6.8 *
    ************************************************

    Add new branches to an existing tree (continuing from the example above):

    >>> b = np.array([(4, 10),
    ...               (3, 5)],
    ...              dtype=[('d', np.int32),
    ...                     ('e', np.int32)])
    >>> array2tree(b, tree=tree)
    <ROOT.TTree object ("tree") at 0x1449970>
    >>> tree.Scan()
    ************************************************************************
    *    Row   *         a *         b *         c *         d *         e *
    ************************************************************************
    *        0 *         1 *       2.5 *       3.4 *         4 *        10 *
    *        1 *         4 *         5 *       6.8 *         3 *         5 *
    ************************************************************************

    """
    import ROOT
    if tree is not None:
        if not isinstance(tree, ROOT.TTree):
            raise TypeError("tree must be a ROOT.TTree")
    else:
        tree  = ROOT.TTree(name, name)

    if arr.dtype.names is None:
        raise NameError("array.dtype.names does not exist")

    br_description = ":".join(name for name in a.dtype.names)
    #for i in range(len(arr.dtype.names)):
    #    branchname = arr.dtype.names[i]
    #    dtype = arr.dtype.fields[branchname]
    arr_str = "\n".join("\t".join( str(x) for x in row[0]) for row in arr)
    istring        = ROOT.istringstream(arr_str)
    tree.ReadStream(istring, br_description)
    #tree.Scan()
    return tree




#########################


import ROOT

tree = ROOT.TChain("t")

tfile = " /Users/taohuang/Documents/DiHiggs/20180126/20180205_20180202_10k_Louvain_ALL/radion_M350_all.root"
#tree.Add("ST_tW_antitop_5f_86887FC7-9313-E811-AC4A-0CC47A13CDA0.root")
tree.Add(tfile)

#branches = [""]
cuts = "isSF>0 && jj_pt>100"
cuts = "1"

FullInputs = ['jj_M','jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j','isSF']
#FullInputs = ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula','mt2', 'jj_M','isSF','hme_h2mass_reco','isElEl','isElMu','isMuEl','isMuMu','met_pt']
#FullInputs = ['jj_pt', 'll_pt', 'll_M', 'll_DR_l_l', 'jj_DR_j_j', 'llmetjj_DPhi_ll_jj', 'llmetjj_minDR_l_j', 'llmetjj_MTformula','mt2', 'jj_M','hme_h2mass_reco','isElEl','isElMu','isMuEl','isMuMu','met_pt']
#FullInputs = ["isSF"]

#print("Input variables ", FullInputs)

tree.Scan()
a = tree2array(tree, branches = FullInputs, selection = cuts)
shape = a.shape
dt = np.dtype({'names': FullInputs, 'formats':[np.double]*len(FullInputs)})
#dt = {'names': FullInputs}
#print("dt ",dt)
a.dtype = dt

#print("old shape ",shape, " after dtype assign ", a.shape, " dtype ", a.dtype, " names ", a.dtype.names, " field ", len(a.dtype.names))
#print "print a type ", type(a)," shape ",a.shape, " len(a) ",len(a), " a.dtype.names ",a.dtype.names, " a[:10,:] ", a[:10]," a.dtype.type ", a.dtype.type

for i in range(len(a.dtype.names)):
    branchname = a.dtype.names[i]
    dtype = a.dtype.fields[branchname]
    #print("i ",i," branchname ", branchname, " type ",dtype)
br_description = ":".join(name for name in a.dtype.names)
#print " a.dtype.names ",a.dtype.names, " isSF ", a["isSF"][:10]

for row in a[:10]:
    print " row ",str(row[0]), "\t".join( str(x) for x in row[0])

#arr_str = np.tostring(a[:10])
#arr_str = "\n".join("\t".join( str(x) for x in row[0]) for row in a[:10] )
#istring        = ROOT.istringstream(arr_str)
#print "array into string type ",type(arr_str), " arr_str ", arr_str
#testtree = ROOT.TTree("test","test")
#testtree.ReadStream(istring, br_description)
testtree = array2tree(a, "outtree") 
testtree.Scan()
#print("Only the content of the branch 'x':\n{}\n".format(np.squeeze(arr)))
