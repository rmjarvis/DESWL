import numpy as np

def match(id1, id2):
    """This is equivalent to Erin's esutil.numpy_util.match, except less efficient.
    """

    arg1 = np.argsort(id1)
    arg2 = np.argsort(id2)
    # Now id1[arg1] is sorted and id2[arg2] is sorted.
    use1 = np.in1d(id1[arg1], id2[arg2])
    use2 = np.in1d(id2[arg2], id1[arg1])

    index1 = arg1[use1]
    index2 = arg2[use2]
    return index1, index2

