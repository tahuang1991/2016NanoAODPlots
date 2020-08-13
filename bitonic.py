def bitonic_sort(up, x):
    if len(x) <= 1:
        return x
    else: 
        first = bitonic_sort(True, x[:len(x) // 2])
        second = bitonic_sort(False, x[len(x) // 2:])
        return bitonic_merge(up, first + second)

def bitonic_merge(up, x): 
    # assume input x is bitonic, and sorted list is returned 
    if len(x) == 1:
        return x
    else:
        bitonic_compare(up, x)
        first = bitonic_merge(up, x[:len(x) // 2])
        second = bitonic_merge(up, x[len(x) // 2:])
        return first + second

def bitonic_compare(up, x):
    dist = len(x) // 2
    print "bitonic_compare ",x," dist ",dist," order ", "increasing " if up else "decreasing "
    for i in range(dist):  
        if (x[i] > x[i + dist]) == up:
            x[i], x[i + dist] = x[i + dist], x[i] #swap
    print "results ",x

inputlist = [8,7,6,5,4,3,2,1]
print "inputlist ",inputlist
print "outputlist increasing ",bitonic_sort(True, inputlist)
print "outputlist decreasing ",bitonic_sort(False, inputlist)
