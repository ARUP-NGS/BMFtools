'''
Code snippets working on uniqueness analysis.


singletons = [i for i in b if i.qname in set(key for key, val in counts.iteritems() if val == 1)]
pairs = list(chunk(sorted([i for i in b if i.qname in set(key for key, val in counts.iteritems() if val == 2)], key=lambda x: x.qname), 2))

def get_olap(read, start, stop):                                                                                   
    return sum([pos >= start and pos < stop for pos in read.get_reference_positions()])

sorted(freq(map(lambda x: pair_olap(x, "10", 89624206, 89624325), pairs)).iteritems(), key = lambda x: x[0])
overlaps = [(pair_olap(x, "10", 89624206, 89624325), x) for x in pairs]

[of.write(i) for i in chain.from_iterable([j for i, j in overlaps if i == 0])]
'''
