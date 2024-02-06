from multiprocessing import Pool
from itertools import islice


def batched(iterable, n):
    ''' Given an iterable creates a generator that return batches of n elements of iterable'''
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError('n must be at least one')
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


def flatten_list_of_lists(l):
    '''l is a list of lists that have to be flattened to get a list of els'''
    return [el for list_ in l for el in list_]


def apply_f_parallelized_batched(f, iterable_, nels_per_batch, workers=1):
    ''' 
    given an iterable_, applies f to a batch of iterable_ of size nels_per_batch
    the batches are treated in parallel over workers
    f must be appliable over a list of elements
    '''
    it = batched(iterable_, nels_per_batch)
    with Pool(processes=workers) as P:
        out = P.map(f, it)
    return flatten_list_of_lists(out)
