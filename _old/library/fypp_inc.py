def ranksuffix(rank):
    '''Returns a Fortran rank specification suffix.

    Args:
        rank (int): Rank of the object (must be >= 0).

    Returns:
        str: Rank specification suffix (e.g. (:)) or empty string for rank = 0.
    '''
    if rank == 0:
        return ''
    else:
        return '(' + ','.join([':'] * rank) + ')'

def ranksuffix_last(rank,last):
    '''Same as ranksuffix but replaces last : with an i
    '''
    if rank == 0:
        return ''
    elif rank == 1:
        return '(' + last + ')'
    else:
        return '(' + ','.join([':'] * (rank-1)) + ',' + last + ')'

def variants(name, prefixes=None, suffixes=None, separator=', '):
    '''Returns all possible variants of a name.

    Args:
        name (str): Name to return the variants of.
        prefixes (list of str): Prefixes to use for building variants.
        suffixes (list of str): Suffixes to use for building variants.
        separator (str): Separator to use between variants.

    Returns:
        str: All combinations of the form <prefix><name><suffix> separated
            by the separator.
    '''
    if prefixes is None:
        prefixes = ['']
    if suffixes is None:
        suffixes = ['']
    variants = [prefix + name + suffix
                for prefix in prefixes for suffix in suffixes]
    return separator.join(variants)




