import argparse

def argparser(*args: str):
    ''' 
    :param *args: ordered list of strings taken from cmd line
    usage: access property of the return of this func u
    args = argparser("path")
    path = args.path
    '''
    parser = argparse.ArgumentParser()
    for arg in args:
        parser.add_argument(arg)
    return parser.parse_args()