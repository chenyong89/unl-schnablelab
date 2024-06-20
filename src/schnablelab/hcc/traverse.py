# -*- coding: UTF-8 -*-

"""
traverse files to avoid purge policy on hcc
"""

import os
import sys
import random
from schnablelab.apps.base import ActionDispatcher, OptionParser


def main():
    actions = (
        ('walk', 'traverse files using python os.walk'),
        ('find', 'traverse fiels using find command'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def walk(args):
    """
    %prog walk dir

    traverse files using python os.walk
    """
    p = OptionParser(walk.__doc__)
    p.add_option("--frac", default=0.1, type=float,
                 help="percentage of files that will be read")
    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())

    folder, = args
    all_fns = []
    for dirpath, dirnames, filenames in os.walk(folder):
        for filename in filenames:
            fn = os.path.join(dirpath, filename)
            all_fns.append(fn)
    part_fns = random.sample(all_fns, int(len(all_fns)*opts.frac))
    for i in part_fns:
        print(i)
        f = open(fn)
        f.readline()
        f.close()


def find(args):
    """
    %prog find dir

    traverse fiels using find command
    """
    p = OptionParser(find.__doc__)
    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())
    folder, = args
    cmd = 'find %s -exec touch {} +' % folder
    print(cmd)


if __name__ == "__main__":
    main()
