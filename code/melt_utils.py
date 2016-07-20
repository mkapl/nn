import sys
import math
import random
import string
import subprocess as sp
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

pfunc_path = '/Users/Matias/nn/library/nupack3.0.6/bin/pfunc'


alphabet = string.ascii_lowercase + string.ascii_uppercase
def rand_string(n=6):
    """
    random string of lowercase characters of length n
    """
    return ''.join(random.choice(alphabet) for _ in range(n))


def get_pf(seq, T=37, constraint=False):
    """
    calculate partition function at a specific temperature
    optionally with a constraint
    """
    # write sequence and constraint to nupack input file
    name = rand_string()
    with open('%s.in' % name, 'w') as f:
        f.write('%s\n' % seq)
        if constraint:
            f.write('%s\n' % constraint)

    # get command as list
    # commented out 7/20/16 because I would get "pfunc: unrecognized option `-onstraint'"
    if constraint:
        command = [pfunc_path, '-material', 'dna', '-constraint', '-T', str(T), name]
    else:
        command = [pfunc_path, '-material', 'dna', '-T', str(T), name]

    # call pfunc and and get partition function value
    result = sp.check_output(command, stderr=sp.STDOUT)
    sp.call(['rm', '%s.in' % name])
    return float(result.split('\n')[-2])


def get_prop_closing_pair(seq, T=37, nopen=1):
    """
    calculate the proportion of RNAs with formation of closing pair
    """
    n = len(seq)

    # calculate partition function for forced closing pair
    constraint = 'x' * nopen + '.' * (n-2*nopen) + 'x' * nopen
    pf_unpaired = get_pf(seq, T, constraint)

    # calculate overall partition function
    pf_full = get_pf(seq, T)

    return pf_unpaired/pf_full


def get_melting_curve(seq, Trange=False):
    """
    get simulated melting curve with NUPACK
    """
    # get temperature range
    if not Trange:
        Trange = range(0, 101, 5)

    # get melting curve
    curve = []
    for T in Trange:
        curve.append(get_prop_closing_pair(seq, T))
    return curve


def plot_melting_curve(seq, out, Trange=False):
    """
    generate plot of melting curve
    """
    if not Trange:
        Trange = range(0, 101, 5)
    curve = get_melting_curve(seq, Trange)

    fig = plt.figure()
    plt.plot(Trange, curve)
    plt.axvline(x=20, color='k', linestyle=':')
    plt.axvline(x=55, color='k', linestyle=':')
    fig.savefig(out)


def estimate_tm(Trange, curve):
    """
    estimate melting temperature based on coarse curve
    """
    if len(Trange) != len(curve):
        raise ValueError('melting curve x and y have different lengths')
    ix = next(i for i,x in enumerate(curve) if x > 0.5)
    diff = (0.5 - curve[ix-1]) / (curve[ix] - curve[ix-1])
    return Trange[ix-1] + (Trange[ix] - Trange[ix-1]) * diff


def plot_melting_curves(seqs, out, Trange=False, overlap=True):
    """
    generate plot of melting curve
    """
    if not Trange:
        Trange = range(0, 101, 5)

    n = len(seqs)
    ncol = int(math.ceil(float(n)/4))
    if overlap:
        fig = plt.figure(figsize=(5, 4))
    else:
        fig, ax = plt.subplots(4, ncol, figsize=(ncol*5, 12))
    colors = ['#0066cc', '#99cc00', '#009900', '#ff9999', '#ff0000', '#ffdb4d',
              '#ff9900', '#d9b3ff']
    colors = colors * int(math.ceil(float(n)/8))
    for n,seq in enumerate(seqs):
        if overlap:
            obj = plt
            plt.xlabel('temperature')
        else:
            obj = ax[n % 4, n / 4]
            obj.set_title(seq.strip())
            obj.set_xlabel('temperature')
        curve = get_melting_curve(seq.strip(), Trange)
        obj.plot(Trange, curve, color=colors[n])
        obj.axvline(x=20, color='k', linestyle=':')
        obj.axvline(x=55, color='k', linestyle=':')
        if curve[0] < 0.5 and curve[-1] > 0.5 and not overlap:
            obj.axvline(x=estimate_tm(Trange, curve), ymax=0.5, color='k',
                        linestyle='--')
    if not overlap:
        fig.subplots_adjust(hspace=0.5)
        fig.text(0.05, 0.5, 'proportion fluorescent (unpaired closing pair)',
                 ha='center', va='center', rotation='vertical')
    else:
        plt.ylabel('proportion fluorescent\n(unpaired closing pair)')
        plt.tight_layout()
    fig.savefig(out)


def main():
    with open(sys.argv[1]) as f:
        seqs = f.readlines()
    plot_melting_curves(seqs, sys.argv[2], overlap=False)


if __name__ == '__main__':
    main()
