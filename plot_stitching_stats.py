import matplotlib
matplotlib.use('Agg')
from   matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot  as plt
import matplotlib.patches as patchs
import matplotlib.image   as mpimg

import gzip
import numpy
import pandas
import seaborn

import sys

def main():
    fastq  = sys.argv[1]
    log    = sys.argv[2]
    output = sys.argv[3]

    data           = gzip.open(fastq)
    count          = 0
    num_overlaps   = []
    num_mismatches = []
    for line in data:
        if count % 4 == 0:
            tokens = line.decode().strip().split("_")
            num_overlaps.append(int(tokens[1]))
            num_mismatches.append(int(tokens[2]))
        count += 1
    data.close()

    num_overlaps   = numpy.array(num_overlaps)
    num_mismatches = numpy.array(num_mismatches)
    df = pandas.DataFrame({'Number of Overlapping Bases':num_overlaps, 'Number of Mismatched Bases':num_mismatches})
    pp = PdfPages(output)
    print("Plotting joint distribution")
    fig = seaborn.jointplot(x="Number of Overlapping Bases", y="Number of Mismatched Bases", data=df, kind="kde")
    pp.savefig(fig.fig)
    pp.close()



if __name__ == "__main__":
    main()
