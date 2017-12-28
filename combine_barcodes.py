import sys
import os.path
import itertools

work_dir = sys.argv[1] 

class FastqRead(object):
    def __init__(self, read):
        self.desc, self.seq, self.qual = read
    
    def __repr__(self):
        return self.desc + "\n" + self.seq + "\n+\n" + self.qual + "\n"

def _grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3) --> ABC DEF
    args = [iter(iterable)] * n
    return itertools.izip(*args)

def parse_fastq(f):
    for desc, seq, _, qual in _grouper(f, 4):
        desc = desc.rstrip()[1:]
        seq = seq.rstrip()
        qual = qual.rstrip()
        yield desc, seq, qual

I1 = os.path.join(work_dir, "Undetermined_S0_L001_I1_001.fastq")
I2 = os.path.join(work_dir, "Undetermined_S0_L001_I2_001.fastq")
I = os.path.join(work_dir, "Undetermined_S0_L001_I12_001.fastq")

I1_handle = open(I1)
I2_handle = open(I2)    
fwds = (FastqRead(x) for x in parse_fastq(I1_handle))
revs = (FastqRead(x) for x in parse_fastq(I2_handle))

with open(I, "w") as f_out:
    for fwd, rev in itertools.izip(fwds, revs):
        f_out.write("@%s\n%s\n+\n%s\n" % (fwd.desc, fwd.seq+rev.seq, fwd.qual+rev.qual))

I1_handle.close()
I2_handle.close()
