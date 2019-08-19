import sys, optparse

"""
groups together lines of the same motif with slightly different footprint coordinates
filter3.py -i input.bed -o output.bed
"""

def parse(l):
    parts = l.strip().split('\t')
    parts[1] = int(parts[1])
    parts[2] = int(parts[2])
    return parts


def fmt(p):
    return '\t'.join([str(i) for i in p])

def __main__():
    parser=optparse.OptionParser()
    parser.add_option( '-i', '--input', dest='inputF', action='store', type="string", help='Input bed file' )
    parser.add_option( '-o', '--output', dest='outputF', action='store', type="string", default=None, help='Output bed file' )
    (options, args) = parser.parse_args()

    f=open(options.outputF, 'w')
    with open(options.inputF) as fp:
        last_parts = parse(next(fp))

        for line in fp:
            parts = line.strip().split('\t')
            parts[1] = int(parts[1])
            parts[2] = int(parts[2])

            # same chromosome, same motif, same motif_name
            if parts[0] == last_parts[0] and parts[-1] == last_parts[-1] and parts[3] == last_parts[3]:
                if parts[1] < last_parts[2]:
                    # within the interval
                    last_parts[2] = parts[2]
                    continue
 
            f.write("%s\n" % fmt(last_parts))
            last_parts = parts
    f.close()

if __name__=="__main__":
        __main__()
