import pysam
import sys

def main(arguments):
    pstrand_in = pysam.AlignmentFile(arguments[1], 'rb')
    mstrand_in = pysam.AlignmentFile(arguments[2], 'rb')
    nostrand_in = pysam.AlignmentFile(arguments[3], 'rb')

    pstrand_out = pysam.AlignmentFile(arguments[4], 'wb', template=pstrand_in)
    mstrand_out = pysam.AlignmentFile(arguments[5], 'wb', template=mstrand_in)
    nostrand_out = pysam.AlignmentFile(arguments[6], 'wb', template=nostrand_in)

    for read in pstrand_in.fetch(until_eof=True):
        if read.has_tag('XS'):
            XS = read.get_tag('XS')
            read.set_tag('ES', XS)
        else:
            read.set_tag('ES', 'Unassigned')
        if read.has_tag('XT'):
            XT = read.get_tag('XT')
            read.set_tag('GE', XT)
        pstrand_out.write(read)
    for read in mstrand_in.fetch(until_eof=True):
        if read.has_tag('XS'):
            XS = read.get_tag('XS')
            read.set_tag('ES', XS)
        else:
            read.set_tag('ES', 'Unassigned')
        if read.has_tag('XT'):
            XT = read.get_tag('XT')
            read.set_tag('GE', XT)
        mstrand_out.write(read)
    for read in nostrand_in.fetch(until_eof=True):
        if read.has_tag('XS'):
            XS = read.get_tag('XS')
            read.set_tag('ES', XS)
        else:
            read.set_tag('ES', 'Unassigned')
        if read.has_tag('XT'):
            XT = read.get_tag('XT')
            read.set_tag('GE', XT)
        nostrand_out.write(read)
    
    pstrand_in.close()
    pstrand_out.close()
    mstrand_in.close()
    mstrand_out.close()
    nostrand_in.close()
    nostrand_out.close()
    
if __name__ == "__main__":
    main(sys.argv)


