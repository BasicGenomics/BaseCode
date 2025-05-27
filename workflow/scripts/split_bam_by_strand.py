import pysam
import sys

def split_bam_by_strand(inbam):
    ### this function parses reads mapped by HISAT2-3N, splits the output bam file by inferred strand while integrating paired information if possible
    ### is also used to parse out information from extract&tag step stored in the readID
    outbams = [
        inbam.replace(".trimmed.aligned.bam",".trimmed.aligned.nostrand.bam"),
        inbam.replace(".trimmed.aligned.bam",".trimmed.aligned.pstrand.bam"),
        inbam.replace(".trimmed.aligned.bam",".trimmed.aligned.mstrand.bam")
    ]
    save = pysam.set_verbosity(0)
    inbamfh = pysam.AlignmentFile(inbam, "rb")
    #inbamfh = pysam.AlignmentFile(inbam, "rb", ignore_truncation = True)
    x = pysam.set_verbosity(save)
    outbamfhs = [pysam.AlignmentFile(f, "wb", template = inbamfh) for f in outbams]
    for read in inbamfh:
        try:
            strand = read.get_tag('YZ') #try to fetch tag and except it, possible values are + and - or tag not set
        except:
            strand = 'no'
            read.set_tag(tag = "YZ", value_type = "Z", value = "NA")
        readidlist = read.query_name.split("-")
        
        if len(readidlist) == 9:
            readidlist[7] = readidlist[7]+'-'+readidlist[8]
        read.query_name = readidlist[0]
        if not read.has_tag('Yf'):
            read.set_tag(tag = 'Yf', value_type = "i", value = 0)

        if not read.is_paired:
            read.set_tag(tag = "BC", value_type = "Z", value = readidlist[1].replace("BC:",""))
            read.set_tag(tag = "XX", value_type = "Z", value = readidlist[2].replace("XX:",""))
            read.set_tag(tag = "UB", value_type = "Z", value = readidlist[3].replace("UB:",""))
            read.set_tag(tag = "DT", value_type = "Z", value = readidlist[4].replace("DT:",""))
            read.set_tag(tag = "SB", value_type = "Z", value = readidlist[5].replace("SB:",""))
            read.set_tag(tag = "SM", value_type = "Z", value = readidlist[6].replace("SM:",""))
            read.set_tag(tag = "RS", value_type = "i", value = int(readidlist[7].replace("RS:","")))

        if read.is_paired and read.is_read1:
            
            read2 = next(inbamfh)
            try:
                strand2 = read2.get_tag('YZ') #try to fetch tag and except it, possible values are + and - or tag not set
            except:
                strand2 = 'no'
                read2.set_tag(tag = "YZ", value_type = "Z", value = "NA")
            
            readidlist2 = read2.query_name.split("-")
            if len(readidlist2) == 9:
                readidlist2[7] = readidlist2[7]+'-'+readidlist2[8]
            read2.query_name = readidlist[0]
            if strand != strand2 and read.query_name == read2.query_name: # if the reads are paired with each other but have conflicting strand, we can combine the information
                if (strand == "no" and strand2 == "+") or (strand2 == "no" and strand == "+"):
                    strand = "+"
                    strand2 = "+"
                elif (strand == "no" and strand2 == "-") or (strand2 == "no" and strand == "-"):
                    strand = "-"
                    strand2 = "-"
                else:
                    strand = "no"
                    strand2 = "no"
                read.set_tag(tag = "YZ", value_type = "Z", value = strand)
                read2.set_tag(tag = "YZ", value_type = "Z", value = strand2)
            if not read.has_tag('Yf'):
                read.set_tag(tag = 'Yf', value_type = "i", value = 0)
            if not read2.has_tag('Yf'):
                read2.set_tag(tag = 'Yf', value_type = "i", value = 0)

            read.set_tag(tag = "BC", value_type = "Z", value = readidlist[1].replace("BC:",""))
            read.set_tag(tag = "XX", value_type = "Z", value = readidlist[2].replace("XX:",""))
            read.set_tag(tag = "UB", value_type = "Z", value = readidlist[3].replace("UB:",""))
            read.set_tag(tag = "DT", value_type = "Z", value = readidlist[4].replace("DT:",""))
            read.set_tag(tag = "SB", value_type = "Z", value = readidlist[5].replace("SB:",""))
            read.set_tag(tag = "SM", value_type = "Z", value = readidlist[6].replace("SM:",""))
            read.set_tag(tag = "RS", value_type = "i", value = int(readidlist[7].replace("RS:","")))

            read2.set_tag(tag = "BC", value_type = "Z", value = readidlist2[1].replace("BC:",""))
            read2.set_tag(tag = "XX", value_type = "Z", value = readidlist2[2].replace("XX:",""))
            read2.set_tag(tag = "UB", value_type = "Z", value = readidlist2[3].replace("UB:",""))
            read2.set_tag(tag = "DT", value_type = "Z", value = readidlist2[4].replace("DT:",""))
            read2.set_tag(tag = "SB", value_type = "Z", value = readidlist2[5].replace("SB:",""))
            read2.set_tag(tag = "SM", value_type = "Z", value = readidlist2[6].replace("SM:",""))
            read2.set_tag(tag = "RS", value_type = "i", value = int(readidlist2[7].replace("RS:","")))
            
            if strand == "+":
                outbamfhs[1].write(read)
            elif strand == "-":
                outbamfhs[2].write(read)
            else:
                outbamfhs[0].write(read)
            if strand2 == "+":
                outbamfhs[1].write(read2)
            elif strand2 == "-":
                outbamfhs[2].write(read2)
            else:
                outbamfhs[0].write(read2)

        else: #singleton reads are done, just write out
            if strand == "+":
                outbamfhs[1].write(read)
            elif strand == "-":
                outbamfhs[2].write(read)
            else:
                outbamfhs[0].write(read)
    inbamfh.close()
    [b.close() for b in outbamfhs]
    return 1


if __name__ == "__main__":
    split_bam_by_strand(sys.argv[1], set([]))