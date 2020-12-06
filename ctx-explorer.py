#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
@author: Lukasz M. Szafron
@email: lukszafron@gmail.com
"""
appname = "CTX-Explorer"

import subprocess, re, os, gzip, sys, getopt, statistics
from itertools import groupby, chain
from termcolor import colored

def usage():
    print(
        "\nWelcome to the",colored(appname,"green")+" app.\n\n"
        "The following options are available:\n\n"
        "\t-h, --help:\t\tprints this help message.\n"
        "\t-t, --threads:\t\tdefines the number of CPU threads that are to be used (default: 1).\n"
        "\t-T, --tmpdir:\t\tdefines the directory where temporary files are to be stored (default: current dir).\n"
        "\t-q, --qual:\t\tdefines the minimal Phred quality score of mappings to be used (default: 20).\n"
        "\t-b, --bamfile:\t\tthe name of a bamfile to be used",colored("(MANDATORY)", "green")+".\n"
        "\t-I, --intractx:\t\tindicates if intra-chromosomal translocations should be evaluated (default: 'no').",colored("This option significantly increases memory usage.","green")+"\n"
        "\t-1, --chrom1:\t\tthe name of the first chromosome (optional) involved in a translocation (e.g., 'chr1').\n"
        "\t-2, --chrom2:\t\tthe name of the second chromosome (optional) involved in a translocation (e.g., 'chr2').\n"
        "\t-l, --tlen:\t\tminimal distance (in nucleotides) between the first and the second read forming a read pair (default: 1,000,000).\n"
        "\t-i, --insert:\t\tthe maximum insert size used for identification of CTX-supporting reads and hit pairs (default: median insert size * 2).\n"
        "\t-n, --nohits:\t\tminimal number of hits for a translocation to be stored in the final report (default: 2).\n"
        "\t-N, --nohits_sec:\tminimal number of hits per the second chromosome for a translocation to be stored in the final report (default: 2).\n"
        "\t-s, --min_size:\t\tminimal size of a hit group for a translocation to be stored in the final report (default: 2).\n"
        "\t-d, --no_filter:\tspecifies whether the filtering by the number of supporting reads should be turned off (default: 'no').\n"
        "\t-p, --prefix:\t\tspecifies a prefix of a file in which the final report will be saved (default: 'output').\n"
        "\t-g, --gzipped:\t\tindicates whether the final report file should be gzipped (default: 'no').\n"
        )
# The next two lines are for debugging purposes only and should be commented in the final program.
#option_list = ["-b", "/workspace/lukasz/NGS-all-in-one/RUNS/BGRYG/MAPPINGS_TRIMMED/Homo_sapiens.GRCh38.p12.genome/STAR/BGRYG-93.sorted.bam", "-t","6","-T","/tmpdir/lukasz","--prefix","outfile3","-1","chr8", "-2", "chr14"]
#opts, args = getopt.getopt(option_list, "ht:T:q:b:I1:2:l:i:n:N:s:dp:g", ["help","threads=","tmpdir=","qual=","bamfile=","intractx","chrom1=","chrom2=","tlen=","insert=","nohits=","nohits_sec=","min_size=","no_filter","prefix=","gzipped"])

try:
        opts, args = getopt.getopt(sys.argv[1:], "ht:T:q:b:I1:2:l:i:n:N:s:dp:g", ["help","threads=","tmpdir=","qual=","bamfile=","intractx","chrom1=","chrom2=","tlen=","insert=","nohits=","nohits_sec=","min_size=","no_filter","prefix=","gzipped"])
        if len(opts) == 0:
                usage()
                sys.exit()
        for o, a in opts:
            if o in ("-h", "--help"):
                usage()
                sys.exit()
            elif o in ("-t", "--threads"):
                threads = a
            elif o in ("-T", "--tmpdir"):
                tmpdir = a
            elif o in ("-q", "--qual"):
                qual = a
            elif o in ("-b", "--bamfile"):
                bamfile = a
            elif o in ("-I", "--intractx"):
                intractx = True
            elif o in ("-1", "--chrom1"):
                chr1 = a
            elif o in ("-2", "--chrom2"):
                chr2 = a
            elif o in ("-l", "--tlen"):
                tlen = a
            elif o in ("-i", "--insert"):
                insert = a
            elif o in ("-n", "--nohits"):
                nohits = a
            elif o in ("-N", "--nohits_sec"):
                nohits_secchr = a
            elif o in ("-s","--min_size"):
                group_size = a
            elif o in ("-d","--no_filter"):
                no_filter = True
            elif o in ("-p", "--prefix"):
                prefix = a
            elif o in ("-g", "--gzipped"):
                compressed_output = True
            else:
                assert False, "Unhandled option: "+o

except getopt.GetoptError as err:
    # print help information and exit:
    print("\n"+colored(str(err),"red")) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

try:
    threads
except:
    threads='1'
try:
    tmpdir
except:
    tmpdir = os.getcwd()
try:
    qual
except:
    qual = 20
try:
    bamfile
except:
    usage()
    raise(Exception("The bam file was not provided. The program will exit now."))
try:
    intractx
except:
    intractx = False
try:
    tlen
except:
    tlen = 1e6
print("Calculating basic BAM file statistics...")
try:
    insert
    all_reads = subprocess.check_output(["samtools","view","-@", threads,"-c", bamfile]).decode().strip()
    print("Maximum insert size defined by the user is:",colored("{}","green").format(insert),"bp.")
except:
     bamfile_stats = subprocess.check_output(["bamtools","stats","-insert","-in", bamfile])
     bamfile_stats = bamfile_stats.strip().decode().splitlines()
     all_reads = bamfile_stats[4].split(":")[1].strip()
     insert_median = round(float(bamfile_stats[-1].split(":")[1].strip()))
     insert = 2*insert_median
     print("Median insert size equals:", colored("{}","green").format(insert_median),"bp, so the maximum insert size is set to:", colored("{}","green").format(insert),"bp.")
try:
    nohits
except:
    nohits = 2
try:
    nohits_secchr
except:
    nohits_secchr = 2
try:
    group_size
except:
    group_size = 2
try:
    no_filter
except:
    no_filter = False
try:
    prefix
except:
    prefix = "output"
try:
    compressed_output
except:
    compressed_output = False

if intractx:
    filterflag = "-F12"
    print(colored("Intra-chromosomal translocations will be included in the analysis.","green"))
else:
    filterflag = "-F14"
    print(colored("Intra-chromosomal translocations will be excluded from the analysis.","green"))

bamfile_header = subprocess.check_output(["samtools","view","-@",threads,"-H",bamfile]).decode()

print("Getting coordinates for all mapped reads...")

mapped_lines = subprocess.check_output(["samtools","view","-@",threads, filterflag, bamfile]).decode().splitlines()
len_mapped_lines = len(mapped_lines)

def sel_element(lst, *indices):
    return (lst[i] for i in indices)

mapped_lines_coords = ['\t'.join(list(sel_element(value.split("\t"),1,2,3,6,5,8,0))) for value in mapped_lines] # keep all bitwise flags, chromosome names, mapping position coordinates, CIGAR strings, TLENs and read names.
del(mapped_lines)

def sort_chroms(value):
    return (value.split("\t")[1],value.split("\t")[3])
mapped_lines_coords.sort(key=sort_chroms) # sort coordinates by chromosomes

def group_chroms(value):
    return value.split("\t")[1]+"\t"+value.split("\t")[3]
mapped_lines_coords = [(name,list(group)) for name,group in groupby(mapped_lines_coords,key=group_chroms)] # group coordinates by chrom. names.

print("Filtering out unmapped and low-quality reads...")

filtered_lines = subprocess.check_output(["samtools","view","-@",threads, filterflag,"-q",str(qual),bamfile]).decode().splitlines()
len_filtered_lines = len(filtered_lines)

print("CIGAR-based filtering of the reads...")

filtered_lines = [value for value in filtered_lines if re.search(pattern="\d+(S|H)", string=value.split("\t")[5])] # keep only mappings with clippings in the CIGAR strings.
filtered_lines = [value for value in filtered_lines if not re.search(pattern="^\d+(S|H).*\d+(S|H)$", string=value.split("\t")[5])] # keep only mappings with one overhang.

len_cigar_filtered = len(filtered_lines)

print("Defining the regions of potential chromosomal translocations...")

filtered_lines = ['\t'.join(flowvalue.split(sep="\t")[0:10]) for flowvalue in filtered_lines] # keep only the first 10 mandatory fields of the bam file (1-10).

for index,flowvalue in enumerate(filtered_lines): # Separates matchings from overhangs.
    if re.search(pattern="^\d+S.*$", string=flowvalue.split(sep="\t")[5]):
        shift = int(flowvalue.split("\t")[5].split("S")[0])
        filtered_lines[index] = flowvalue+"\t"+flowvalue.split("\t")[9][shift:]+"\t"+flowvalue.split("\t")[9][:shift]+"\t"+str(shift)
    elif re.search(pattern="^.*\d+S$", string=flowvalue.split(sep="\t")[5]):
        shift = int(re.findall(pattern="\d+",string=flowvalue.split(sep="\t")[5].split("S")[0])[-1])
        filtered_lines[index] = flowvalue+"\t"+flowvalue.split("\t")[9][:-shift]+"\t"+flowvalue.split("\t")[9][-shift:]+"\t"+str(shift)
    elif re.search(pattern="\d+H", string=flowvalue.split(sep="\t")[5]):
        filtered_lines[index] = flowvalue+"\t"+flowvalue.split("\t")[9]+"\t\t0"
    else: raise Exception("CIGAR string analysis failed.")

filtered_lines.sort(key=lambda x:x.split("\t")[2]) # sort reads by chr1

def group_chr(value):
    return value.split(sep="\t")[2]

filtered_lines = [[name,list(group)] for name,group in groupby(filtered_lines,key=group_chr)] # group by chr1

locations = [[] for chrom in filtered_lines]
for i1,v1 in enumerate(filtered_lines):
    chrom = v1[0]
    for i_flowvalue,flowvalue in enumerate(v1[1]):
        if flowvalue.split(sep="\t")[2] == chrom:
            locations[i1].append([v for v in sel_element(flowvalue.split(sep="\t"),3,5)])
        else:
            raise Exception("Values assignment error occurred.")
        if re.search(pattern="^\d+(S|H).*$", string=locations[i1][i_flowvalue][1]):
            locations[i1][i_flowvalue].append(locations[i1][i_flowvalue][0])
        elif re.search(pattern="^.*\d+(S|H)$", string=locations[i1][i_flowvalue][1]):
            shift = sum(shift for shift in list(map(int,re.findall(string=''.join(re.findall(pattern="(\d+M|\d+D|\d+N|\d+P)", string=locations[i1][i_flowvalue][1])),pattern="\d+")))) # Defines the exact ctx position.
            locations[i1][i_flowvalue].append(int(locations[i1][i_flowvalue][0])+shift)
        else:
            raise Exception("Determination of chromosomal translocation position failed.")

for i1,v1 in enumerate(filtered_lines):
    for i_flowvalue,flowvalue in enumerate(v1[1]):
        if flowvalue.split(sep="\t")[3] == locations[i1][i_flowvalue][0] and flowvalue.split(sep="\t")[5] == locations[i1][i_flowvalue][1]:
            filtered_lines[i1][1][i_flowvalue] = flowvalue+"\t"+str(locations[i1][i_flowvalue][2]) # add CTX location at the end of the read.
        else:
            raise Exception("Values assignment error occurred.")
del(locations)

def loc_sort(value):
    return int(value.split(sep="\t")[13]) # sort by ctx
for chrom in filtered_lines:
    chrom[1].sort(key=loc_sort)
try:
    del(chrom)
except:
    print(colored("The initial filtering step returned no reads.", "red"))

print("Calculating the number of hits for each translocation...")

for i1,v1 in enumerate(filtered_lines):
    filtered_lines[i1][1] = [[name,list(group)] for name,group in groupby(v1[1],key=lambda x:x.split(sep="\t")[13])] # group by ctx

filtered_lines = [v2 for v1 in filtered_lines for v2 in v1[1] if len(v2[1]) >= int(nohits)] # filter by nohits

def chrom_regroup(lst):
    for value in lst[1]:
        return value.split("\t")[2]
filtered_lines = [[name,list(group)] for name,group in groupby(filtered_lines,key=chrom_regroup)] # regroup data by chromosomes

for i1,v1 in enumerate(filtered_lines):
    for i2,v2 in enumerate(v1[1]):
        count = len(v2[1])
        for index,value in enumerate(v2[1]):
            filtered_lines[i1][1][i2][1][index] = "\t".join([value,str(count)]) # add total number of hits for each ctx.

counter = 0
for v1 in filtered_lines:
    for v2 in v1[1]:
        counter = counter + len(v2[1])
len_nohits = counter

for i1,v1 in enumerate(filtered_lines):
    for i2,v2 in enumerate(v1[1]):
        filtered_lines[i1][1][i2][1].sort(key=lambda x:x.split("\t")[6]) # sort by chr2

for i1,v1 in enumerate(filtered_lines):
    for i2,v2 in enumerate(v1[1]):
        filtered_lines[i1][1][i2][1] = [[name,list(group)] for name,group in groupby(filtered_lines[i1][1][i2][1],key=lambda x:x.split("\t")[6])] # group by chr2

final_hits = []
for i1,v1 in enumerate(filtered_lines):
    final_hits.append([v3 for v2 in v1[1] for v3 in v2[1] if len(v3[1]) >= int(nohits_secchr)]) # filter reads by secchr
del(filtered_lines)
final_hits = list(chain.from_iterable(final_hits)) # flatten the list
final_hits = [[name,list(group)] for name,group in groupby(final_hits,key=chrom_regroup)] # regroup data by chrom1.

def ctx_regroup(lst):
    for value in lst[1]:
        return value.split("\t")[13]
for i1,v1 in enumerate(final_hits):
    final_hits[i1][1] = [[name,list(group)] for name,group in groupby(v1[1],key=ctx_regroup)] # regroup by ctx position

for i1,v1 in enumerate(final_hits):
    for i2,v2 in enumerate(v1[1]):
        for i3,v3 in enumerate(v2[1]):
            count = str(len(v3[1]))
            for index,value in enumerate(v3[1]):
                final_hits[i1][1][i2][1][i3][1][index] = "\t".join([value,count]) # add number of hits for each ctx per secchr.

counter = 0
for v1 in final_hits:
    for v2 in v1[1]:
        for v3 in v2[1]:
            counter = counter + len(v3[1])
len_nohits_secchr = counter

final_hits = [value for v1 in final_hits for v2 in v1[1] for v3 in v2[1] for value in v3[1] if (int(value.split("\t")[8]) >= int(tlen) or int(value.split("\t")[8]) <= -int(tlen)) or (int(value.split("\t")[8]) == 0 and value.split("\t")[6] != "=")] # filter reads by tlen and flatten the list.
len_tlen = len(final_hits)

final_hits.sort(key=lambda x:(x.split("\t")[2],x.split("\t")[13])) # sort hits by chr1 and ctx position.
final_hits = [(name,list(group)) for name,group in groupby(final_hits,key=lambda x:(x.split("\t")[2],x.split("\t")[13]))] # group hits by chr1 and ctx position.
final_hits = ["\t".join([v2,str(len(v1[1]))]) for v1 in final_hits for v2 in v1[1] if len(v1[1]) >= int(group_size)] # add sizes of the hit groups after tlen filtering to each hit, filter hits by the group size and then flatten the list.
len_group_size = len(final_hits)

final_hits = [[name,list(group)] for name,group in groupby(final_hits,key=lambda x:x.split("\t")[2])] # regroup hits by chrom1.
for i1,v1 in enumerate(final_hits):
    final_hits[i1][1] = [[name,list(group)] for name,group in groupby(v1[1],key=lambda x:x.split("\t")[13])] # regroup hits by ctx.
for i1,v1 in enumerate(final_hits):
    for i2,v2 in enumerate(v1[1]):
        final_hits[i1][1][i2][1] = [[name,list(group)] for name,group in groupby(v2[1],key=lambda x:x.split("\t")[6])] # regroup hits by chrom2.

print("Defining the number of reads supporting each translocation...")

forward_flags = [65, 67, 129, 131, 97, 99, 161, 163, 321, 323, 385, 387, 353, 355, 417, 419, 1089, 1091, 1153, 1155, 1121, 1123, 1185, 1187, 1345, 1347, 1409, 1411, 1377, 1379, 1441, 1443, 2113, 2115, 2177, 2179, 2145, 2147, 2209, 2211, 2369, 2371, 2433, 2435, 2401, 2403, 2465, 2467, 3137, 3139, 3201, 3203, 3169, 3171, 3233, 3235, 3393, 3395, 3457, 3459, 3425, 3427, 3489, 3491]
reverse_flags = [81, 83, 145, 147, 113, 115, 177, 179, 337, 339, 401, 403, 369, 371, 433, 435, 1105, 1107, 1169, 1171, 1137, 1139, 1201, 1203, 1361, 1363, 1425, 1427, 1393, 1395, 1457, 1459, 2129, 2131, 2193, 2195, 2161, 2163, 2225, 2227, 2385, 2387, 2449, 2451, 2417, 2419, 2481, 2483, 3153, 3155, 3217, 3219, 3185, 3187, 3249, 3251, 3409, 3411, 3473, 3475, 3441, 3443, 3505, 3507]
primary_alignments = [65, 67, 129, 131, 97, 99, 161, 163, 81, 83, 145, 147, 113, 115, 177, 179, 1089, 1091, 1153, 1155, 1121, 1123, 1185, 1187, 1105, 1107, 1169, 1171, 1137, 1139, 1201, 1203, 2113, 2115, 2177, 2179, 2145, 2147, 2209, 2211, 2129, 2131, 2193, 2195, 2161, 2163, 2225, 2227, 3137, 3139, 3201, 3203, 3169, 3171, 3233, 3235, 3153, 3155, 3217, 3219, 3185, 3187, 3249, 3251]
secondary_alignments = [321, 323, 385, 387, 353, 355, 417, 419, 337, 339, 401, 403, 369, 371, 433, 435, 1345, 1347, 1409, 1411, 1377, 1379, 1441, 1443, 1361, 1363, 1425, 1427, 1393, 1395, 1457, 1459, 2369, 2371, 2433, 2435, 2401, 2403, 2465, 2467, 2385, 2387, 2449, 2451, 2417, 2419, 2481, 2483, 3393, 3395, 3457, 3459, 3425, 3427, 3489, 3491, 3409, 3411, 3473, 3475, 3441, 3443, 3505, 3507]

def supporting_reads(lst,ctx):
    counter = 0
    for value in lst[1]:
        if lst[0].split("\t")[1] != "=": # ctx on different chromosomes
            if int(value.split("\t")[0]) in forward_flags: # forward bitwise flags
                if not re.search(pattern="^\d+(S|H).*$",string=value.split("\t")[4]):
                    shift = sum(shift for shift in list(map(int,re.findall(string=''.join(re.findall(pattern="(\d+M|\d+D|\d+N|\d+P)",string=value.split("\t")[4])),pattern="\d+")))) # Defines the exact end position of the read +1.
                    if int(value.split("\t")[2])+shift <= int(ctx) and int(value.split("\t")[2])+shift >= int(ctx) - int(insert):
                        counter += 1
                else:
                    if int(value.split("\t")[2]) == int(ctx):
                        counter += 1
            elif int(value.split("\t")[0]) in reverse_flags: # reverse bitwise flags
                if not re.search(pattern="^.*\d+(S|H)$",string=value.split("\t")[4]):
                    if int(value.split("\t")[2]) >= int(ctx) and int(value.split("\t")[2]) <= int(ctx) + int(insert):
                        counter += 1
                else:
                    shift = sum(shift for shift in list(map(int,re.findall(string=''.join(re.findall(pattern="(\d+M|\d+D|\d+N|\d+P)",string=value.split("\t")[4])),pattern="\d+"))))
                    if int(value.split("\t")[2])+shift == int(ctx):
                        counter += 1
            else:
                raise Exception(''.join (["Bitwise flags evaluation failed. The following flag is missing: ", value.split("\t")[0],". Please, check if all the reads are paired."]))
        else: # ctx on the same chromosome
            if int(value.split("\t")[0]) in forward_flags and (int(value.split("\t")[5]) >= int(tlen) or int(value.split("\t")[5]) <= -int(tlen)): # forward bitwise flags
                if not re.search(pattern="^\d+(S|H).*$",string=value.split("\t")[4]):
                    shift = sum(shift for shift in list(map(int,re.findall(string=''.join(re.findall(pattern="(\d+M|\d+D|\d+N|\d+P)",string=value.split("\t")[4])),pattern="\d+")))) # Defines the exact end position of the read.
                    if int(value.split("\t")[2])+shift <= int(ctx) and int(value.split("\t")[2])+shift >= int(ctx) - int(insert):
                        counter += 1
                else:
                    if int(value.split("\t")[2]) == int(ctx):
                        counter += 1
            elif int(value.split("\t")[0]) in reverse_flags and (int(value.split("\t")[5]) >= int(tlen) or int(value.split("\t")[5]) <= -int(tlen)): # reverse bitwise flags
                if not re.search(pattern="^.*\d+(S|H)$",string=value.split("\t")[4]):
                    if int(value.split("\t")[2]) >= int(ctx) and int(value.split("\t")[2]) <= int(ctx) + int(insert):
                        counter += 1
                else:
                    shift = sum(shift for shift in list(map(int,re.findall(string=''.join(re.findall(pattern="(\d+M|\d+D|\d+N|\d+P)",string=value.split("\t")[4])),pattern="\d+"))))
                    if int(value.split("\t")[2])+shift == int(ctx):
                        counter += 1
            elif int(value.split("\t")[0]) not in forward_flags and int(value.split("\t")[0]) not in reverse_flags:
                raise Exception(''.join (["Bitwise flags evaluation failed. The following flag is missing: ", value.split("\t")[0], ". Please, check if all the reads are paired."]))
    return counter

for i1,v1 in enumerate(final_hits):
    for i2,v2 in enumerate(v1[1]):
        for i3,v3 in enumerate(v2[1]):
            for c1 in mapped_lines_coords:
                if v1[0] == c1[0].split("\t")[0] and v3[0] == c1[0].split("\t")[1]:
                    sup_reads = supporting_reads(lst=c1,ctx=v2[0])
                    final_hits[i1][1][i2][1][i3][1] = ["\t".join([hit,str(sup_reads)]) for hit in v3[1]] # add number of supporting reads to each hit.
                    break

final_hits = [value for v1 in final_hits for v2 in v1[1] for v3 in v2[1] for value in v3[1]] # flatten the list.

def hit_sort(value): # sort hits by hit group size, chr1, ctx, no_hits, no_hits_secchr, the number of supporting reads, and overhang length in the reverse order.
    return (int(value.split(sep="\t")[16]),value.split(sep="\t")[2],list(map(int,value.split(sep="\t")[13:16])),int(value.split(sep="\t")[17]),int(value.split(sep="\t")[12]))
final_hits.sort(key=hit_sort,reverse=True)

for i1,v1 in enumerate(final_hits):
    l = v1.split("\t")
    l[-1],l[-2] = l[-2],l[-1]
    final_hits[i1]= "\t".join(l) # Reorder the hits placing the hit group size in the last column.

if no_filter is False:
    final_hits = [v for v in final_hits if int(v.split("\t")[16]) > int(v.split("\t")[15])] # keep only those hits where the number of supporting reads is higher than nohits_secchr.
len_final_filter = len(final_hits)

final_hits = [value for value in final_hits if (int(value.split("\t")[1]) in forward_flags and re.search(pattern="^.*\d+(S|H)$",string=value.split("\t")[5])) or (int(value.split("\t")[1]) in reverse_flags and re.search(pattern="^\d+(S|H).*$",string=value.split("\t")[5]))] # keep only those hits where CIGAR strings support the ctx position.
no_hits_with_supp_cigars = len(final_hits)

final_hits = [list(group) for name,group in groupby(final_hits, key=lambda x:(x.split("\t")[2],x.split("\t")[13]))] # regroup hits by chr1 and ctx.

no_final_hit_groups = len(final_hits)

final_hits = ["\t".join([str(i1+1),value]) for i1,v1 in enumerate(final_hits) for value in v1] # index and flatten the list.

print("Discovering hit pairs...")

for i1,v1 in enumerate(final_hits):
    if v1.split("\t")[7] == "=":
        l1 = v1.split("\t")
        l1[7] = l1[3]
        final_hits[i1] = "\t".join(l1) # replace '=' in the RNEXT field with the chromosome name.

final_hits = [[name,list(group)] for name,group in groupby(final_hits,key=lambda x:(x.split("\t")[0],x.split("\t")[3],x.split("\t")[14],x.split("\t")[7]))] # regroup hits by hit index, chr1, ctx position and chr2.

for i1,v1 in enumerate(final_hits): # discover hit pairs
    pnext_list = []
    for coords in mapped_lines_coords:
        coords_list = coords[0].split("\t")
        if coords_list[1] == "=":
            coords_list[1] = coords_list[0] # replace "=" with the chromosome name
        if v1[0][3] == coords_list[0] and v1[0][1] == coords_list[1]:
            for value1 in v1[1]:
                value1 = value1.split("\t")
                if int(value1[2]) in primary_alignments:
                    for value2 in coords[1]:
                        value2 = value2.split("\t")
                        if value1[1] == value2[6] and value1[8] == value2[2]:
                            if int(value2[0]) in forward_flags: # forward bitwise flags
                                shift = sum(shift for shift in list(map(int,re.findall(string=''.join(re.findall(pattern="(\d+M|\d+D|\d+N|\d+P)",string=value2[4])),pattern="\d+"))))-1 # Defines the exact end position of the mate read.
                                pnext = int(value2[2]) + shift
                            elif int(value2[0]) in reverse_flags: # reverse bitwise flags
                                pnext = int(value2[2])
                            else:
                                raise Exception(''.join (["Bitwise flags evaluation failed. The following flag is missing: ", value2[0], ". Please, check if all the reads are paired."]))
                            pnext_list.append(pnext)
                elif int(value1[2]) in secondary_alignments:
                    pnext_list.append(int(value1[8]))
                else:
                    raise Exception("Unable to determine whether the alignment is primary.")
    if len(pnext_list) == 0:
        raise Exception("PNEXT evaluation failed.")
    pnext_median = statistics.median(pnext_list)
    pairs = []
    for v2 in final_hits:
        if v2[0][1] == v1[0][3] and v1[0][1] == v2[0][3]:
            sec_idx = v2[0][0]
            distance = (pnext_median-int(v2[0][2]))
            counter_forward = 0
            counter_reverse = 0
            for value3 in v2[1]:
                value3 = value3.split("\t")
                if int(value3[2]) in forward_flags and re.search(pattern="^.*\d+(S|H)$",string=value3[6]):
                    counter_forward += 1
                elif int(value3[2]) in reverse_flags and re.search(pattern="^\d+(S|H).*$",string=value3[6]):
                    counter_reverse += 1
                else:
                    raise Exception("Identification of CTX orientation failed.")
            if counter_forward > 0 and counter_reverse == 0:
                if distance >= -int(insert) and distance < 0: # Distance has to be lower than 0, since CTX position is increased by 1.
                    res = ":".join([sec_idx,str(round(distance))])
                    pairs.append(res)
            elif counter_reverse > 0 and counter_forward == 0:
                if distance >= 0 and distance <= int(insert): # Distance may equal 0, since CTX position is not increased by 1.
                    res = ":".join([sec_idx,str(round(distance))])
                    pairs.append(res)
            else:
                if abs(distance) <= int(insert):
                    res = ":".join([sec_idx,str(round(distance)),"abs"])
                    pairs.append(res)
    if len(pairs) == 0:
        pairs.append("NA")
    final_hits[i1][1] = ["\t".join([hit,"; ".join(pairs)]) for hit in v1[1]]
del(mapped_lines_coords)

final_hits = [value for v1 in final_hits for value in v1[1]] # flatten the list

for index,value in enumerate(final_hits): # reorder each hit by placing paired hits next to the hit number.
    l = value.split("\t")
    l.insert(1,l[-1])
    l.pop()
    final_hits[index] = "\t".join(l)

# Hit_index is generated to determine the number of hit pairs.
def count_hit_pairs(lst):
    hit_index = ["\t".join(hit.split("\t")[0:2]) for hit in lst if hit.split("\t")[1] != "NA"]
    hit_index = [[name,list(group)] for name,group in groupby(hit_index,key=lambda x:x.split("\t")[0])]
    
    for i1,v1 in enumerate(hit_index):
        for value in v1[1]:
            value = value.split("\t")
            value.pop(0)
            hit_index[i1][1] = value
    
    for i1,v1 in enumerate(hit_index):
        for value in v1[1]:
            hit_index[i1][1] = re.findall(pattern="\d+\:[^a]", string=value)
    
    for i1,v1 in enumerate(hit_index):
        for i2,v2 in enumerate(v1[1]):
            hit_index[i1][1][i2] = v2.split(":")[0]
    
    hit_index = list(set(["\t".join([v1[0],v2]) for v1 in hit_index for v2 in v1[1]]))
    
    counter = 0
    for i1,v1 in enumerate(hit_index):
        for v2 in hit_index:
            if v1.split("\t")[0] == v2.split("\t")[1] and v1.split("\t")[1] == v2.split("\t")[0]:
                counter += 1
                lst[i1] = "\t".join([lst[i1], "TRUE"])
                break
    
    return int(counter/2)
no_final_hit_pairs = count_hit_pairs(final_hits)

for i1,v1 in enumerate(final_hits):
    if re.search(pattern="^\d+(S|H).*$",string=v1.split("\t")[7]):
        final_hits[i1] = v1
    elif re.search(pattern="^.*\d+(S|H)$",string=v1.split("\t")[7]):
        v1 = v1.split("\t")
        v1[15] = str(int(v1[15])-1) # Subtract 1 to get the exact CTX position.
        final_hits[i1] = '\t'.join(v1)
    else:
        raise Exception("CIGAR string evaluation failed.")

try:
    chr1, chr2
    print("Finding all hits for the specified pair of chromosomes...")
    final_hits_chr = [value for value in final_hits if (value.split("\t")[4] == chr1 and value.split("\t")[8] == chr2) or (value.split("\t")[4] == chr2 and value.split("\t")[8] == chr1)] # keep only mappings for the selected pair of chromosomes
    no_final_hit_groups_chr = len(set([hit.split("\t")[0] for hit in final_hits_chr]))
    no_final_hit_pairs_chr = count_hit_pairs(final_hits_chr)
except:
    None

header = "HIT_NO\tPAIRED_HITS\tQNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tRNEXT\tPNEXT\tTLEN\tSEQ\tSEQ_MATCHING\tSEQ_OVERHANG\tOVERHANG_LENGTH\tCTX_POS\tNO_HITS\tNO_HITS_PER_SEC_CHR\tNO_SUPP_READS\tHIT_GROUP_SIZE\tPAIRED_HIT"

if compressed_output:
    suffix = "_"+appname+"_results.tsv.gz"
else:
    suffix = "_"+appname+"_results.tsv"
outfile = prefix+suffix
print("Saving the results to a file: "+colored(outfile,"green"))
def print_output(out):
    out.write(appname+" analysis report\n\n")
    out.write("User-defined options used:\n")
    for opt,arg in opts:
        print(opt,arg, file=out, sep="\t")
    out.write("\nSummary:\n")
    print("Total number of reads:",str(all_reads),sep="\t",file=out)
    print("Number of reads after initial filtering:",str(len_mapped_lines),sep="\t",file=out)
    print("Number of reads after filtering out low-quality reads:",str(len_filtered_lines),sep="\t",file=out)
    print("Number of read groups after the CIGAR-based filtering:", str(len_cigar_filtered),sep="\t",file=out)
    print("Number of hits surviving the first filter (nohits):",str(len_nohits),sep="\t",file=out)
    print("Number of hits surviving the second filter (nohits_secchr):",str(len_nohits_secchr),sep="\t",file=out)
    print("Number of hits after the TLEN filtering:",str(len_tlen),sep="\t",file=out)
    print("Number of hits after filtering by the hit group size:",str(len_group_size),sep="\t",file=out)
    print("Number of hits after filtering by the number of supporting reads:",str(len_final_filter),sep="\t",file=out)
    print("Number of hits after filtering out hits with overhangs non-adjacent to the CTX:",str(no_hits_with_supp_cigars),sep="\t",file=out)
    print("Number of hit groups:",str(no_final_hit_groups),sep="\t",file=out)
    print("Number of hit pairs:",str(no_final_hit_pairs),sep="\t",file=out)

    try:
        chr1, chr2
        print("Number of hit groups after chromosome-based filtering:",str(no_final_hit_groups_chr),sep="\t",file=out)
        print("Number of hit pairs after chromosome-based filtering:",str(no_final_hit_pairs_chr),sep="\t",file=out)
    except:
        None
    out.write("\nSamtools header:\n")
    out.write(bamfile_header)
    try:
        chr1, chr2
        print("\n"+appname+" results ("+chr1+", "+chr2+"):",file=out)
        out.write(header+"\n")
        out.writelines('\n'.join(final_hits_chr)+"\n")
    except:
        None
    print("\n"+appname+" results (full table):",file=out)
    out.write(header+"\n")
    out.writelines('\n'.join(final_hits))

if compressed_output:
    with gzip.open(outfile, mode="wt") as out:
        print_output(out)
else:
    with open(outfile, mode="wt") as out:
        print_output(out)

print(colored("The analysis is complete.", "green"))
