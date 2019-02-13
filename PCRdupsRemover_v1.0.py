#!/usr/bin/env python2


import itertools
import gzip
import regex
import sys
import os
from optparse import OptionParser


""" Parsing command line options """
parser = OptionParser(prog="PCRdups_remover", usage="%prog -1 filename_R1 -2 filename_R2 [-n] [-d motif] [-m num] [-s num] [-l num] [-o foldername]"
                      , version="%prog 1.1")
parser.add_option('-1', '--r1', action="store", type="string", dest="file_r1",
                  help="Name of the R1 file\n")
parser.add_option('-2', '--r2', action="store", type="string", dest="file_r2",
                  help="Name of the R2 file\n")
parser.add_option('-n', '--no-dbr-match', action="store_true", dest="no_dbr_match", default=False,
                  help="Skip DBR motif matching, if specified, options -m and -s are ignored")
parser.add_option('-d', '--dbr-motif', action="store", type="string", dest="dbr_motif", default="NNNNNHMMGG",
                  help="DBR motif, default=NNNNNHMMGG\n")
parser.add_option('-m', '--mismatches', action="store", type="string", dest="dbr_mm", default="1",
                  help="Number of mismatches allowed in the DBR, default=1\n")
parser.add_option('-s', '--shift', action="store", type="int", dest="dbr_shift", default=2,
                  help="Number of positions to shift the DBR to the right to attempt matching, default=2\n")
parser.add_option('-l', '--dna-length', action="store", type="int", dest="dna_len", default=50,
                  help="Number of basepairs (excluding the DBR) to compare and look for "
                       "PCR duplications, default=50\n")
parser.add_option('-o', '--outfolder', action="store", type="string", dest="outfolder", default="output",
                  help="Name of the folder to write the results to, default=output\n")
(options, args) = parser.parse_args()


""" Translation of ambiguities to regular expressions for DBR comparison """
ambiguities = {"A":"A",
               "T":"T",
               "C":"C",
               "G":"G",
               "R":"[AG]",
               "Y":"[CT]",
               "K":"[GT]",
               "M":"[AC]",
               "S":"[CG]",
               "W":"[AT]",
               "B":"[CGT]",
               "D":"[AGT]",
               "H":"[ACT]",
               "V":"[ACG]",
               "N":"[ACGTN]"}


""" Checking for mandatory options """
if not all([options.file_r1, options.file_r2]):
    print "\n\tMust include -1 (--r1), and -2 (--r2) options, see -h for help\n"
    sys.exit()


""" Pass options from parser to variables """
file_r1 = options.file_r1
file_r2 = options.file_r2
file_r1_basename = options.file_r1.split(".")[0]
file_r2_basename = options.file_r2.split(".")[0]
no_dbr_match = options.no_dbr_match
dbr_motif = list(options.dbr_motif)
dbr_len = len(options.dbr_motif)
dbr_mm = options.dbr_mm
dbr_shift = options.dbr_shift
dna_len = options.dna_len
outfolder = "./" + options.outfolder + "/"


""" Create output folder if not already created """
if not os.path.exists(outfolder):
    os.makedirs(outfolder)


""" Load fastq files 4 lines at a time """
if ".gz" in file_r1: # if files are compressed
    f1 = gzip.open(file_r1, "rb")
    f2 = gzip.open(file_r2, "rb")
else: # if files are not compressed
    f1 = open(file_r1, "r")
    f2 = open(file_r2, "r")
k1 = itertools.izip(*[iter(f1)]*4)
k2 = itertools.izip(*[iter(f2)]*4)


""" Compile regular expression for the DBR motif, regex library needed """
dbr_regex = "".join([ambiguities[nuc] for nuc in dbr_motif])
if dbr_mm == "0": # when 0 mismatches in DBR allowed
    dbr = regex.compile(dbr_regex)
else: # when dbr_mm mismatches are allowed in the DBR
    dbr_regex = "("+dbr_regex+"){s<="+dbr_mm+"}"
    dbr = regex.compile(dbr_regex)


""" Start processing of files """
print "\n\n\t====================================="
print     "\t Remover of PCR-duplicates for ddRAD"
print     "\t=====================================\n"

print "\n\tProcessing "+ file_r1 + " and " + file_r2 + "..."

pcr_dups = {} # Dictionary to contain duplicated portions of reads
dbrs = {}     # Dictionary to contain DBR motifs
seq_len = dbr_len + dna_len
r1_good = []
r2_good = []
r1_dups = []
r2_dups = []
r1_bad = []
r2_bad = []

while 1:
    try: r2 = list(k2.next())
    except StopIteration: break
    r1 = list(k1.next())
    dbr_seq = r2[1][:dbr_len]
    seqout = r2[1][dbr_len:]
    qualout = r2[3][dbr_len:]
    if no_dbr_match == False:
        if dbr.match(dbr_seq):
            seq_uniq = r2[1][:seq_len]
            if not seq_uniq in pcr_dups:
                if not dbr_seq in dbrs:
                    dbrs[dbr_seq] = 1
                else:
                    dbrs[dbr_seq] += 1
                pcr_dups[seq_uniq] = 1
                r2[1] = seqout
                r2[3] = qualout
                r2_good.append("".join(r2))
                r1_good.append("".join(r1))
            else:
                pcr_dups[seq_uniq] += 1
                r2[1] = seqout
                r2[3] = qualout
                r1_dups.append("".join(r1))
                r2_dups.append("".join(r2))
            continue
        else:
            if dbr_shift == 0:
                r1_bad.append("".join(r1))
                r2_bad.append("".join(r2))
                continue
            else:
                if dbr_shift > 0:
                    for s in range(1,dbr_shift+1):
                        dbr_seq = r2[1][s:dbr_len+s]
                        seqout = r2[1][dbr_len+s:]
                        qualout = r2[3][dbr_len+s:]
                        if dbr.match(dbr_seq):
                            seq_uniq = r2[1][s:seq_len+s]
                            if not seq_uniq in pcr_dups:
                                if not dbr_seq in dbrs:
                                    dbrs[dbr_seq] = 1
                                else:
                                    dbrs[dbr_seq] += 1
                                pcr_dups[seq_uniq] = 1
                                r2[1] = seqout
                                r2[3] = qualout
                                r2_good.append("".join(r2))
                                r1_good.append("".join(r1))
                            else:
                                pcr_dups[seq_uniq] += 1
                                r2[1] = seqout
                                r2[3] = qualout
                                r1_dups.append("".join(r1))
                                r2_dups.append("".join(r2))
                            break
                        else:
                            if s == dbr_shift:
                                r1_bad.append("".join(r1))
                                r2_bad.append("".join(r2))
                            continue
    elif no_dbr_match == True:
        seq_uniq = r2[1][:seq_len]
        if not seq_uniq in pcr_dups:
            pcr_dups[seq_uniq] = 1
            r2[1] = seqout
            r2[3] = qualout
            r2_good.append("".join(r2))
            r1_good.append("".join(r1))
        else:
            pcr_dups[seq_uniq] += 1
            r2[1] = seqout
            r2[3] = qualout
            r1_dups.append("".join(r1))
            r2_dups.append("".join(r2))
        if dbr.match(dbr_seq):
            if not dbr_seq in dbrs:
                dbrs[dbr_seq] = [1,1]
            else:
                dbrs[dbr_seq][0] += 1
        else:
            if not dbr_seq in dbrs:
                dbrs[dbr_seq] = [1,0]
            else:
                dbrs[dbr_seq][0] += 1


""" Create output files """
out_r1_good = gzip.open(outfolder+file_r1_basename+".fq.gz", "wb")
out_r1_good.write("".join(r1_good))
out_r1_good.close()
out_r2_good = gzip.open(outfolder+file_r2_basename+".fq.gz", "wb")
out_r2_good.write("".join(r2_good))
out_r2_good.close()

out_r1_dups = gzip.open(outfolder+file_r1_basename+"_duplicates.fq.gz", "wb")
out_r1_dups.write("".join(r1_dups))
out_r1_dups.close()
out_r2_dups = gzip.open(outfolder+file_r2_basename+"_duplicates.fq.gz", "wb")
out_r2_dups.write("".join(r2_dups))
out_r2_dups.close()

if no_dbr_match == True:
    out_r1_bad_dbr = gzip.open(outfolder+file_r1_basename+"_bad_DBRs.fq.gz", "wb")
    out_r1_bad_dbr.write("".join(r1_bad))
    out_r1_bad_dbr.close()
    out_r2_bad_dbr = gzip.open(outfolder+file_r2_basename+"_bad_DBRs.fq.gz", "wb")
    out_r2_bad_dbr.write("".join(r2_bad))
    out_r2_bad_dbr.close()

out_report = open(outfolder+file_r2_basename+"_report.txt", "w")
if no_dbr_match == True:
    report = ['Accepted reads:\t'+str(len(r1_good))+' including unmatched DBRs']
else:
    report = ['Accepted reads:\t'+str(len(r1_good))]

report.append('Duplicated reads:\t'+str(len(r1_dups)))

if no_dbr_match == False:
    report.append('Reads with unmatched DBR:\t'+str(len(r1_bad)))

report.append('\nPCR duplicate\tTotal')
pcr_dups_sorted = sorted(pcr_dups.items(), key=lambda x:x[1], reverse=True)
for dup in pcr_dups_sorted:
    if dup[1] > 4:
        report.append(dup[0] + "\t" + str(dup[1]))

report.append('\nPotential Repeat\tTotal')
repeats = {}
for dup in pcr_dups:
    seq = dup[dbr_len:]
    if not seq in repeats:
        repeats[seq] = pcr_dups[dup]
    else:
        repeats[seq] += pcr_dups[dup]
repeats_sorted = sorted(repeats.items(), key=lambda x:x[1], reverse=True)
for repeat in repeats_sorted:
    if repeat[1] > 4:
        report.append(repeat[0] + "\t" + str(repeat[1]))

if no_dbr_match == False:
    report.append('\nNo.\tDBR\tTotal')
    dbrs_sorted = sorted(dbrs.items(), key=lambda x:x[1], reverse=True)
    DBRnumber = 1
    for motif in dbrs_sorted:
        report.append(str(DBRnumber) + "\t" + motif[0] + "\t" + str(motif[1]))
        DBRnumber += 1
else:
    report.append('\nNo.\tDBR\tTotal\tUnmatched')
    dbrs_sorted = sorted(dbrs.items(), key=lambda x:x[1][0], reverse=True)
    DBRnumber = 1
    for motif in dbrs_sorted:
        report.append(str(DBRnumber) + "\t" + motif[0] + "\t" + str(motif[1][0]) + ("\t*" if motif[1][1]==0 else "\t"))
        DBRnumber += 1

out_report.write("\n".join(report))
out_report.close()

print "\tProcessing of "+ file_r1 + " and " + file_r2 + " complete"
