#!/usr/bin/env python3

import sys,re,gzip

#Making a function that returns the relative abundance of the contig the subtype is found on
def findAbundance(contigName):
        abundance=0
        #Opening profiles made with msamtools:
        with gzip.open(sys.argv[1],'rt') as infile:
                for abundanceProfile in infile:
                        #If the contigname is found in the profile, the relative abundance is saved:
                        if contigName in abundanceProfile:
                                abundance=float(re.search(r'\s+([-\w.]+)\n', abundanceProfile).group(1))
                                return(abundance)


#Defining each CRISPR type as a variable
(IA, IB, IC, ID, IE, IF, IFT, IG,
IIA, IIB, IIC,
IIIA, IIIB, IIIC, IIID, IIIE, IIIF,
IVA1, IVA2, IVA3, IVB, IVC, IVD, IVE,
VA, VB, VC, VD, VE, VF, VG, VH, VI, VJ, VK,
VIA, VIB, VIC, VID)=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Looking for each subtype and adding to its relative abundance with regular expressions:
for line in sys.stdin:
        if re.search(r'\s+(I\WA\n|Hybrid\WI\WA,|\S+,I\WA\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IA+=foundAbundance
        if re.search(r'\s+(I\WB\n|Hybrid\WI\WB,|\S+,I\WB\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IB+=foundAbundance
        if re.search(r'\s+(I\WC\n|Hybrid\WI\WC,|\S+,I\WC\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IC+=foundAbundance
        if re.search(r'\s+(I\WD\n|Hybrid\WI\WD,|\S+,I\WD\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                ID+=foundAbundance
        if re.search(r'\s+(I\WE\n|Hybrid\WI\WE,|\S+,I\WE\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IE+=foundAbundance
        if re.search(r'\s+(I\WF\n|Hybrid\WI\WF,|\S+,I\WF\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IF+=foundAbundance
        if re.search(r'\s+(I\WF_T\n|Hybrid\WI\WF_T,|\S+,I\WF_T\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IFT+=foundAbundance
        if re.search(r'\s+(I\WG\n|Hybrid\WI\WG,|\S+,I\WG\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IG+=foundAbundance
        if re.search(r'\s+(II\WA\n|Hybrid\WII\WA,|\S+,II\WA\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIA+=foundAbundance
        if re.search(r'\s+(II\WB\n|Hybrid\WII\WB,|\S+,II\WB\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIB+=foundAbundance
        if re.search(r'\s+(II\WC\n|Hybrid\WII\WC,|\S+,II\WC\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIC+=foundAbundance
        if re.search(r'\s+(III\WA\n|Hybrid\WIII\WA,|\S+,III\WA\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIIA+=foundAbundance
        if re.search(r'\s+(III\WB\n|Hybrid\WIII\WB,|\S+,III\WB\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIIB+=foundAbundance
        if re.search(r'\s+(III\WC\n|Hybrid\WIII\WC,|\S+,III\WC\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIIC+=foundAbundance
        if re.search(r'\s+(III\WD\n|Hybrid\WIII\WD,|\S+,III\WD\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIID+=foundAbundance
        if re.search(r'\s+(III\WE\n|Hybrid\WIII\WE,|\S+,III\WE\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIIE+=foundAbundance
        if re.search(r'\s+(III\WF\n|Hybrid\WIII\WF,|\S+,III\WF\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IIIF+=foundAbundance
        if re.search(r'\s+(IV\WA1\n|Hybrid\WIV\WA1,|\S+,IV\WA1\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IVA1+=foundAbundance
        if re.search(r'\s+(IV\WA2\n|Hybrid\WIV\WA2,|\S+,IV\WA2\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IVA2+=foundAbundance
        if re.search(r'\s+(IV\WA3\n|Hybrid\WIV\WA3,|\S+,IV\WA3\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IVA3+=foundAbundance
        if re.search(r'\s+(IV\WB\n|Hybrid\WIV\WB,|\S+,IV\WB\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IVB+=foundAbundance
        if re.search(r'\s+(IV\WC\n|Hybrid\WIV\WC,|\S+,IV\WC\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IVC+=foundAbundance
        if re.search(r'\s+(IV\WD\n|Hybrid\WIV\WD,|\S+,IV\WD\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IVD+=foundAbundance
        if re.search(r'\s+(IV\WE\n|Hybrid\WIV\WE,|\S+,IV\WE\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                IVE+=foundAbundance
        if re.search(r'\s+(V\WA\n|Hybrid\WV\WA,|\S+,V\WA\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VA+=foundAbundance
        if re.search(r'\s+(V\WB\n|Hybrid\WV\WB,|\S+,V\WB\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VB+=foundAbundance
        if re.search(r'\s+(V\WC\n|Hybrid\WV\WC,|\S+,V\WC\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VC+=foundAbundance
        if re.search(r'\s+(V\WD\n|Hybrid\WV\WD,|\S+,V\WD\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VD+=foundAbundance
        if re.search(r'\s+(V\WE\n|Hybrid\WV\WE,|\S+,V\WE\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VE+=foundAbundance
        if re.search(r'\s+(V\WF\n|Hybrid\WV\WF,|\S+,V\WF\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VF+=foundAbundance
        if re.search(r'\s+(V\WG\n|Hybrid\WV\WG,|\S+,V\WG\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VG+=foundAbundance
        if re.search(r'\s+(V\WH\n|Hybrid\WV\WH,|\S+,V\WH\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VH+=foundAbundance
        if re.search(r'\s+(V\WI\n|Hybrid\WV\WI,|\S+,V\WI\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VI+=foundAbundance
        if re.search(r'\s+(V\WJ\n|Hybrid\WV\WJ,|\S+,V\WJ\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VJ+=foundAbundance
        if re.search(r'\s+(V\WK\n|Hybrid\WV\WK,|\S+,V\WK\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VK+=foundAbundance
        if re.search(r'\s+(VI\WA\n|Hybrid\WVI\WA,|\S+,VI\WA\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VIA+=foundAbundance
        if re.search(r'\s+(VI\WB\n|Hybrid\WVI\WB,|\S+,VI\WB\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VIB+=foundAbundance
        if re.search(r'\s+(VI\WC\n|Hybrid\WVI\WC,|\S+,VI\WC\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VIC+=foundAbundance
        if re.search(r'\s+(VI\WD\n|Hybrid\WVI\WD,|\S+,VI\WD\W)',line):
                contig=re.search(r'^([-\w.]+)\s', line).group(1)
                foundAbundance=findAbundance(contig)
                VID+=foundAbundance


#Setting name/ID of sample
SampleName=0
SampleName=re.search(r'results/(\S+)',sys.argv[2]).group(1)

#Printing line with results:
print("S"+SampleName+"\t",IA,"\t",IB,"\t",IC,"\t",ID,"\t",IE,"\t",IF,"\t",IFT,"\t",IG,"\t",
IIA,"\t",IIB,"\t",IIC,"\t",
IIIA,"\t",IIIB,"\t",IIIC,"\t",IIID,"\t",IIIE,"\t",IIIF,"\t",
IVA1,"\t",IVA2,"\t",IVA3,"\t",IVB,"\t",IVC,"\t",IVD,"\t",IVE,"\t",
VA,"\t",VB,"\t",VC,"\t",VD,"\t",VE,"\t",VF,"\t",VG,"\t",VH,"\t",VI,"\t",VJ,"\t",VK,"\t",
VIA,"\t",VIB,"\t",VIC,"\t",VID)
