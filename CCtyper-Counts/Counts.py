#!/usr/bin/env python3
import sys,re

#Defining each CRISPR-Cas subtype as a variable
(IA, IB, IC, ID, IE, IF, IFT, IG,
IIA, IIB, IIC,
IIIA, IIIB, IIIC, IIID, IIIE, IIIF,
IVA1, IVA2, IVA3, IVB, IVC, IVD, IVE,
VA, VB, VC, VD, VE, VF, VG, VH, VI, VJ, VK,
VIA, VIB, VIC, VID, unk)=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

#Using regular expressions to look for each subtype:
for line in sys.stdin:
        if re.search(r'(^|\W)I\WA',line):
                IA+=1
        if re.search(r'(^|\W)I\WB',line):
                IB+=1
        if re.search(r'(^|\W)I\WC',line):
                IC+=1
        if re.search(r'(^|\W)I\WD',line):
                ID+=1
        if re.search(r'(^|\W)I\WE',line):
                IE+=1
        if re.search(r'(^|\W)I\WF',line):
                IF+=1
        if re.search(r'(^|\W)I\WF_T',line):
                IFT+=1
        if re.search(r'(^|\W)I\WG',line):
                IG+=1
        if re.search(r'(^|\W)II\WA',line):
                IIA+=1
        if re.search(r'(^|\W)II\WB',line):
                IIB+=1
        if re.search(r'(^|\W)II\WC',line):
                IIC+=1
        if re.search(r'(^|\W)III\WA',line):
                IIIA+=1
        if re.search(r'(^|\W)III\WB',line):
                IIIB+=1
        if re.search(r'(^|\W)III\WC',line):
                IIIC+=1
        if re.search(r'(^|\W)III\WD',line):
                IIID+=1
        if re.search(r'(^|\W)III\WE',line):
                IIIE+=1
        if re.search(r'(^|\W)III\WF',line):
                IIIF+=1
        if re.search(r'(^|\W)IV\WA1',line):
                IVA1+=1
        if re.search(r'(^|\W)IV\WA2',line):
                IVA2+=1
        if re.search(r'(^|\W)IV\WA3',line):
                IVA3+=1
        if re.search(r'(^|\W)IV\WB',line):
                IVB+=1
        if re.search(r'(^|\W)IV\WC',line):
                IVC+=1
        if re.search(r'(^|\W)IV\WD',line):
                IVD+=1
        if re.search(r'(^|\W)IV\WE',line):
                IVE+=1
        if re.search(r'(^|\W)V\WA',line):
                VA+=1
        if re.search(r'(^|\W)V\WB',line):
                VB+=1
        if re.search(r'(^|\W)V\WC',line):
                VC+=1
        if re.search(r'(^|\W)V\WD',line):
                VD+=1
        if re.search(r'(^|\W)V\WE',line):
                VE+=1
        if re.search(r'(^|\W)V\WF',line):
                VF+=1
        if re.search(r'(^|\W)V\WG',line):
                VG+=1
        if re.search(r'(^|\W)V\WH',line):
                VH+=1
        if re.search(r'(^|\W)V\WI',line):
                VI+=1
        if re.search(r'(^|\W)V\WJ',line):
                VJ+=1
        if re.search(r'(^|\W)V\WK',line):
                VK+=1
        if re.search(r'(^|\W)VI\WA',line):
                VIA+=1
        if re.search(r'(^|\W)VI\WB',line):
                VIB+=1
        if re.search(r'(^|\W)VI\WC',line):
                VIC+=1
        if re.search(r'(^|\W)VI\WD',line):
                VID+=1
        if re.search(r'(^|\W)Unknown',line):
                unk+=1

#Setting name of sample
SampleName=0
SampleName=re.search(r'results/(\S+)',sys.argv[1]).group(1)

#Printing results in a line:
print("S"+SampleName+"\t",IA,"\t",IB,"\t",IC,"\t",ID,"\t",IE,"\t",IF,"\t",IFT,"\t",IG,"\t",
IIA,"\t",IIB,"\t",IIC,"\t",
IIIA,"\t",IIIB,"\t",IIIC,"\t",IIID,"\t",IIIE,"\t",IIIF,"\t",
IVA1,"\t",IVA2,"\t",IVA3,"\t",IVB,"\t",IVC,"\t",IVD,"\t",IVE,"\t",
VA,"\t",VB,"\t",VC,"\t",VD,"\t",VE,"\t",VF,"\t",VG,"\t",VH,"\t",VI,"\t",VJ,"\t",VK,"\t",
VIA,"\t",VIB,"\t",VIC,"\t",VID,"\t", unk)

