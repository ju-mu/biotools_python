

from collections import OrderedDict
from operator import itemgetter


gtf_file="/media/BIOFS/Genomes/mm9/ensembl/gtf/Mus_musculus.NCBIM37.57.gtf"

def gtf2introns(gtf):
    intron_dict={}#container of the introns
    exon_dict={}#holds all exon related infos
    exoncount_dict={}#total number of exons per gene
    gord=OrderedDict()#stores the keys to prevent resorting
    lchrom=""
    for ec,exon in enumerate(open(gtf,"r")):
        exon=exon.strip("\n").split("\t")
        if(exon[2]!="exon"):continue
        start=int(exon[3])
        end=int(exon[4])
        exinfos=exon[8].split(" ")
        if(ec==0):gpos=[xc for xc,x in enumerate(exinfos) if x=='transcript_id'][0]
        gname=exinfos[gpos+1].strip("\"|;")
        chrom=exon[0]
        if(lchrom != chrom):
            intron_dict[chrom]={}
            exon_dict[chrom]=[]
            if(ec>0):exon_dict[lchrom]=sorted(exon_dict[lchrom], key=itemgetter(0))
            lchrom=chrom
            gord[chrom]=[]
        if gname not in exoncount_dict:
            exoncount_dict[gname]=0
            intron_dict[chrom][gname]={"start":[],"end":[],"info":exon}
            
        exoncount_dict[gname]+=1
        exon_dict[chrom].append([start,end,gname])

    cgenes=set()
    
    for chrom in exon_dict:
        for exon in exon_dict[chrom]:
            gname=exon[2]
            
            if gname not in cgenes :
                if exoncount_dict[gname]==1:
                    exoncount_dict[gname]=0
                    continue
                gord[chrom].append(gname)
                intron_dict[chrom][gname]["start"].append(exon[1]+1)#correct for opened interval               
                exoncount_dict[gname]-=1
                
            for cgene in cgenes:
                if(gname == cgene):
                    exoncount_dict[cgene]-=1
                    if(exoncount_dict[cgene]==0):
                        intron_dict[chrom][cgene]["end"].append(exon[0]-1)
                        continue
                tstarts=len(intron_dict[chrom][cgene]["start"])
                if exon[0]<=intron_dict[chrom][cgene]["start"][tstarts-1]:#correct for overlaps
                    intron_dict[chrom][cgene]["start"][tstarts-1]=max(intron_dict[chrom][cgene]["start"][tstarts-1],exon[1]+1)
                    continue
                intron_dict[chrom][cgene]["end"].append(exon[0]-1)
                intron_dict[chrom][cgene]["start"].append(exon[1]+1)
            cgenes.add(gname)
            if exoncount_dict[gname]==0:cgenes.remove(gname)        
    assert len(cgenes)==0,"Error in parsing!"
    print("Conversion completed")
    return intron_dict,gord
        

idict,gord=gtf2introns(gtf_file)

igtf=gtf_file.rstrip(".gtf")+"_intron.gtf"
igtf_fh=open(igtf,"w")
igl=0

for chrom in gord:
    for gene in gord[chrom]:
        gtfrow=idict[chrom][gene]["info"]
        gtfrow[2]="intron"
        exi=gtfrow[8].split(" ")
        exd=-1
        if "exon_number" in exi:
            exd=exi.index("exon_number")
            exi[exd]="intron_number"
            
        for exc in range(len(idict[chrom][gene]["start"])):
            if exd>-1:exi[exd+1]="\""+str(exc+1)+"\";"
            gtfrow[8]=" ".join(exi)
            gtfrow[3]=str(idict[chrom][gene]["start"][exc])
            gtfrow[4]=str(idict[chrom][gene]["end"][exc])
            igtf_fh.write("\t".join(gtfrow)+"\n")
            igl+=1

igtf_fh.close()
print("Successfully wrote %d lines to %s"%(igl,igtf))





