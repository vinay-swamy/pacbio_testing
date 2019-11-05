

def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        file.readline()#skip the first line because it has a header
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'paired':True if info[1]=='y' else False, 'tissue':info[2],'subtissue':info[3]}
    return(res)

def readPacBioSamples(sfile):
    res={}
    with open(sfile) as file:
        file.readline()#skip the first line because it has a header
        for line in file:
            info=line.strip('\n').split('\t')
            if info[0] not in res.keys():
                res[info[0]]={'tissue':info[1], 'bax':[info[2]]}
            else:
                res[info[0]]['bax'].append(info[2])
    return(res)


def subtissue_to_tissue(tissue, sample_dict):
    res=[]
    for sample in sample_dict.keys():
        pgtf='data/filtered_tissue_gtfs/{}.gtf'.format(sample_dict[sample]['subtissue'])
        if sample_dict[sample]['tissue'] == tissue and pgtf not in res:
            res.append(pgtf)
    return res
stdir=config['st_brain_outdir']
sample_file=config['sample_file']
tissue_file=config['tissueFile']
subtissue_file=config['subtissueFile']
sample_dict=readSampleFile(sample_file)
pb_samplefile=config['pacbiosamples']
pb_sd=readPacBioSamples(pb_samplefile)
pb_samplenames=pb_sd.keys()
with open(tissue_file) as tf, open(subtissue_file) as sf:
    tissues= [line.strip('\n') for line in tf]
    subtissues= [line.strip('\n') for line in sf]
ref_GTF=config['ref_GTF']
working_dir=config['working_dir']
rule all:
    input:expand('results/{tissue}.combined.gtf', tissue=tissues), expand('results/{tissue}.all_size_hg38.gff', tissue=tissues)#, expand('data/fasta/flnc/{sample}_flnc.fa', sample=pb_samplenames)

rule filter_gtf:
    input: gtf=stdir + '{st}.combined.gtf', track=stdir + '{st}.tracking'
    output: 'data/filtered_tissue_gtfs/{st}.gtf'
    shell:
        '''
        module load R
        Rscript filter_gtf.R {working_dir} {input.gtf} {ref_GTF} {input.track} {wildcards.st} {sample_file} {output}
        '''

rule merge_all_brain:
    input:lambda wildcards: subtissue_to_tissue(wildcards.tissue, sample_dict)
    output:'results/{tissue}.combined.gtf'
    shell:
        '''
        module load gffcompare
        gffcompare -r {ref_GTF} -p {wildcards.tissue} -o results/{wildcards.tissue} {input}

        '''

rule bax2bam:
    input:lambda wildcards: expand('data/raw/{path}', path=pb_sd[wildcards.sample]['bax'])
    output: 'data/bams/{sample}.subreads.bam'
    shell:
        '''
        module load smrtanalysis
        bax2bam -o data/bams/{wildcards.sample} {input}
        '''

rule bam2ccs:
    input:'data/bams/{sample}.subreads.bam'
    output:css='data/ccs/{sample}_css.bam', xml='data/ccs/{sample}_css.xml'
    shell:
        '''
        module load smrtanalysis
        ccs --noPolish --minLength=200 --minPasses=1 --minZScore=-999 \
          --maxDropFraction=0.8 --minPredictedAccuracy=0.8 --minSnr=4 {input} {output.css}
         dataset create --type ConsensusReadSet {output.xml} {output.css}
        '''



rule css2fasta:
    input: xml='data/ccs/{sample}_css.xml'
    output:draft='data/fasta/draft/{sample}_draft.fa', flnc='data/fasta/flnc/{sample}_flnc.fa', nfl='data/fasta/nfl/{sample}_nfl.fa'
    shell:
        '''
        module load smrtanalysis
        pbtranscript classify {input.xml} {output.draft}  --flnc={output.flnc} --nfl={output.nfl}

        '''



rule liftOver_pacbio_gff:
    input:'data/{tissue}.all_size.5merge.collapsed.gff'
    output:'results/{tissue}.all_size_hg38.gff'
    shell:
        '''
        module load crossmap
        crossmap gff ref/hg19ToHg38.over.chain.gz {input} {output}
        '''
