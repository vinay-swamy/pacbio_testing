

def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        file.readline()#skip the first line because it has a header
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'tissue':info[1],'subtissue':info[2],'origin':info[3], 'filename':info[4]}
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

def build_to_gtf(tissue, build):
    if build == 'gc':
        gtf='ref/gencode_comp_ano.gtf'
    else:
        gtf=f'results/{tissue}.combined.gtf'
    return gtf


stdir=config['st_brain_outdir']
sample_file=config['sample_file']
tissue_file=config['tissueFile']
subtissue_file=config['subtissueFile']
sample_dict=readSampleFile(sample_file)
sample_names=sample_dict.keys()
fastq_path=config['fastq_path']
with open(tissue_file) as tf, open(subtissue_file) as sf:
    tissues= [line.strip('\n') for line in tf]
    subtissues= [line.strip('\n') for line in sf]
ref_GTF=config['ref_GTF']
ref_genome=config['ref_genome']
working_dir=config['working_dir']
variants=config['variants']
chr_names=config['chr_names']

rule all:
    # , expand('data/fasta/flnc/{sample}_flnc.fa', sample=pb_samplenames)
    input: #expand('results/all_{tissue}.merged.gtf', tissue=tissues)
        expand('results/{subtissue}.{build}_talon.gtf', subtissue=subtissues, build = ['raw', 'clean'])


# rule clean_genome:
#     input: ref_genome
#     output: 'ref/gencode_genome_clean.fa'
#     shell:
#         '''
#         python3 scripts/filterFasta.py {input} {chr_names} {output}
#         '''


rule build_minimap_index:
    input: ref_genome
    output: 'ref/minimap2_index.mmi'
    shell:
        '''
        module load minimap2
        minimap2 -d {output} {input}
        '''


rule build_talon_db:
    input: 
        gtf=ref_GTF, 
        genome=ref_genome
    params: 
        anno_name=lambda wildcards: f'db_{wildcards.tissue}', 
        prefix= lambda wildcards: f'{wildcards.tissue}', 
        outdir= lambda wildcards: f'data/talon_db/{wildcards.tissue}'
    output: 'data/talon_db_{build}/{tissue}.db'
    shell:
        '''
        talon_initialize_database --f {input.gtf} --g {input.genome} --a {params.anno_name} --idprefix {params.prefix} --o {params.outdir} 
        '''

rule align_minimap2:
    input:
        fa=lambda  wildcards: fastq_path + sample_dict[wildcards.sample]['filename'], 
        idx='ref/minimap2_index.mmi', 
        genome=ref_genome
    output:
        sam='data/raw_sams/{tissue}/{sample}.sam'
    shell:
        '''
        module load minimap2 
        module load samtools 
        minimap2 -ax splice -uf --secondary=no -C5 {input.genome} {input.fa} | samtools sort -O sam - > {output.sam}     
        '''

rule make_spliceJN_file:
    input:
        gtf=ref_GTF, 
        genome='ref/gencode_genome_clean.fa'
    output:
        'data/SJ_tab.txt'
    shell:
        '''
        module load samtools 
        module load bedtools
        python TranscriptClean/accessory_scripts/get_SJs_from_gtf.py --f {input.gtf} --g {input.genome} --o {output}

        '''



rule TranscriptClean: # pull from github repo
    input: 
        sam = 'data/raw_sams/{tissue}/{sample}.sam', 
        sj = 'data/SJ_tab.txt',
        genome = 'ref/gencode_genome_clean.fa'
    params:
        out_dir=lambda wildcards: f'data/clean_sam/{wildcards.tissue}'
    output:
        'data/clean_sams/{tissue}/{sample}.sam'
    shell:
        '''
        python TranscriptClean/TranscriptClean.py --sam {input.sam} \
            --genome {input.genome}  \
            --spliceJns {input.sj}\
            --variants {variants}\
            --maxLenIndel 5 \
            --maxSJOffset 5 \
            -m true \
            -i true \
            --correctSJs true \
            --primaryOnly \
            --outprefix {params.out_dir} \
            --threads 16 
        '''


'''
JK build is the genome
'''



rule run_talon:
    input: 
        sam = lambda wildcards: [f'data/{wildcards.build}_sams/{wildcards.tissue}/{sample}.sam' for sample in sample_names if sample_dict[sample]['subtissue'] == wildcards.tissue], 
        db = 'data/talon_db_{build}/{tissue}.db', 
        genome = 'ref/gencode_genome_clean.fa'
    params: 
        build='gencode', 
        outdir=lambda wildcards: f'data/talon_results_{wildcards.build}/{wildcards.tissue}'
    output: 
        anno = 'data/talon_results_{build}/{tissue}_talon_read_annot.tsv', 
        qc = 'data/talon_results_{build}/{tissue}_QC.log'
    shell:
        '''
        tlncfg=/tmp/config.csv
        echo pac_bio_{wildcards.tissue},whole_{wildcards.tissue},pacbio,{input.sam} > $tlncfg
        talon --f $tlncfg --db {input.db} --build {input.genome} --threads 16 --tmpdir {wildcards.build}_{wildcards.tissue}_talon_tmp/ --o {params.outdir}
    
        rm -rf {wildcards.build}_{wildcards.tissue}
        '''



rule make_gtf_from_talon:
    input: 
        anno='data/talon_results_{build}/{tissue}_talon_read_annot.tsv', 
        db='data/talon_db_{build}/{tissue}.db', 
        genome=ref_genome
    params:
        prefix= lambda wildcards: f'{wildcards.tissue}_{wildcards.build}', 
        outfix=lambda wildcards: f'results/{wildcards.tissue}.{wildcards.build}'
    output:
        'results/{tissue}.{build}_talon.gtf'
    shell:
        '''
        talon_create_GTF --db {input.db} --build {input.genome} --annot {params.prefix} --o {params.outfix}
        '''


# rule merge_pacbio_stringtie_gtfs:
#     input: pb=expand('results/{{tissue}}.{build}_talon.gtf', build=['dn', 'gc']), st='results/{tissue}.combined.gtf'
#     params: outdir=lambda wildcards: f'results/all_{wildcards.tissue}'
#     output: 'results/all_{tissue}.merged.gtf'
#     shell:
#         '''
#         module load gffcompare 
#         gffcompare -r {ref_GTF} -o {params.outdir} {input.st} {input.pb}
#         mv results/all_{wildcards.tissue}.combined.gtf {output}
#         '''




'''
These are rules to try and convert the raw files to fasta, but it wasnt working and I dont care enough to try 

'''




# rule bax2bam:
#     input:lambda wildcards: expand('data/raw/{path}', path=pb_sd[wildcards.sample]['bax'])
#     output: 'data/bams/{sample}.subreads.bam'
#     shell:
#         '''
#         module load smrtanalysis
#         bax2bam -o data/bams/{wildcards.sample} {input}
#         '''

# rule bam2ccs:
#     input:'data/bams/{sample}.subreads.bam'
#     output:css='data/ccs/{sample}_css.bam', xml='data/ccs/{sample}_css.xml'
#     shell:
#         '''
#         module load smrtanalysis
#         ccs --noPolish --minLength=200 --minPasses=1 --minZScore=-999 \
#           --maxDropFraction=0.8 --minPredictedAccuracy=0.8 --minSnr=4 {input} {output.css}
#          dataset create --type ConsensusReadSet {output.xml} {output.css}
#         '''



# rule css2fasta:
#     input: xml='data/ccs/{sample}_css.xml'
#     output:draft='data/fasta/draft/{sample}_draft.fa', flnc='data/fasta/flnc/{sample}_flnc.fa', nfl='data/fasta/nfl/{sample}_nfl.fa'
#     shell:
#         '''
#         module load smrtanalysis
#         pbtranscript classify {input.xml} {output.draft}  --flnc={output.flnc} --nfl={output.nfl}

#         '''
