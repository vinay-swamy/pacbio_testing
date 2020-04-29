def readPacBioSamples(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res={}
    with open(samplefile) as file:
        file.readline()#skip the first line because it has a header
        for line in file:
            info=line.strip('\n').split('\t')
            res[info[0]]={'tissue':info[1],'subtissue':info[2],'origin':info[3], 'filename':info[4]}
    return(res)


def readSampleFile(samplefile):
    # returns a dictionary of dictionaries where first dict key is sample id and second dict key are sample  properties
    res = {}
    with open(samplefile) as file:
        file.readline()  # skip the first line because it has a header
        for line in file:
            info = line.strip('\n').split('\t')
            res[info[0]] = {'paired': True if info[1] == 'y' else False, 'tissue': info[2],
                            'subtissue': info[3], 'origin': info[4], 'body_location': info[5], 'lib_type': info[6]}
    return(res)


def make_talon_config(wc, sample_dict):
    res = [f'{sample},{wc},pacbio,data/clean_sams/RPE_Fetal.Tissue/{sample}_clean.sam' for sample in sample_dict.keys()
           if sample_dict[sample]['subtissue'] == wc]
    return('\n'.join(res))


def get_stringtie_libflag(sample, sample_dict):
    if sample_dict[sample]['lib_type'] == 'first':
        return '--rf'
    elif sample_dict[sample]['lib_type'] == 'second':
        return '--fr'
    else:
        return ''


def salmon_input(id,sample_dict,fql):
    paired=sample_dict[id]['paired']
    id= fql + 'fastq_files/' + id
    if paired:
        return('-1 {s}_1.fastq.gz -2 {s}_2.fastq.gz'.format(s=id))
    else:
        return('-r {}.fastq.gz'.format(id))


def build2gtf(wc):
    if wc == 'gencode_comp_ano':
        return ref_GTF
    elif wc == 'all_RPE_loose':
        return 'data/combined_gtfs/all_RPE_loose.combined.gtf'
    elif wc == 'all_RPE_strict':
        return 'data/combined_gtfs/all_RPE_strict.combined.gtf'   

def build2fa(wc):
    if wc == 'gencode_comp_ano':
        return 'data/seqs/gencode_comp_ano_tx.fa'
    elif wc == 'all_RPE_loose':
        return 'data/seqs/all_RPE_loose_tx.fa'    


########## long read sample setup
working_dir = config['working_dir']
lr_sample_file = config['lr_sample_file']
tissue_file=config['tissueFile']
subtissue_file=config['subtissueFile']
lr_sample_dict = readPacBioSamples(lr_sample_file)
lr_sample_names=lr_sample_dict.keys()
lr_fql=config['lr_fql']
with open(tissue_file) as tf, open(subtissue_file) as sf:
    tissues= [line.strip('\n') for line in tf]
    subtissues= [line.strip('\n') for line in sf]
ref_GTF=config['ref_GTF']
ref_genome=config['ref_genome']

working_dir=config['working_dir']
variants=config['variants']
chr_names=config['chr_names']
uniprot_file=config['uniprot_file']
############ shot read sample set up 
sr_sample_file = config['sr_sampleFile']
sr_sample_dict = readSampleFile(config['sr_sampleFile'])
sr_sample_names = sr_sample_dict.keys()
sr_fql = config['sr_fql']
bam_path = config['bam_path']
STARindex = bam_path + 'ref/STARindex'

########## software
salmon_version = config['salmon_version']
stringtie_version = config['stringtie_version']
#R_version = config['R_version']
samtools_version = config['samtools_version']
gffcompare_version = config['gffcompare_version']
bedtools_version = config['bedtools_version']
scallop_version = config['scallop_version']
minimap2_version = config['minimap2_version']


rule all:
    # , expand('data/fasta/flnc/{sample}_flnc.fa', sample=pb_samplenames)
    #input:'data/gtf_info/all_RPE_loose_detdf.tsv.gz'
    #input:expand('data/stringtie_superread/{sample}.gtf', sample=[sample for sample in sr_sample_names if sr_sample_dict[sample]['origin'] == 'true_RPE'])
    input: #expand('results/all_{tissue}.merged.gtf', tissue=tissues)
         'data/combined_gtfs/all_RPE_loose.combined.gtf'
        # expand(
        #     'data/loose_set/model_results/all_RPE_dm-{dm}_wc-{wc}_kmers-{size}_dims-{dims}/model_results.csv', 
        #     dm=['0','1'],
        #     wc=['3','15','30'],
        #     size=['10', '12', '16', '20'], 
        #     dims=['100', '200', '300']),
        # expand('data/salmon_quant/all_RPE_loose/{sampleID}/quant.sf', sampleID=sr_sample_names),
        #  expand('data/gtf_info/all_RPE_{type}_target_tx.tsv', type=['strict', 'loose']),
         
        

#### build transcriptomes 
rule run_stringtie:
    input: 
        bam_path + 'STARbams/{sample}/Sorted.out.bam'
    params: 
        lib_type = lambda wildcards: get_stringtie_libflag(wildcards.sample, sr_sample_dict)
    output: 'data/stringtie/{sample}.gtf'
    shell:
        '''
        module load {stringtie_version}
        stringtie {input[0]} -o {output[0]} {params.lib_type} -l st_{wildcards.sample} -p 8 -G {ref_GTF}
        '''

def get_sr_fq_cmd(sample,fql, sample_dict):
    if sample_dict[sample]['paired']:
        l = fql + f'fastq_files/{sample}_1.fastq.gz' 
        r = fql + f'fastq_files/{sample}_2.fastq.gz'
        cmd= f' -1 {l} -2 {r}'
    else:
        up= fql + f'fastq_files/{sample}.fastq.gz'
        cmd = f'-U {up}'
    return cmd 

# rule make_super_read_gtfs:
#     input:
#         fastqs=lambda wildcards: [sr_fql + f'fastq_files/{wildcards.sample}_1.fastq.gz', sr_fql + f'fastq_files/{wildcards.sample}_2.fastq.gz'] if sr_sample_dict[wildcards.sample]['paired'] else sr_fql + f'fastq_files/{wildcards.sample}.fastq.gz'
#     params:
#         fq_cmd= lambda wildcards: get_sr_fq_cmd(wildcards.sample, sr_fql, sr_sample_dict)
#     output: 
#         super_read_bam='bams/{sample}/sr_sort.bam', 
#         gtf = 'data/stringtie_superread/{sample}.gtf'
#     shell:
#         '''
#         module load singularity
#         singularity exec --bind /data/swamyvs/,/data/OGVFB_BG/ /data/swamyvs/singularity_images/stringtie2.sif \
#             python2 /stringtie/SuperReads_RNA/create_rna_sr.py \
#                 -p 32  \
#                 {params.fq_cmd}  \
#                 -H ref/hisat2_index \
#                 -g /data/swamyvs/pacbio_testing/ \
#                 -G ref/gmap_index/gencode/  \
#                 --hisat2-cmd "/hisat2-2.2.0/hisat2" \
#                 --out-dir /data/swamyvs/pacbio_testing/bams/{wildcards.sample}/ 
#         module load {stringtie_version}
#         stringtie {output.super_read_bam} -o {output.gtf} -l stsr_{wildcards.sample} -p 32 -G {ref_GTF}
#         '''



# rule basic_filter_stringtie:
#     input:'data/stringtie/{sample}.gtf'
#     output:'data/stringtie_filtered/{sample}.gtf'
#     shell:
#         '''
#         module load {stringtie_version}
#         stringtie --merge -o {output} -F 1 -T 1 -l stf_{wildcards.sample} {input}
#         '''

# rule run_scallop:
#     input: bam_path + 'STARbams/{sample}/Sorted.out.bam'
#     params: pfx = lambda wildcards: f'data/scallop/{wildcards.sample}'
#     output: 'data/scallop/{sample}.gtf'
#     shell:
#         '''
#         module load {scallop_version}
#         scallop -i {input} -o /tmp/{wildcards.sample}
#         module load gffcompare
#         gffcompare -p sp_{wildcards.sample} -o {params.pfx} /tmp/{wildcards.sample}
#         '''



### get reference transcript quantification
# rule make_tx_fasta:
#     input:
#         gtf=lambda wildcards: build2gtf(wildcards.build)
#     output: 
#         'data/seqs/{build}_tx.fa'
#     shell:
#         '''
#         {bam_path}/gffread/gffread -w {output} -g {ref_genome}  {input.gtf}
#         '''



# rule build_salmon_index:
#     input: 
#         'data/seqs/{build}_tx.fa'
#     params:
#         tmp_dir = lambda wildcards: f'testing/salmon_tmp_{wildcards.build}/'
#     output: 
#         directory('data/salmon_indices/{build}')
#     shell:
#         '''
#         rm -rf {params.tmp_dir}
#         mkdir -p {params.tmp_dir}
#         module load {salmon_version}
#         salmon index -t {input} --tmpdir {params.tmp_dir} -i {output} 
#         '''

# rule run_salmon:
#     input: 
#         fastqs=lambda wildcards: [sr_fql+'fastq_files/{}_1.fastq.gz'.format(wildcards.sampleID),sr_fql+'fastq_files/{}_2.fastq.gz'.format(wildcards.sampleID)] if sr_sample_dict[wildcards.sampleID]['paired'] else sr_fql+'fastq_files/{}.fastq.gz'.format(wildcards.sampleID),
#         index='data/salmon_indices/{build}'
#     params: 
#         cmd=lambda wildcards: salmon_input(wildcards.sampleID,sr_sample_dict,sr_fql),\
#         outdir=lambda wildcards: f'data/salmon_quant/{wildcards.build}/{ wildcards.sampleID}/'
#     output: 
#         'data/salmon_quant/{build}/{sampleID}/quant.sf'
#     shell:
#         '''
#         id={wildcards.sampleID}
#         module load {salmon_version}
#         salmon quant -p 8 -i {input.index} -l A --gcBias --seqBias  --validateMappings {params.cmd} -o {params.outdir}
#         '''





#####long read stuff
rule clean_genome:
    input: ref_genome
    output: 'ref/gencode_genome_clean.fa'
    shell:
        '''
        python3 scripts/filterFasta.py {input} {chr_names} {output}
        '''


rule build_minimap_index:
    input: 'ref/gencode_genome_clean.fa'
    output: 'ref/minimap2_index.mmi'
    shell:
        '''
        module load {minimap2_version}
        minimap2 -d {output} {input}
        '''


rule build_talon_db:
    input: 
        gtf=ref_GTF, 
        genome = 'ref/gencode_genome_clean.fa'
    params: 
        anno_name= lambda wildcards: f'db_{wildcards.tissue}', 
        prefix= lambda wildcards: f'{wildcards.tissue}', 
        outdir= lambda wildcards: f'data/talon_db/{wildcards.tissue}'
    output: 'data/talon_db/{tissue}.db'
    shell:
        '''

        talon_initialize_database \
            --f {input.gtf} --g {input.genome} \
            --a {params.anno_name} \
            --5p 50 \
            --3p 50 \
            --idprefix {params.prefix} \
            --o {params.outdir} 

        '''

rule align_minimap2:
    input:
        fa=lambda  wildcards: lr_fql + lr_sample_dict[wildcards.sample]['filename'], 
        idx='ref/minimap2_index.mmi', 
        genome = 'ref/gencode_genome_clean.fa'
    output:
        sam='data/raw_sams/{tissue}/{sample}.sam'
    shell:
        '''
        module load {minimap2_version}
        module load {samtools_version}
        minimap2 -ax splice -uf --secondary=no -C5 {input.genome} {input.fa} | samtools sort -O sam - > {output.sam}     
        '''

rule make_spliceJN_file:
    input:
        gtf=ref_GTF, 
        genome = 'ref/gencode_genome_clean.fa'
    output:
        'data/SJ_tab.txt'
    shell:
        '''
        module load {samtools_version}
        module load {bedtools_version}
        python TranscriptClean/accessory_scripts/get_SJs_from_gtf.py --f {input.gtf} --g {input.genome} --o {output}

        '''



rule TranscriptClean: # pull from github repo
    input: 
        sam = 'data/raw_sams/{tissue}/{sample}.sam', 
        sj = 'data/SJ_tab.txt',
        genome = 'ref/gencode_genome_clean.fa'
    params:
        out_pref=lambda wildcards: f'data/clean_sams/{wildcards.tissue}/{wildcards.sample}',
        tmp_dir=lambda wildcards: f'data/tmp/{wildcards.tissue}_{wildcards.sample}/'
    output:
        sam='data/clean_sams/{tissue}/{sample}_clean.sam',
        sam_link='data/clean_sams/{tissue}/{sample}.sam'
    shell:
        '''
        rm -rf {params.tmp_dir}
        mkdir -p {params.tmp_dir}
        module load samtools
        python TranscriptClean/TranscriptClean.py --sam {input.sam} \
            --genome {input.genome}  \
            --spliceJns {input.sj}\
            --variants {variants}\
            --maxLenIndel 5 \
            --maxSJOffset 5 \
            --correctMismatches true \
            --correctIndels true \
            --correctSJs true \
            --primaryOnly \
            --outprefix {params.out_pref} \
            --tmpDir {params.tmp_dir} \
            --deleteTmp \
            --threads 64 
        ln -s {output.sam} > {output.sam_link}
        '''


'''
JK build is the genome
'''


rule run_talon:
    input: 
        sam= lambda wildcards: [f'data/clean_sams/{wildcards.tissue}/{sample}_clean.sam' for sample in lr_sample_dict.keys() if lr_sample_dict[sample]['subtissue'] == wildcards.tissue],
        db = 'data/talon_db/{tissue}.db', 
        genome = 'ref/gencode_genome_clean.fa'
    params: 
        talon_config= lambda wildcards: make_talon_config(wildcards.tissue, lr_sample_dict),
        build='gencode', 
        res_outdir=lambda wildcards: f'data/talon_results/{wildcards.tissue}/',
        prefix= lambda wildcards: f'db_{wildcards.tissue}', 
        outdir_pf=lambda wildcards: f'data/talon_results/{wildcards.tissue}/{wildcards.tissue}'
    output: 
        anno = 'data/talon_results/{tissue}/talon_read_annot.tsv', 
        qc = 'data/talon_results/{tissue}/QC.log',
        gtf='data/talon_results/{tissue}/{tissue}_talon_observedOnly.gtf',
        abundance = 'data/talon_results/{tissue}/{tissue}_talon_abundance.tsv'
    shell:
        '''
        tlncfg=/tmp/config.{wildcards.tissue}.csv
        echo "{params.talon_config}" >  $tlncfg
        talon \
            --f $tlncfg \
            --db {input.db} \
            --build {input.genome} \
            --threads 16 \
            --tmpdir {wildcards.tissue}_talon_tmp/ \
            --o {params.res_outdir}
        
        talon_create_GTF \
            --db {input.db} \
            --build {input.genome} \
            --annot {params.prefix} \
            --observed \
            --o {params.outdir_pf}
        
        talon_abundance \
            --db {input.db} \
            --build {input.genome} \
            --annot {params.prefix} \
            --o {params.outdir_pf}
        
        rm -rf {wildcards.tissue}_talon_tmp/ 
        '''


rule merge_all_gtfs:
    input: 
        stringtie_gtfs=[f'data/stringtie/{sample}.gtf'for sample in sr_sample_names if sr_sample_dict[sample]['origin'] == 'true_RPE'],
        #stringtie_sr = [f'data/stringtie_superread/{sample}.gtf'for sample in sr_sample_names if sr_sample_dict[sample]['origin'] == 'true_RPE'],
        #stringtie_filt = [f'data/stringtie_filtered/{sample}.gtf'for sample in sr_sample_names if sr_sample_dict[sample]['origin'] == 'true_RPE'],
        #scallop_gtfs = [f'data/scallop/{sample}.combined.gtf' for sample in sr_sample_names if sr_sample_dict[sample]['origin'] == 'true_RPE'],
        lr_gtfs = expand('data/talon_results/{tissue}/{tissue}_talon_observedOnly.gtf', tissue=subtissues)
    params: 
            prefix = lambda wildcards: f'data/combined_gtfs/all_RPE_{wildcards.type}',
            flag= lambda wildcards:'--strict-match' if wildcards.type == 'strict' else '' 
    output: 
        gtf = 'data/combined_gtfs/all_RPE_{type}.combined.gtf', 
        track = 'data/combined_gtfs/all_RPE_{type}.tracking'
    shell:
        '''
        module load {gffcompare_version}
        gffcompare {params.flag}  -r {ref_GTF} -o {params.prefix} {input.stringtie_gtfs} {input.scallop_gtfs} {input.stringtie_filt} {input.stringtie_sr} {input.lr_gtfs}
         '''

# moved this 
        # '''
        # module load {gffcompare_version}
        # gffcompare {params.flag}  -r {ref_GTF} -o {params.prefix} {input.stringtie_gtfs} {input.scallop_gtfs} {input.stringtie_filt} {input.stringtie_sr} {input.lr_gtfs}
        #  '''
# rule merge_all_with_filt:
#     input: 
#         all='data/combined_gtfs/all_RPE_loose.combined.gtf', 
#         filt= [f'data/stringtie_filtered/{sample}.gtf'for sample in sr_sample_names if sr_sample_dict[sample]['origin'] == 'true_RPE']
#     output: 
#         'data/combined_gtfs/all_RPE_loose_wfilter.combined.gtf'
#     shell:
#         '''
#         module load {gffcompare_version}
#         gffcompare  -r {ref_GTF} -o data/combined_gtfs/all_RPE_loose_wfilter {input.all} {input.filt}
#         '''



# rule extract_gtf_info:
#     input:
#         gtf='data/combined_gtfs/all_RPE_{type}.combined.gtf',
#         track_file='data/combined_gtfs/all_RPE_{type}.tracking',
#         abundance_file='data/talon_results/RPE_Fetal.Tissue/RPE_Fetal.Tissue_talon_abundance.tsv'
#     params:
#         out_prefix='data/gtf_info/all_RPE_{type}'
#     output:
#         convtab = 'data/gtf_info/all_RPE_{type}_convtab.tsv.gz',
#         det_df = 'data/gtf_info/all_RPE_{type}_detdf.tsv.gz',
#         target_tx = 'data/gtf_info/all_RPE_{type}_target_tx.tsv',
#         validaion_tx = 'data/gtf_info/all_RPE_{type}_validation_tx.tsv'
#     shell:
#         '''
#         module load R
#         Rscript scripts/clean_track_file.R \
#             --workingDir {working_dir} \
#             --rawTrackFile {input.track_file} \
#             --refGtfFile {ref_GTF} \
#             --gtfFile {input.gtf} \
#             --uniprotFile {uniprot_file} \
#             --talonAbFile {input.abundance_file} \
#             --outPfx {params.out_prefix}

#         '''



# rule kmerize_transcripts:
#     input: 
#         fasta='data/seqs/all_RPE_{type}_tx.fa',
#         target_tx = 'data/gtf_info/all_RPE_{type}_target_tx.tsv',
#         validaion_tx = 'data/gtf_info/all_RPE_{type}_validation_tx.tsv'
#     params:
#         out_prefix = lambda wildcards: f'data/{wildcards.type}_set/raw_model_data/all_RPE_kmers_{wildcards.size}'
#     output:  
#         full_lsf = 'data/{type}_set/raw_model_data/all_RPE_kmers_{size}_train.lsf', 
#         full_txids = 'data/{type}_set/raw_model_data/all_RPE_kmers_{size}_train_txids.txt',
#         val_lsf = 'data/{type}_set/raw_model_data/all_RPE_kmers_{size}_val.lsf', 
#         val_txids = 'data/{type}_set/raw_model_data/all_RPE_kmers_{size}_val_txids.txt'
#     shell:
#         '''
#         python3  scripts/kmerize_fasta_low_mem.py \
#             --infasta {input.fasta} \
#             --kmerSize {wildcards.size} \
#             --trainTx {input.target_tx} \
#             --valTx {input.validaion_tx} \
#             --outPfx {params.out_prefix}
           
#         '''


# '''
# 03/11/20 changed so we are training on only the target tx, maybe this will help thigns 
# q
# '''

# rule train_doc2vec_and_infer_vectors:
#     input: 
#         full_lsf = 'data/{type}_set/raw_model_data/all_RPE_kmers_{size}_train.lsf', 
#         full_txids = 'data/{type}_set/raw_model_data/all_RPE_kmers_{size}_train_txids.txt',
#         val_lsf = 'data/{type}_set/raw_model_data/all_RPE_kmers_{size}_val.lsf', 
#         val_txids = 'data/{type}_set/raw_model_data/all_RPE_kmers_{size}_val_txids.txt'
#     params:
#         model= lambda wildcards:  f'data/{wildcards.type}_set/embedding_models/all_RPE_doc2vec_ep-15_dm-{wildcards.dm}_wc-{wildcards.wc}_kmers-{wildcards.size}_dims-{wildcards.dims}.pymodel'
#     output: 
#         out_matrix = 'data/{type}_set/embedded_model_data/all_RPE_dm-{dm}_wc-{wc}_kmers-{size}_dims-{dims}.csv.gz',
#         val_matrix = 'data/{type}_set/embedded_model_data/all_RPE_validation-tx_dm-{dm}_wc-{wc}_kmers-{size}_dims-{dims}.csv.gz'
        
#     shell:
#         '''
#         python3 scripts/train_doc2vec.py \
#             --corpusFile {input.full_lsf} \
#             --corpusTxIDs {input.full_txids} \
#             --valTxLsf {input.val_lsf} \
#             --valTxIDs {input.val_txids}\
#             --edim {wildcards.dims} \
#             --wc {wildcards.wc} \
#             --dm {wildcards.dm} \
#             --trainedModel {params.model} \
#             --loadModel \
#             --outTrainMatrix {output.out_matrix} \
#             --outValMatrix {output.val_matrix}
#         '''


# rule run_experiment:
#     input:
#         out_matrix = 'data/{type}_set/embedded_model_data/all_RPE_dm-{dm}_wc-{wc}_kmers-{size}_dims-{dims}.csv.gz'
#     params:
#         out_dir = lambda wildcards: f'data/{wildcards.type}_set/model_results/all_RPE_dm-{wildcards.dm}_wc-{wildcards.wc}_kmers-{wildcards.size}_dims-{wildcards.dims}/'
#     output:
#         out_csv = 'data/{type}_set/model_results/all_RPE_dm-{dm}_wc-{wc}_kmers-{size}_dims-{dims}/model_results.csv'
#     shell:
#         '''
#         python3 scripts/run_models.py \
#             --workingDir {working_dir} \
#             --inputFile {input.out_matrix} \
#             --labFile data/loose_set/all_RPE_loose_labels.csv  \
#             --inputType skl \
#             --model rf,rf_weighted \
#             --nproc 32 \
#             --outputDir {params.out_dir}
    
#         '''



