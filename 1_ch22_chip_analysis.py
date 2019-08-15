#!/usr/bin/python


'''
The MIT License (MIT)

Copyright (c) 2018 Charles Lin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
'''


#Main method run script for processing of slam seq analysis from Muhar et al., 2018




#==========================================================================
#=============================DEPENDENCIES=================================
#==========================================================================


import sys, os
# Get the script's full local path
whereAmI = os.path.dirname(os.path.realpath(__file__))

pipeline_dir = '/storage/cylin/bin/pipeline/'
py27_path = '/storage/cylin/anaconda3/envs/py27_anaconda/bin/python2'

sys.path.append(whereAmI)
sys.path.append(pipeline_dir)

import pipeline_dfci
import utils
import string
import numpy
import os
import re
from collections import defaultdict
import subprocess
#==========================================================================
#============================PARAMETERS====================================
#==========================================================================



projectName = 'chordoma_ch22_chip'
genome ='hg19'
annotFile = '%s/annotation/%s_refseq.ucsc' % (pipeline_dir,genome)

#project folders
projectFolder = '/storage/cylin/grail/projects/%s' % (projectName) #PATH TO YOUR PROJECT FOLDER


projectFolder = utils.formatFolder(projectFolder,True)
#standard folder names
gffFolder ='%sgff/' % (projectFolder)
macsFolder = '%smacsFolder/' % (projectFolder)
macsEnrichedFolder = '%smacsEnriched/' % (projectFolder)
mappedEnrichedFolder = '%smappedEnriched/' % (projectFolder)
mappedFolder = '%smappedFolder/' % (projectFolder)
wiggleFolder = '%swiggles/' % (projectFolder)
metaFolder = '%smeta/' % (projectFolder)
metaRoseFolder = '%smeta_rose/' % (projectFolder)
roseFolder = '%srose/' % (projectFolder)
fastaFolder = '%sfasta/' % (projectFolder)
bedFolder = '%sbed/' % (projectFolder)
figuresFolder = '%sfigures/' % (projectFolder)
geneListFolder = '%sgeneListFolder/' % (projectFolder)
bedFolder = '%sbeds/' % (projectFolder)
signalFolder = '%ssignalTables/' % (projectFolder)
tableFolder = '%stables/' % (projectFolder)

#mask Files
maskFile = '%smasks/hg19_encode_blacklist.bed' % (projectFolder)

#genomeDirectory #select your genome
genomeDirectory = '/storage/cylin/grail/genomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/'
#making folders
folderList = [gffFolder,macsFolder,macsEnrichedFolder,mappedEnrichedFolder,mappedFolder,wiggleFolder,metaFolder,metaRoseFolder,roseFolder,fastaFolder,figuresFolder,geneListFolder,bedFolder,signalFolder,tableFolder]

for folder in folderList:
    pipeline_dfci.formatFolder(folder,True)



#==========================================================================
#============================LIST OF DATAFILES=============================
#==========================================================================

#this project will utilize multiple datatables
#data tables are organized largely by type/system
#some data tables overlap for ease of analysis

#ChIP-Seq
chip_data_file = '%sdata_tables/CH22_CHIP_TABLE.txt' % (projectFolder)




#==========================================================================
#===========================MAIN METHOD====================================
#==========================================================================


def main():


    print('main analysis for project %s' % (projectName))

    print('changing directory to project folder')
    os.chdir(projectFolder)

    print('\n\n')
    print('#======================================================================')
    print('#======================I. LOADING DATA ANNOTATION======================')
    print('#======================================================================')
    print('\n\n')

    #This section sanity checks each data table and makes sure both bam and .bai files are accessible

    #for data file
    pipeline_dfci.summary(chip_data_file)


    #use macs1.4 to make wiggles
    macs14Folder = utils.formatFolder('%smacs14/' % (projectFolder),True)
    macs14EnrichedFolder = utils.formatFolder('%smacs14Enriched/' % (projectFolder),True)
    pipeline_dfci.run_macs(chip_data_file,projectFolder,macs14Folder,macs14EnrichedFolder,wiggleFolder,useBackground=True)

    print('\n\n')
    print('#======================================================================')
    print('#=======================II. DEFINING ENHANCERS=========================')
    print('#======================================================================')
    print('\n\n')

    #running ROSE2_meta on the chordoma k27ac

    chip_data_dict = pipeline_dfci.loadDataTable(chip_data_file)
    # for name in chip_data_dict.keys():
    #     print(name)
    #     print(chip_data_dict[name]['enrichedMacs'])


    # #run rose2 wrapper for both
    # bashFileName,region_map_path,names_list = define_enhancer_landscape(projectFolder,pipeline_dir,chip_data_file,analysis_name = 'CH22_H3K27AC')
    # print(bashFileName,region_map_path,names_list)

    # #runs only if no output detected                                                                      
    # if not utils.checkOutput(enhancer_region_map_path,0,0):                                               
    #     print(enhancer_bashFileName)                                                                      
    #     os.system('bash %s' % (enhancer_bashFileName))   


    # #=========
    # #=========
    # #=========

    # #sanity check debug

    # #run rose2 meta for one dataset as a test w/ stitch at 500 just for a control
    # bashFileName,region_map_path,names_list = define_enhancer_landscape(projectFolder,pipeline_dir,chip_data_file,analysis_name = 'CH22_H3K27AC_1_TEST',names_list = ['CH22_H3K27AC_1'],stitch = '500')
    # print(bashFileName,region_map_path,names_list)

    # #running regular rose2
    # rose2_parent_folder = utils.formatFolder('%srose2' % (projectFolder),True)
    # rose2_bash = pipeline_dfci.callRose2(chip_data_file,macsEnrichedFolder,rose2_parent_folder,namesList=['CH22_H3K27AC_1'],extraMap = [],inputFile='',tss=2500,stitch='500',bashFileName ='',mask=maskFile,useBackground=True,py27_path =py27_path)

    # print(rose2_bash)

    # #ok, ROSE2 META and ROSE2 still produce same result when run on single dataset
    # #whew

    # #=========
    # #=========
    # #=========

    print('\n\n')
    print('#======================================================================')
    print('#=======================III. DEFINE T LANDSCAPE========================')
    print('#======================================================================')
    print('\n\n')


    # #since this is an HA chip we need to remove HA background
    
    # #get the pos regions (T)
    # data_dict = pipeline_dfci.loadDataTable(chip_data_file)
    # t_list = ['%s%s' % (macsEnrichedFolder,data_dict[name]['enrichedMacs']) for name in data_dict.keys() if name.count('dTag_T') == 1 and name.count('WCE') == 0]
    # print(t_list)

    # #get the negative_control list from IRF2 project
    # ha_ctl_list = ['%sHK_CTL_HA_1_peaks.bed' % (macsEnrichedFolder),'%sHK_CTL_HA_2_peaks.bed' % (macsEnrichedFolder)]
    
    # t_bed_path_intersect = '%sCH22_T_INTERSECT.bed' % (bedFolder)
    # #merge_regions(pos_list = t_list,neg_list = ha_ctl_list,analysis_name = 'CH22_T_INTERSECT',output_path=t_bed_path_intersect,merge_type = 'INTERSECT')

    # t_bed_path_union = '%sCH22_T_UNION.bed' % (bedFolder)
    # #merge_regions(pos_list = t_list,neg_list = ha_ctl_list,analysis_name = 'CH22_T_UNION',output_path=t_bed_path_union,merge_type = 'UNION')

    # #for k27ac

    # h3k27ac_list = ['%s%s' % (macsEnrichedFolder,data_dict[name]['enrichedMacs']) for name in data_dict.keys() if name.count('H3K27AC') == 1 and name.count('WCE') == 0]
    # h3k27ac_bed_path_union = '%sCH22_H3K27AC_UNION.bed' % (bedFolder)
    # merge_regions(pos_list = h3k27ac_list,neg_list = [],analysis_name = 'CH22_H3K27AC_UNION',output_path=h3k27ac_bed_path_union,merge_type = 'UNION')


    print('\n\n')
    print('#======================================================================')
    print('#=========================IV. T MOTIF FINDING==========================')
    print('#======================================================================')
    print('\n\n')

    #use the T union and then calculate signal
    
    data_dict = pipeline_dfci.loadDataTable(chip_data_file)
    map_list = [name for name  in data_dict.keys() if name.count('dTag_T') == 1]
    print(map_list)
    gffList = ['%sCH22_T_UNION.bed' % (bedFolder)]
    #signal_table_list = pipeline_dfci.map_regions(chip_data_file,gffList,mappedFolder,signalFolder,names_list=map_list,medianNorm=False,output='',extendReadsTo=200)

    #column order =  CH22_dTag_T_MUT_HA      CH22_dTag_T_MUT_WCE     CH22_dTag_T_WT_HA       CH22_dTag_T_WT_WCE

    signal_table_path= '%sCH22_T_UNION_CH22_CHIP_TABLE_SIGNAL.txt' % (signalFolder)

    top =1000
    fasta_path = make_T_top_regions(signal_table_path,top)

    #now run meme
    analysis_name = 'HG19_CH22_T_UNION_TOP'    
    meme_bash_path = wrap_meme(analysis_name)

    print('\n\n')
    print('#======================================================================')
    print('#===================IV. MAKE HEATMAPS OF T LANDSCAPE===================')
    print('#======================================================================')
    print('\n\n')

    #use the union for T
    t_union_path = '%sbeds/CH22_T_UNION.bed' % (projectFolder)


    print('\n\n')
    print('#======================================================================')
    print('#=======================IV. DEFINE ACTIVE GENES========================')
    print('#======================================================================')
    print('\n\n')

    # #make the relevant gffs
    # pipeline_dfci.makeGeneGFFs(annotFile,gffFolder,species='HG19')

    # # #Making a list of all active genes
    # tss_gff = utils.parseTable('%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder),'\t')
    # start_dict = utils.makeStartDict(annotFile)

    # all_gene_table = []
    # ticker = 1
    # for line in tss_gff:
    #     new_line = [ticker,line[1],start_dict[line[1]]['name']]
    #     all_gene_table.append(new_line)
    #     ticker+=1
        
    # utils.unParseTable(all_gene_table,'%sHG19_UCSC_REFSEQ_ALL.txt' % (geneListFolder),'\t')
    # sys.exit()

    # setName = 'CH22_H3K27AC'
    # cellTypeList = ['CH22']
    # map_list = ['CH22_H3K27AC_1','CH22_H3K27AC_2']
    # gffList = ['%sHG19_TSS_ALL_-1000_+1000.gff' % (gffFolder)]
    # pipeline_dfci.mapEnrichedToGFF(chip_data_file,setName,gffList,cellTypeList,macsEnrichedFolder,mappedEnrichedFolder,macs=True,namesList=map_list,useBackground=True)


    # setList = [['CH22_H3K27AC_1'],['CH22_H3K27AC_2']] #bound by either
    # output = '%sHG19_CHORDOMA_CH22_H3K27AC_ACTIVE.txt' % (geneListFolder)
    # mappedEnrichedFile = '%sHG19_TSS_ALL_-1000_+1000/HG19_TSS_ALL_-1000_+1000_CH22_H3K27AC.txt' % (mappedEnrichedFolder)
    # pipeline_dfci.makeGFFListFile(mappedEnrichedFile,setList,output,annotFile)

    print('\n\n')
    print('#======================================================================')
    print('#====================V. RUNNING ENHANCER PROMOTER======================')
    print('#======================================================================')
    print('\n\n')

    # data_dict = pipeline_dfci.loadDataTable(chip_data_file)

    # #need to run enhancer promoter code on both k27ac and T


    # #for T at active genes
    # activity_path = '%sHG19_CHORDOMA_CH22_H3K27AC_ACTIVE.txt' % (geneListFolder)
    # input_path = '%sCH22_T_UNION.bed' % (bedFolder)
    # analysis_name = 'CH22_T_UNION'
    # t_list = [name for name in data_dict.keys() if name.count('dTag_T') == 1 and name.count('WCE') == 0]

    # #wrap_enhancer_promoter(chip_data_file,input_path,activity_path,analysis_name,names_list = t_list,useBackground=True)
    
    # #for T at all genes
    # activity_path = '%sHG19_UCSC_REFSEQ_ALL.txt' % (geneListFolder)
    # input_path = '%sCH22_T_UNION.bed' % (bedFolder)
    # analysis_name = 'CH22_T_UNION_ALL_GENES'
    # t_list = [name for name in data_dict.keys() if name.count('dTag_T') == 1 and name.count('WCE') == 0]

    # #wrap_enhancer_promoter(chip_data_file,input_path,activity_path,analysis_name,names_list = t_list,useBackground=True)



    # #for H3K27AC at active genes
    # # activity_path = '%sHG19_CHORDOMA_CH22_H3K27AC_ACTIVE.txt' % (geneListFolder)

    # # input_path = '%sCH22_H3K27AC_UNION.bed' % (bedFolder)
    # # analysis_name = 'CH22_H3K27C_UNION'
    # # h3k27ac_list = [name for name in data_dict.keys() if name.count('H3K27AC') == 1 and name.count('WCE') == 0]

    # # wrap_enhancer_promoter(chip_data_file,input_path,activity_path,analysis_name,names_list = h3k27ac_list,useBackground=True)


    # #for H3K27AC at all genes
    # activity_path = '%sHG19_UCSC_REFSEQ_ALL.txt' % (geneListFolder)

    # input_path = '%sCH22_H3K27AC_UNION.bed' % (bedFolder)
    # analysis_name = 'CH22_H3K27C_UNION_ALL_GENES'
    # h3k27ac_list = [name for name in data_dict.keys() if name.count('H3K27AC') == 1 and name.count('WCE') == 0]

    # wrap_enhancer_promoter(chip_data_file,input_path,activity_path,analysis_name,names_list = h3k27ac_list,useBackground=True)


    print('\n\n')
    print('#======================================================================')
    print('#=================VI. LINKING CHROMATIN TO EXPRESSION==================')
    print('#======================================================================')
    print('\n\n')


    #a gene counts if it is expressed above cut in at least one sample
    #may need to collapse NM IDs per genes

    #first check that expression table doesn't have duplicates

    def merge_rna(exp_path,exp_cutoff = 1,output_path =''):

        '''
        just a wrapper for combining expression data w/ gene level h3k27ac and T data
        '''
        

        
        exp_table = utils.parseTable(exp_path,'\t')
        
        exp_dict= defaultdict(list)
        for line in exp_table[1:]:
            #here's where we can filter for an expression cutoff
            exp_line = [float(x) for x in  line[1:]]
            if max(exp_line) > exp_cutoff:
                exp_dict[line[0]] = exp_line
        
        #now figure out genes w/ T binding
        t_gene_path = '%senhancerPromoter/CH22_T_UNION_ALL_GENES/CH22_T_UNION_ALL_GENES_GENE_TABLE.txt' % (projectFolder)
        t_table  = utils.parseTable(t_gene_path,'\t')

        t_dict = defaultdict(list)
        
        for line in t_table[1:]:
            t_dict[line[0]] = [float(x) for x in line[1:]]

        #now figure out genes w/ H3K27AC binding
        h3k27ac_gene_path = '%senhancerPromoter/CH22_H3K27C_UNION_ALL_GENES/CH22_H3K27C_UNION_ALL_GENES_GENE_TABLE.txt' % (projectFolder)
        h3k27ac_table  = utils.parseTable(h3k27ac_gene_path,'\t')

        h3k27ac_dict = defaultdict(list)
        
        for line in h3k27ac_table[1:]:
            h3k27ac_dict[line[0]] = [float(x) for x in line[1:]]


        #now set up the output
        gene_table = []
        gene_table_header = ['GENE','T_PROMOTER','T_DISTAL','H3K27AC_PROMOTER','H3K27AC_DISTAL'] + exp_table[0]
        gene_table.append(gene_table_header)

        #anchor analysis on genes w/ detectable expr
        exp_gene_list = exp_dict.keys()
        exp_gene_list.sort()
        for gene in  exp_gene_list:
            
            if gene in t_dict:
                t_line = t_dict[gene]
            else:
                t_line = [0.0,0.0]

            if gene in h3k27ac_dict:
                h3k27ac_line = h3k27ac_dict[gene]
            else:
                h3k27ac_line = [0.0,0.0]

            new_line = [gene] + t_line + h3k27ac_line + exp_dict[gene]
            
            gene_table.append(new_line)


        utils.unParseTable(gene_table,output_path,'\t')



    #for norm data
    rna_project_folder = '/storage/cylin/grail/projects/chordoma_ch22_rna/'

    #this table is just in alphabetical order
    exp_path = '%s190612_rna_seq/cuffnorm_output/cuffnorm_all_fpkm_exprs_norm.txt' % (rna_project_folder)

            
    exp_cutoff = 1
    output_path = '%stables/HG19_CHORDOMA_CH22_GENE_TABLE_NORM.txt' % (projectFolder)
    #merge_rna(exp_path,exp_cutoff,output_path)

    #for raw data
    exp_path = '%s190612_rna_seq/cuffnorm_output/cuffnorm_all_fpkm_exprs_raw.txt' % (rna_project_folder)
    output_path = '%stables/HG19_CHORDOMA_CH22_GENE_TABLE_RAW.txt' % (projectFolder)
    #merge_rna(exp_path,exp_cutoff,output_path)


#==========================================================================
#==========================SCRIPT FUNCTIONS================================
#==========================================================================


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~DEFINING H3K27AC ENHANCER LANDSCAPE~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def define_enhancer_landscape(projectFolder,pipeline_dir,data_file,analysis_name = '',names_list = [],stitch='',use_background = True):

    '''
    defines the NB enhancer baseline using H3K27ac chips 
    enhancers defined using auto optimized stitching of nearby regions
    w/ a 2.5kb tss exclusion
    uses the meta rose code and writes out a .sh file for reproducibility
    '''

    #For H3K27AC
    #with TSS exclusion and auto stitching

    dataDict = pipeline_dfci.loadDataTable(data_file)

    if  analysis_name == '':
        #get the name from the data_file
        analysis_name = data_file.split('/')[-1].split('.')[0]
    
    print("RUNNING ANALYSIS: '%s'" % (analysis_name))


    if names_list == []:
        names_list = [name for name in dataDict.keys() if name.upper().count('H3K27AC') == 1]
    print('FOR H3K27AC DATASETS USING:')
    print(names_list)


    bamFileList = [dataDict[name]['bam'] for name in names_list]
    bamString = string.join(bamFileList,',')

    controlBams = [dataDict[name]['background'] for name in names_list]
    controlFileList = [dataDict[name]['bam'] for name in controlBams]
    controlBamString = string.join(controlFileList,',')

    bedFileList = [macsEnrichedFolder + dataDict[name]['enrichedMacs'] for name in names_list]
    bedString = string.join(bedFileList,',')

    roseFolder = '%smeta_rose/' % (projectFolder)
    roseFolder = utils.formatFolder(roseFolder,True)

    outputFolder = '%s%s/' % (roseFolder,analysis_name)
    bashFileName = '%s%s_meta_rose.sh' % (roseFolder,analysis_name)

    bashFile = open(bashFileName,'w')
    bashFile.write('#!/usr/bin/bash\n\n')
    bashFile.write('cd %s\n' % (pipeline_dir))

    metaRoseCmd = '%s %sROSE2_META.py -g hg19 -i %s -r %s -c %s -o %s -n %s -t 2500 --mask %s' % (py27_path,pipeline_dir,bedString,bamString,controlBamString,outputFolder,analysis_name,maskFile)
    if stitch != '':
        metaRoseCmd += ' -s %s' % (stitch)

    bashFile.write(metaRoseCmd + '\n')
    bashFile.close()


    #getting the region_map_path as a way to know if it's done
    region_map_path = '%s%s/%s_AllEnhancers.table.txt' % (roseFolder,analysis_name,analysis_name)
    return bashFileName,region_map_path,names_list



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~DEFINING CONSENSUS CHIP REGIONS~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    


def merge_regions(pos_list = [],neg_list = [],analysis_name = '',output_path='',merge_type = 'INTERSECT'):

    '''
    merges positive regions and filters out negative regions
    takes in as arguments the actual paths  of the regions
    '''

    #for positive regions
    pos_loci = []
    peak_dict = {}
    for peak_bed_path in pos_list:
        #get a name from the peak bed path
        peak_name = peak_bed_path.split('/')[-1].split('.')[0]
        peak_collection  = utils.importBoundRegion(peak_bed_path,peak_name)
        pos_loci += peak_collection.getLoci()
        peak_dict[peak_bed_path] = peak_collection

    pos_collection = utils.LocusCollection(pos_loci)
    pos_stitched_collection = pos_collection.stitchCollection()


    #for negative regions
    neg_loci = []
    for peak_bed_path in neg_list:
        #get a name from the peak bed path
        peak_name = peak_bed_path.split('/')[-1].split('.')[0]
        peak_collection  = utils.importBoundRegion(peak_bed_path,peak_name)
        neg_loci += peak_collection.getLoci()

    neg_collection = utils.LocusCollection(neg_loci)
    neg_stitched_collection = neg_collection.stitchCollection()

    filtered_loci = []
    if merge_type == 'INTERSECT':

        for locus in pos_stitched_collection.getLoci():
            #for each pos region first make sure it overlaps original loci
            overlap_test_list= []
            for peak_bed_path in pos_list:        
                overlap_test_list.append(peak_dict[peak_bed_path].getOverlap(locus))

            if min(overlap_test_list) != [] and neg_stitched_collection.getOverlap(locus) == []:
                filtered_loci.append(locus)

    if merge_type == 'UNION':
        #only filter for negative regions
 
        for locus in pos_stitched_collection.getLoci():

            if neg_stitched_collection.getOverlap(locus) == []:
                filtered_loci.append(locus)
       

    print('FOR %s, THERE ARE %s ORIGINAL PEAKS' % (analysis_name,len(pos_stitched_collection)))
    print('%s PEAKS REMAIN AFTER FILTERING FOR NEGATIVE REGIONS AND MERGING BY %s' % (len(filtered_loci),merge_type))

    #convert the filtered regions into a bed
    filtered_collection = utils.LocusCollection(filtered_loci)
    filtered_bed = utils.locusCollectionToBed(filtered_collection)

    if output_path == '':
        output_path = './%s.bed' % (analysis_name)

    print('WRITING BED OUTPUT TO %s' % (output_path))
    utils.unParseTable(filtered_bed,output_path,'\t')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~DEFINING TOP T REGIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def make_T_top_regions(signal_table_path,top=1000):

    '''
    makes a gff and fasta of the top N T regions based off the signal table
    '''

    signal_table = utils.parseTable(signal_table_path,'\t')

    signal_dict = defaultdict(float)
    for line in signal_table[1:]:

        signal = (max(float(line[2]) - float(line[3]),0) + max(float(line[4]) - float(line[5]),0))/2
        signal_dict[line[1]] = signal


    signal_vector = [signal_dict[line[1]] for line in signal_table[1:]]
    signal_order = utils.order(signal_vector,decreasing=True)

    t_top_gff_path = '%sCH22_T_UNION_TOP_%s_-0_+0.gff' % (gffFolder,str(top))
    print(t_top_gff_path)
    t_top_gff = []
    for i in range(top):
        signal_row = signal_order[i] + 1
        line = signal_table[signal_row]
        region_id = line[1]
        chrom = region_id.split('(')[0]
        coords = region_id.split(':')[-1].split('-')
        gff_line = [chrom,region_id,'',coords[0],coords[1],'','.','',region_id]
        t_top_gff.append(gff_line)

    utils.unParseTable(t_top_gff,t_top_gff_path,'\t')

    t_top_fasta = utils.gffToFasta('HG19',genomeDirectory,t_top_gff)
    t_top_fasta_path = '%sHG19_CH22_T_UNION_TOP_%s_-0_+0.fasta' % (fastaFolder,top)
    utils.unParseTable(t_top_fasta,t_top_fasta_path,'')
    
    return t_top_fasta_path


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~WRAPPING MEME~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def wrap_meme(analysis_name):

    '''
    wrapper to run meme-chip w/ a pwm
    '''
    meme_folder = utils.formatFolder('%smeme/' % (projectFolder),True)

    output_folder = utils.formatFolder('%s%s' % (meme_folder,analysis_name),True)

    meme_bash_path = '%s%s_%s_meme.sh' % (meme_folder,analysis_name,top)
    meme_path = '/storage/cylin/bin/meme/bin/meme-chip'
    pwm_path = '/storage/cylin/bin/pipeline/crc/annotation/VertebratePWMs.txt'

    meme_bash = open(meme_bash_path,'w')
    meme_bash.write('#!/usr/bin/bash\n')
    meme_bash.write('#SBATCH -n 32\n')
    meme_bash.write('#SBATCH -p short\n')

    meme_cmd = '%s -meme-nmotifs 5 -spamo-skip -oc %s -db %s %s' % (meme_path,output_folder,pwm_path,fasta_path)

    meme_bash.write(meme_cmd)

    meme_bash.close()

    return meme_bash_path




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~RUNNING ENHANCER PROMOTER ANALYSIS~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def wrap_enhancer_promoter(data_file,input_path,activity_path,analysis_name,names_list = [],useBackground=True):

    '''
    runs enhancer promoter on everybody with the conserved regions and union of active genes
    '''
    
    #hard coded paths
    tads_path ='%shESC_domains_hg19.bed' %(bedFolder)

    #setting the output folder
    ep_folder = utils.formatFolder('%senhancerPromoter/' % (projectFolder),True)
    

    dataDict = pipeline_dfci.loadDataTable(data_file)
    if len(names_list) == 0:
        names_list = [name for name in dataDict.keys() if name.count('WCE') == 0] #basically everybody except WCE
        names_list.sort()

    print('RUNNING ENHANCER PROMOTER ANALYSIS ON:')
    print(names_list)

    bams_list = [dataDict[name]['bam'] for name in names_list]
    bams_string = ' '.join(bams_list)
    
    background_names = [dataDict[name]['background'] for name in names_list]
    background_list = [dataDict[background_name]['bam'] for background_name in background_names]
    background_string = ' '.join(background_list)

    ep_bash_path = '%s%s_enhancer_promoter.sh' % (ep_folder,analysis_name)
    ep_bash = open(ep_bash_path,'w')

    ep_bash.write('#!/usr/bin/bash\n\n\n')
    
    ep_bash.write('#enhancer promoter analysis for %s\n\n' % (analysis_name))

    if useBackground:
        python_cmd = '%s %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 2000\n\n' % (py27_path,pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)

        ep_bash.write(python_cmd)

        python_cmd = '%s %senhancerPromoter.py -b %s -c %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 5000\n\n' % (py27_path,pipeline_dir,bams_string,background_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)

        ep_bash.write(python_cmd)
    else:
        python_cmd = '%s %senhancerPromoter.py -b %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 2000\n\n' % (py27_path,pipeline_dir,bams_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)

        ep_bash.write(python_cmd)

        python_cmd = '%s %senhancerPromoter.py -b %s -g %s -i %s -o %s -a %s --name %s --tads %s --top 5000\n\n' % (py27_path,pipeline_dir,bams_string,genome.upper(),input_path,ep_folder,activity_path,analysis_name,tads_path)

        ep_bash.write(python_cmd)

    ep_bash.close()
    
    return(ep_bash_path)


#==========================================================================
#==================================THE END=================================
#==========================================================================

    
if __name__=="__main__":
    main()
