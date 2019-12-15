# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 17:24:26 2019
Python 3.6
@author: woutv
Requirements: this file must be run in a folder that also contains all output folders from SNPeffects
"""
import pandas as pd
from pathlib import Path
import os
import requests
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import re
import csv


# SET THIS PARAMETER TO THE ORGANISM THAT CORRESPONDS WITH YOUR DATA 
organism = 'human'

print("////////////// Absolute Mutation Position Filter \\\\\\\\\\\\\\\\\\\\\\\\\\\\")


# ----------------------------------- Data preprocessing -----------------------------------------------


def find_files():  # goes through all folders in the same dir and returns a list of finalreports
    report = {}
    dir_name = os.path.dirname(os.path.realpath(__file__))
    print("Auto detecting files in " + dir_name)
    path_list = Path(dir_name).glob('**/finalreport.tab')
    for path in path_list:
        folder_name = os.path.split(os.path.split(str(path))[0])[1]
        report[folder_name] = pd.read_csv(str(path), sep='\t')
    return report, dir_name


#grouping duplicates together
def group_duplicates(df_in, folder):
    checklist = {}
    groups = []  # this is a list of lists
    df_in.columns = \
        df_in.columns.str.strip().str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '')
    for row in df_in.iterrows():
        duplicationkey = repr([row[1]['pdb_file'], row[1]['mutationstring']])
        if checklist.get(duplicationkey) is None:
            checklist[duplicationkey] = [row]
        else:
            checklist[duplicationkey].append(row)
    for duplicationkey in checklist:
        groups.append(checklist[duplicationkey])
    return groups


def find_refseq(groups):
    print('Obtaining reference sequences from Uniprot:')
    refDict = {}
    uniprotIDlist = []
    queryErrors = []
    i = 0
    printProgressBar(0, len(groups), prefix='', suffix='Complete', length=50)
    for group in groups:
        
        geneName = group[0][1]['gene_name']
        # entering the geneName into the reference dictionary with a default '' as sequence.
        # if the query is successful this will be replaced with a reference sequence.
        refDict[geneName] = ''
        
        #there may be several gene_names concatenated in the geneNames object therefore we split on the separator '-'
        for name in geneName.split('-'):
            # we only want to query agani if the previous runs did not find anything.
            if refDict[geneName] == '':
                #query parameters
                BASE = 'http://www.uniprot.org'
                KB_ENDPOINT = '/uniprot/'
                payload = {'query': 'gene_exact:' + name + ' ' + organism, 'sort': 'score', 'format': 'fasta'}
                
                #query to uniprot
                results = requests.get(BASE + KB_ENDPOINT, params=payload)
                if results.ok:
                    # if the query returned something we extract the uniprotID and the reference sequence from it 
                    if results.text != '':
                        results = results.text
                        resultsLIST = results.split('>')
                        # adding the uniprot ID to the uniprotIDlist
                        uniprotID = str(re.sub('[A-Z]+$', '', re.sub('\n', '', resultsLIST[1])).split('|')[1])  # get the uniprot ID
                        for var in group:
                            uniprotIDlist.append(uniprotID)
                        # adding the reference sequence for a given geneName to the dictionary
                        # note eventhough geneName is not a unique key in the original dataset,
                        # we can still use it here since the query result will be the same given a geneName 
                        res = re.sub(".*\d", '', resultsLIST[1]).split('\n')
                        resSTR = ''
                        for line in res:
                            resSTR += line                  
                        refDict[geneName] = resSTR
                
                # if the results of the query were not ".ok" for what ever reason we want to report it
                else:
                    queryErrors.append('query for ' + geneName + ' name: ' + name + ' failed.\n' + results.status_code)
        
        # if after the queries for all names in geneNames the reference dictionary still
        # has no reference sequence for this geneName we leave it empty and will give uniprotIDlist 
        # an NA for each variant in this group and register it in the query error list.
        if refDict[geneName] == '':
            for var in group:
                uniprotIDlist.append('NA')
            queryErrors.append('query for ' + geneName + ' returned no results')
        
        #update the progress bar
        printProgressBar(i + 1, len(groups), prefix='', suffix='Complete', length=50)
        i += 1
    
    return [refDict, uniprotIDlist, queryErrors]


# find mutated sequence of each variant
def find_mutseq(groups, datasetname):
    recordDict = {}
    groupVarDict = {}
    seq_records = SeqIO.parse(headfolder + '/' + datasetname + '/' + 'mutated_sequences.fa', 'fasta')
    for record in seq_records:
        recordDict[str(record.id)] = str(record.seq)
    for group in groups:
        groupID = repr([group[0][1]['pdb_file'], group[0][1]['mutationstring']])
        geneName = group[0][1]['gene_name']
        varDict = {}
        for var in group:
            variant = var[1]['variant']
            varDict[variant] = recordDict[variant]  # {var_XX : MHTOHADZRACWGHWRY, ...}
            groupVarDict[groupID] = [geneName, varDict]
    return groupVarDict

#alignment of variants with uniprot reference sequence.
#incase an out of memory exception happens because the sequence is too long, the 
def align_groups(refDict, groupVarDict, x):
    print('Aligning variant sequences to their uniprot reference sequence:')
    alignmentfile = headfolder + '/' + x + '/' + 'Alignments' + x + '.txt'
    out = open(alignmentfile, 'a')
    group_bestmatch = []
    scoreList = []
    alignErrors = []
    i = 0
    printProgressBar(0, len(groupVarDict), prefix='', suffix='Complete', length=50)
    for groupID in groupVarDict:
        geneName = groupVarDict[groupID][0]
        varDict = groupVarDict[groupID][1]
        group_alignments = []
        #reasonably low base alignment score to surpass by the best alignment. 
        #we expect good alignments since we align a sequence with its uniprot reference
        #but setting the base alignment score low can only be beneficial.
        best_score = -999999
        best_match = []
        #if the uniprot query obtained a propper sequence it will align each variant within the group with 
        #this reference sequence 
        if refDict[geneName] != '':
            for var in varDict:
                header = 'Alignment of: ' + var + ' with reference sequence for ' + groupID + '\n'
                try:
                    alignment = pairwise2.align.globalms(refDict[geneName], varDict[var], 1, -4, -0.5, -0.1, one_alignment_only=True)[0]
                    out.writelines([header, format_alignment(*alignment), '\n\n'])
                    score = alignment[2]
                    normscore = round((score / len(alignment[0])) * 1000) # max alignment score = 1000 since (len(alignment)*1 / len(alignment)) = 1000 (1 for all perfect matches)
                    scoreList.append(normscore )
                    group_alignments.append([score, geneName, var, alignment])
                # if an out of memory exception is thrown because for example the sequences are too long, then
                # it will handle the alignment this way and return an NA for score and give for alignment 
                # the variant sequence with itself which is a one to one match and will thus return the original mutation position.
                # the variant for which the exception ocurred is registered
                except MemoryError:
                    scoreList.append('NA')
                    alignment = [varDict[var],varDict[var]]
                    group_alignments.append(['NA', geneName, var, alignment])
                    alignErrors.append(var)
                    
        #if the query to uniprot returned no match the reference sequence will be an empty string '' and you enter the else case:
        #alignment scores will be NA's since there's nothing to align to.
        else:
            for var in varDict:
                try:
                    alignment = pairwise2.align.globalms(varDict[var], varDict[var], 1, -4, -0.5, -0.1, one_alignment_only=True)[0] # gap penalty -4 gap extention -0.5 match +5
                    out.writelines([header, format_alignment(*alignment), '\n\n'])
                    scoreList.append('NA')
                    group_alignments.append(['NA', geneName, var, alignment])
                #if an out of memory exception is thrown because for example the sequences are too long, then
                #it will handle the alignment this way and return an NA for score 
                except MemoryError:
                    scoreList.append('NA')
                    alignment = [varDict[var], varDict[var]]
                    group_alignments.append(['NA', geneName, var, alignment])
        
        # determining which alignment was the best and should be representative for the group
        for a in group_alignments:
            if a[0] == 'NA':
              pass  # ignore NA's
            elif a[0] > best_score:
                best_score = a[0]
                best_match = [a[2], a[3]]
        # if at this point the best_match list of a group is still empty, this means all scores in the group were NA.
        # This can happen when the uniprot query failed to return anything leaving the reference sequence '' 
        # or when a memory exception ocurred during the alignment. 
        # In this case we fill it with the variant name and alignment of the last variant of the group.
        if best_match == []:
            best_match = [a[2], a[3]]
        
        group_bestmatch.append(best_match)    
        printProgressBar(i + 1, len(groupVarDict), prefix='', suffix='Complete'+str(normscore), length=50)
        i += 1
    out.close
    return [group_bestmatch, scoreList, alignErrors]


def getAbsMutpos(group_bestmatch, groups, refDict):
    mutposList = []
    for bestmatch in group_bestmatch:
        var = bestmatch[0]
        alignment = bestmatch[1]
        # find mutation position of the variant
        for group in groups:
            for row in group:
                if str(var) == str(row[1]['variant']):
                    var_mutpos = int(row[1]['position_mutation'])
                    break
        align_ref = str(alignment[0])
        align_var = str(alignment[1])
        i = 0
        mutpos = var_mutpos
        while i < mutpos:
            if align_var[i] == "-":
                mutpos += 1
                i += 1
            else:
                i += 1
        new_mutpos = mutpos

        gapcount = 0
        for i in range(
                new_mutpos - 1):  # -1 because we go up to but not including the mutation position in the variant sequence.
            if align_ref[i] == '-':
                gapcount += 1
        if ('-' in align_ref):
            gapcount += 1  # increase gap count by 1 to adjust for the gap introduced by the mutation.
        # There are still exceptional cases in the final report where the variant IS the uniprot reference which means there are no gaps at all        
        abs_mutpos = new_mutpos - gapcount

        for group in groups:
            for row in group:
                if str(var) == str(row[1]['variant']):
                    for i in range(len(group)):
                        mutposList.append(abs_mutpos)
                    break
    return mutposList


# printing a progress bar
def printProgressBar(iteration, total, prefix='', suffix='', decimals=1, length=100, fill='█', printEnd="\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end=printEnd)
    # Print New Line on Complete
    if iteration == total:
        print('\n')
    return


# the main loop
reports, headfolder = find_files()

for x in reports:
    print('\nStarting '+x)
    current_report = reports[x]
    groups = group_duplicates(current_report, x)
    groupVarDict = find_mutseq(groups, x)
    
    #querying uniprot
    refseqs = find_refseq(groups)
    queryErrors = refseqs[2]
    if queryErrors == []: 
        print('completed without errors')
    else: 
        for error in queryErrors: print(error)
    refDict = refseqs[0]
    current_report['uniprot_id'] = refseqs[1]  # add uniprot_id column
    
    #aligning
    alignments = align_groups(refDict, groupVarDict, x)
    alignErrors = alignments[2]
    if alignErrors == []: 
        print('completed without errors')
    else: 
        for var in alignErrors:
            print('failed alignment for ' + var + ': MemoryException')
    group_bestmatch = alignments[0]
    scoreList = alignments[1]
    current_report['alignment_score'] = scoreList
    
    mutposList = getAbsMutpos(group_bestmatch, groups, refDict)
    current_report['abs_mutpos'] = mutposList

    current_report.to_csv(headfolder + '/' + x + '/' + 'finalreport2.tab', sep="\t", quoting=csv.QUOTE_NONE,
                          header=True, index=False)

print('Done')