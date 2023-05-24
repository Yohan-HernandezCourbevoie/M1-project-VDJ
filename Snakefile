import os
import glob
import json
import csv


configfile: "config.json"

DATA_DIR=config["data_dir"]
ACCESSION_FILES=config['sra_accession']
FLASH_DIR=config["flash_dir"]
VIDJIL_DIR=config["vidjil_dir"]
CSV_DIR=config["csv_dir"]
FASTA_DIR=config["fasta_dir"]
CHAIN_ALIGNMENT_DIR=config["chain_alignment_dir"]
STAT_DIR=config["stat_dir"]
RESULTS_DIR=config["results_dir"]
ACCESSION_Project=config["accession_project"]

SAMPLES=[accession.rstrip() for accession in open(ACCESSION_FILES).readlines()]
CHAINS = ['TRA','TRB','TRG','TRD','IGH','IGK','IGL']


rule all:
    input:
        results = expand(RESULTS_DIR + "/results_{chain}.txt", chain=CHAINS),
        json= expand(RESULTS_DIR + "/results_{chain}.json", chain=CHAINS),
        stat = STAT_DIR + "/stat.txt"

def sra_url(srr):
    prefix=srr[:6]
    url=''
    n = len(srr)
    if n == 9:
        url = srr+"/"
    elif n > 9:
        end=srr[9:].zfill(3)
        url=end+"/"+srr+"/"
    return "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"+prefix+"/"+url

def vidjil2csv(filepath, accession, csv_name, file_project):
    projects = {}
    with open(file_project, "r") as file_in:
        for line in file_in:
            link = line.split()
            projects[link[0]] = link[1]

    jsonstr = ""
    with open(filepath, "r") as fileIn:
        for line in fileIn.readlines():
            jsonstr += line.strip()
    pythonstr = json.loads(jsonstr)

    with open(csv_name, 'w', newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
        if pythonstr['clones'] is not None:
            for clone in pythonstr['clones']:
                if accession in projects:
                    spamwriter.writerow([accession, projects[accession], clone['germline'], clone['id'],
                                         clone['reads'][0], clone['sequence']])
                else:
                    spamwriter.writerow([accession, 'not in list', clone['germline'], clone['id'], clone['reads'][0],
                                         clone['sequence']])
        else:
            pass

def recup_chain(chain, file_in):
    consensus_list = []
    prev_window = None
    prev_project = []
    temp_consensus = []
    accessions_list = []
    window_list = []
    temp_window = []
    temp_accessions = []
    with open(file_in) as fileIn:
        for line in fileIn:
            information = line.split(',')
            if prev_window is None:
                if information[2] == chain:
                    prev_window = information[3]
                    temp_consensus.append(information[5])
                    prev_project.append(information[1])
                    temp_accessions.append(information[0])
                    temp_window.append(information[3])
            else:
                if prev_window == information[3]:
                    if information[2] == chain:
                        if information[1] not in prev_project:
                            temp_consensus.append(information[5])
                            temp_accessions.append(information[0])
                            temp_window.append(information[3])
                            prev_project.append(information[1])
                else: # prev_window != fenetre
                    if len(temp_consensus) > 1:
                        consensus_list += temp_consensus
                        accessions_list += temp_accessions
                        window_list += temp_window
                    prev_window = information[3]
                    temp_consensus = [information[5]]
                    prev_project = [information[1]]
                    temp_accessions = [information[0]]
                    temp_window= [information[3]]
    if len(temp_consensus) > 1:
        consensus_list += temp_consensus
        accessions_list += temp_accessions
        window_list += temp_window
    return consensus_list, accessions_list, window_list

def consensus_to_fasta(consensus_list, accessions_list, window_list, fasta_out, id_window_out):
    dico = {}
    with open(fasta_out,'w') as out:
        for i in range(len(consensus_list)):
            out.write(f'> {accessions_list[i]}_{i}\n')
            out.write(consensus_list[i])
    for i in range(len(consensus_list)):
        if window_list[i] not in dico:
            dico[window_list[i]] = [i]
        else:
            dico[window_list[i]].append(i)
    with open(id_window_out, 'w') as out:
        for key in dico.keys():
            str_out = f''
            for i in dico[key]:
                str_out += f'{i}   '
            out.write(f'{key}   {str_out}\n')

def all_fasta(file_in):
    for type in CHAINS:
        list_consens, list_access , list_window= recup_chain(type, file_in)
        consensus_to_fasta(list_consens, list_access, list_window , f'{FASTA_DIR}/{type}.fasta', f'{FASTA_DIR}/{type}.txt')

def align2result(filein):
    jsonstr = ""
    with open(filein, "r") as file_in:
        for line in file_in:
            jsonstr += line.strip()
        pythonstr = json.loads(jsonstr)
    V, D, J = {}, {}, {}
    VDJ = {}
    vdj_indel = {}
    for i in range(len(pythonstr['clones'])):
        str_recombi = ''
        recombi_indel = ''
        if pythonstr['clones'][i]['germline'] != "not analyzed":
            recup_indel = pythonstr['clones'][i]['name'].split(' ')
            for part in recup_indel:
                recombi_indel += part.split('*')[0] + ' '
            if recombi_indel[:-1] not in vdj_indel:
                vdj_indel[recombi_indel[:-1]]={'view':1, 'window_index':[int(pythonstr['clones'][i]['id'])-1]}
            else:
                vdj_indel[recombi_indel[:-1]]['view'] += 1
                vdj_indel[recombi_indel[:-1]]['window_index'].append(int(pythonstr['clones'][i]['id'])-1)
            for j in range(3, 6):
                if str(j) in pythonstr['clones'][i]['seg']:
                    gene = pythonstr['clones'][i]['seg'][str(j)]['name'].split('*')[0]
                    if j == 3:
                        if gene in J:
                            J[gene] += 1
                        else:
                            J[gene] = 1
                    if j == 4:
                        if gene in D:
                            D[gene] += 1
                        else:
                            D[gene] = 1
                    if j == 5:
                        if gene in V:
                            V[gene] += 1
                        else:
                            V[gene] = 1
    sort_v = dict(sorted(V.items(), key=lambda x: x[1], reverse=True))
    sort_d = dict(sorted(D.items(), key=lambda x: x[1], reverse=True))
    sort_j = dict(sorted(J.items(), key=lambda x: x[1], reverse=True))
    sort_vdj_indel = dict(sorted(vdj_indel.items(),key=lambda x: x[1]['view'],reverse=True))
    return sort_v, sort_d, sort_j, sort_vdj_indel

def result2json(fileout, v, d, j, vdj_indel, id2window):
    id_window = {}
    with open(id2window, 'r') as fileIn:
        for line in fileIn:
            info = line.split()
            temp_window = info[0]
            del info[0]
            id_window[temp_window] = info
            id_window[temp_window]= list(map(int, id_window[temp_window]))
    prepare_json = {'recombination':[], 'V':[], 'D':[], 'J': []}
    for key in vdj_indel.keys():
        window = []
        for i in vdj_indel[key]['window_index']:
            for index in id_window.items():
                if i in index[1] and index[0] not in window :
                    window.append(index[0])
        prepare_json['recombination'].append({'name': key, 'view':vdj_indel[key]['view'], 'window' : window})
    for key_v in v.keys():
        prepare_json['V'].append({key_v: v[key_v]})
    for key_d in d.keys():
        prepare_json['D'].append({key_d: d[key_d]})
    for key_j in j.keys():
        prepare_json['J'].append({key_j: j[key_j]})
    with open(fileout, 'w') as fileOut:
        fileOut.write(json.dumps(prepare_json, indent=3))

def write_result(fileout, v, d, j, vdj_indel):
    with open(fileout, 'w') as file_out:
        file_out.write(f'Gene V ({len(v)}) :\n')
        line = ''
        for key in v.keys():
            line += f"'{key}' : {v[key]}\n"
        file_out.write(f'{line}\nGene D ({len(d)}) :\n')
        line = ''
        for key in d.keys():
            line += f"'{key}' : {d[key]}\n"
        file_out.write(f'{line}\nGene J ({len(j)}) :\n')
        line = ''
        for key in j.keys():
            line += f"'{key}' : {j[key]}\n"
        file_out.write(f'{line}\nV(D)J recombinations with indel ({len(vdj_indel)}) :\n')
        line = ''
        for key in vdj_indel.keys():
            line += f"'{key}' : {vdj_indel[key]['view']}\n"
        file_out.write(line)


rule get_data:
    resources:
        tmpdir = "/data/tmp"
    output:
        sample1=DATA_DIR+"/{accession}_1.fastq.gz", sample2=DATA_DIR+"/{accession}_2.fastq.gz"
    run:
        shell("wget -P {DATA_DIR} -N "+sra_url(wildcards.accession)+"/{wildcards.accession}_{{1,2}}.fastq.gz")

rule run_flash:
    resources:
        tmpdir = "/data/tmp"
    input:
        sample1=DATA_DIR+"/{accession}_1.fastq.gz", sample2=DATA_DIR+"/{accession}_2.fastq.gz"
    output:
        temp(FLASH_DIR+"/{accession}.notCombined_1.fastq.gz"), temp(FLASH_DIR+"/{accession}.hist"), temp(FLASH_DIR+"/{accession}.histogram"),
        extended=FLASH_DIR+"/{accession}.extendedFrags.fastq.gz", notcombined=FLASH_DIR+"/{accession}.notCombined_2.fastq.gz"
    shell:
        '''cd ~/FLASH-1.2.11; ./flash -z -d {FLASH_DIR} {DATA_DIR}/{wildcards.accession}_1.fastq.gz {DATA_DIR}/{wildcards.accession}_2.fastq.gz -o {wildcards.accession}; cd ~'''

rule concat_flash:
    resources:
        tmpdir="/data/tmp"
    input:
        extended=FLASH_DIR+"/{accession}.extendedFrags.fastq.gz", notcombined=FLASH_DIR+"/{accession}.notCombined_2.fastq.gz"
    output:
        regroup=FLASH_DIR+"/{accession}_regroup.fastq.gz"
    shell:
        '''zcat {input.extended} {input.notcombined} | gzip -c -> {FLASH_DIR}/{wildcards.accession}_regroup.fastq.gz; rm -f {input.extended} {input.notcombined}'''
        
rule run_vidjil:
    resources:
        tmpdir = "/data/tmp"
    input:
        regroup=FLASH_DIR+"/{accession}_regroup.fastq.gz"
    output:
        vidjil_regroup = VIDJIL_DIR+"/{accession}.regroup.vidjil"
    shell:
        '''cd ~/vidjil-algo-2021.04; ./vidjil-algo -c clones -g germline/homo-sapiens.g:TRA,TRB,TRG,TRD,IGH,IGK,IGL -o {VIDJIL_DIR} -b {wildcards.accession}.regroup {input.regroup} -y all; cd ~'''

rule vidjil2csv:
    resources:
        tmpdir = "/data/tmp"
    input:
        vidjil_regroup = VIDJIL_DIR+"/{accession}.regroup.vidjil"
    output:
        csv_regroup = CSV_DIR + "/{accession}.regroup.csv"
    run:
        vidjil2csv(f"{input.vidjil_regroup}", f"{wildcards.accession}", f"{CSV_DIR}/{wildcards.accession}.regroup.csv", f'{ACCESSION_Project}')

rule concat_csv:
    resources:
        tmpdir="/data/tmp"
    input:
        csv_regroup = expand(CSV_DIR + "/{sample}.regroup.csv", sample=SAMPLES)
    output:
        data = CSV_DIR + "/data.csv"
    shell:
        '''cat {input.csv_regroup} > {CSV_DIR}/data.csv'''

rule sort_data_windows:
    resources:
        tmpdir="/data/tmp"
    input:
        data = CSV_DIR + "/data.csv"
    output:
        data_sort = CSV_DIR + "/data_sort.csv"
    shell:
        '''sort -t, -k4,4 {input.data} > {output.data_sort} ; rm {input.data}'''

rule statistic:
    resources:
        tmpdir="/data/tmp"
    input:
        data = CSV_DIR + "/data_sort.csv"
    output:
        stat =  STAT_DIR + "/stat.txt"
    shell:
        '''python3 stats.py -f {input.data} -o {output.stat}'''

rule create_multifasta:
    resources:
        tmpdir="/data/tmp"
    input:
        data = CSV_DIR + "/data_sort.csv"
    output:
        fasta_file = FASTA_DIR + "/{chain}.fasta",
        txt_file = FASTA_DIR + "/{chain}.txt"
    run:
        all_fasta(f"{input.data}")

rule chain_alignment:
    resources:
        tmpdir="/data/tmp"
    input:
        data = FASTA_DIR + "/{chain}.fasta"
    output:
        chain_alignment = CHAIN_ALIGNMENT_DIR + "/{chain}.vidjil"
    shell:
        '''cd ~/vidjil-algo-2021.04; ./vidjil-algo -c designations -g germline/homo-sapiens.g:TRA,TRB,TRG,TRD,IGH,IGK,IGL -o {CHAIN_ALIGNMENT_DIR} {input.data}; cd ~'''

rule result_recombi:
    resources:
        tmpdir="/data/tmp"
    input:
        data = CHAIN_ALIGNMENT_DIR + "/{chain}.vidjil"
    output:
        results = RESULTS_DIR + "/results_{chain}.txt",
        json = RESULTS_DIR + "/results_{chain}.json"
    run:
        v,d,j,vdj_indel = align2result(f'{input.data}')
        write_result(f'{output.results}',v,d,j,vdj_indel)
        result2json(f'{output.json}',v,d,j,vdj_indel,f'{FASTA_DIR}/{wildcards.chain}.txt')