import os, sys, gzip

sample_file = sys.argv[1]

def rev_comp(seq): # reverse complement
    r_table = {'C':'G', 'G':'C', 'A':'T', 'T':'A', 'c':'g', 'g':'c', 'a':'t', 't':'a', 'Y':'R', 'R':'Y', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'D':'H', 'H':'D', 'V':'B', 'N':'N'}

    temp = ''
    for i in range(1, len(seq)+1):
        nt = seq[-i]
        if nt in r_table:
            temp += r_table[seq[-i]]
        else:
            temp += nt

    return temp

# read data info
f = open(sample_file)
text = f.read()
#print [text]

run_name = text.split('Run Name\t')[1].split()[0]
root_folder = text.split('Data Folder\t')[1].split()[0]
file_r1 = text.split('R1 File Name\t')[1].split()[0]
file_r2 = text.split('R2 File Name\t')[1].split()[0]
file_index = text.split('Index File Name\t')[1].split()[0]

if not os.path.exists(root_folder + '/' + run_name + '_demultiplexed'):
    os.mkdir(root_folder + '/' + run_name + '_demultiplexed')

# read barcode info
samples = {}
barcode_table = {}
des_table = {}
project_table = {}
for line in text.split('Linker/Primer\tDescription')[1].replace('\r', '\n').split('\n'):
    if line <> '':
        fields = line.split('\t')
        project, sample, barcode, primer, des = fields[:5]

        sample = sample.split()[0]
        barcode = rev_comp(barcode)
        
        samples[project + '\t' + sample] = barcode
        barcode_table[barcode] = project + '\t' + sample
        des_table[project + '\t' + sample] = des

        if project not in project_table:
            project_table[project] = {}

        project_table[project][sample] = True

samples = samples.keys()
samples.sort()

project_list = project_table.keys()
project_list.sort()

print len(project_list), '\tprojects'
print len(samples), '\ttotal samples'
print

# open files
sample_count_table = {}
R1_file_table = {}
R2_file_table = {}
for sample_text in samples:
    project, sample = sample_text.split('\t')
    #print project

    if not os.path.exists(root_folder + '/' + run_name + '_demultiplexed/' + project):
        os.mkdir(root_folder + '/' + run_name + '_demultiplexed/' + project)
        
    R1_file_table[sample_text] = gzip.open(root_folder + '/' + run_name + '_demultiplexed/' + project + '/' + sample + '_L001_R1_001.fastq.gz', 'wb')
    R2_file_table[sample_text] = gzip.open(root_folder + '/' + run_name + '_demultiplexed/' + project + '/' + sample + '_L001_R2_001.fastq.gz', 'wb')
    sample_count_table[sample_text] = 0
        
no_barcode_file_R1 = gzip.open(root_folder + '/' + run_name + '_demultiplexed/' + run_name + '_no_barcode_R1.fastq.gz', 'wb')
no_barcode_file_R2 = gzip.open(root_folder + '/' + run_name + '_demultiplexed/' + run_name + '_no_barcode_R2.fastq.gz', 'wb')
no_barcode_file_index = gzip.open(root_folder + '/' + run_name + '_demultiplexed/' + run_name + '_no_barcode_index.fastq.gz', 'wb')

# demultiplexing
if file_r1.endswith('.gz'):
    f1 = gzip.open(root_folder + '/' + file_r1, 'rb')
else:
    f1 = open(root_folder + '/' + file_r1)

if file_r2.endswith('.gz'):
    f2 = gzip.open(root_folder + '/' + file_r2, 'rb')
else:
    f2 = open(root_folder + '/' + file_r2)

if file_index.endswith('.gz'):
    f = gzip.open(root_folder + '/' + file_index, 'rb')
else:
    f = open(root_folder + '/' + file_index)


i = 1
text_r1 = ''
text_r2 = ''
text_index = ''
sample_found = False
unmapped_table = {}
total = 0
unmapped = 0
for line in f:
    text_r1 += f1.readline()
    text_r2 += f2.readline()
    text_index += line
 
    if i % 4 == 0:
        total += 1
        if sample_found:
            R1_file_table[sample_found].write(text_r1)
            R2_file_table[sample_found].write(text_r2)

        else:
            no_barcode_file_R1.write(text_r1)
            no_barcode_file_R2.write(text_r2)
            no_barcode_file_index.write(text_index)

            unmapped += 1

            if barcode not in unmapped_table:
                unmapped_table[barcode] = 1
            else:
                unmapped_table[barcode] += 1

        text_r1 = ''
        text_r2 = ''
        text_index = ''
        sample_found = False
        
    elif i % 4 == 2:
        barcode = line[:-1]
        
        if barcode in barcode_table:
            sample_found = barcode_table[barcode]
            sample_count_table[sample_found] += 1

    i += 1

    if i <> 4 and (i-4) % 2000000 == 0:
        print '%.1fM'%(total/1000000.0) + '\treads\t' '%.1f'%((total-unmapped)/float(total)*100) + '% mapped'

print '%.1fM'%(total/1000000.0) + '\treads\t' '%.1f'%((total-unmapped)/float(total)*100) + '% mapped'

# cloes files
for sample in samples:
    R1_file_table[sample].close()
    R2_file_table[sample].close()

no_barcode_file_R1.close()
no_barcode_file_R2.close()
no_barcode_file_index.close()

f1.close()
f2.close()
f.close()


# write demultiplexing stats
fw = open(root_folder + '/demultiplexing_stats.txt', 'w')
fw.write(str(total) + '\tTotal Reads\n')
fw.write(str(total - unmapped) + '\tDemultiplexed Reads\n\n')
fw.write('Project\tSample\tDescription\tReads\t% of Total\t% of Demultiplexed\n')
for sample in samples:
    fw.write(sample + '\t' + des_table[sample] + '\t' + str(sample_count_table[sample]) + '\t%.2f'%(sample_count_table[sample]/float(total)*100) + '%'
             + '\t%.2f'%(sample_count_table[sample]/float(total-unmapped)*100) + '%\n')
fw.close()

# write unmapped stats
fw = open(root_folder + '/unmapped_stats.txt', 'w')
fw.write(str(total) + '\tTotal Reads\n')
fw.write(str(unmapped) + '\tUnmapped Reads\n\n')

sort_table = {}
for barcode, count in unmapped_table.iteritems():
    if count not in sort_table:
        sort_table[count] = {}

    sort_table[count][barcode] = True

count_list = sort_table.keys()
count_list.sort()
count_list.reverse()

for count in count_list:
    if count >= 1000:
        for barcode in sort_table[count]:
            fw.write(barcode + '\t' + str(count) + '\n')

    else:
        break

fw.close()
