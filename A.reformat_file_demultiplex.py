import re
data_table = {}
run_num='001' # change everytime

f1=open('barcodes_samples_mb16s_'+run_num+'.txt') # read original barcodes file
fw=open('reformat_barcodes_samples_mb16s_'+run_num+'.txt','w') # write reformatted file
header=f1.readline()
fw.write('Run Name\tmb16s_'+run_num+'\n')
fw.write('Data Folder\t.\n')
fw.write('R1 File Name\tUndetermined_S0_L001_R1_001.fastq.gz\n')
fw.write('R2 File Name\tUndetermined_S0_L001_R2_001.fastq.gz\n')
fw.write('Index File Name\tUndetermined_S0_L001_I1_001.fastq.gz\n')
fw.write('\n')
fw.write('Project\tSampleID\tBarcode\tLinker/Primer\tDescription\n')

i=0
for line in f1:
	#print [line]
	barcode=line[:-1].split("\t")[0].replace(" ","")
	sampleID=line[:-1].split("\t")[1].replace(" ","")
	sampleID=re.sub('[^a-zA-Z0-9_-]','',sampleID)
	i+=1
	sampleID_index=sampleID+'_'+str(i)
	#print sampleID_index
	ilabID=line[:-1].split("\t")[2].replace(" ","")
	ilabID=re.sub('[^a-zA-Z0-9_-]','',ilabID)
	#print ilabID
	fw.write(ilabID+'\t'+sampleID_index+'\t'+barcode+'\t\t'+ilabID+'\n')
	#data_table[sampleID]=[ilabID,sampleID,barcode]

#print len(data_table)
print i	
f1.close()
fw.close()
# read data file and add information


'''# with \n

f = open(file )
fw=open('des_barcodes_samples_mb16s_'+file,'w')
print file
	#print [f.read()[:500]]
	#a

header1 = f.readline()
fw.write('ID\t'+header1[:-1]+'\tgene\tlocus\n')
count=0
n=0
for line in f:
	fw.write(line[:-1])
	id = line[:-1].split('\t')[0]
	if id in data_table:
		count +=1
		fw.write('\t' + '\t'.join(data_table[id]))
		name=data_table[id][0]	

	fw.write('\n')
print count
print n

fw.close()
f.close()'''
