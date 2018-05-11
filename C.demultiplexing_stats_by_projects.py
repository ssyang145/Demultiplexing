import re

f1=open('demultiplexing_stats.txt') # demultiplexing stats file
f1.readline()
f1.readline()
f1.readline()
header=f1.readline()

i=0
data_table={}
project_list=[]
project=False
data_table[project]={}
for line in f1:
	i+=1
	project=line[:-1].split("\t")[0].replace(" ","")
	sample=line[:-1].split("\t")[1].replace(" ","")
	#print sample
	reads=line[:-1].split("\t")[3].replace(" ","")
	if project not in project_list:
		data_table[project]={}
		data_table[project][sample]=[reads]
		project_list.append(project)
	else:
		data_table[project][sample]=[reads]
				

#print data_table

project_list=set(project_list)
print project_list

for project in project_list:
	fw=open('demultiplexing_stats_'+project+'.txt','w') # write demultiplexing stats by project
	fw.write('Project\tSample\tReads\n')
	for sample in data_table[project]:
		fw.write(project+'\t'+sample+'\t'+'\t'.join(data_table[project][sample])+'\n')
	#data_table[sampleID]=[ilabID,sampleID,barcode]
	fw.close()

#print len(data_table)
print i	
f1.close()




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
