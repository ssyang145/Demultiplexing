f=open('16S_library_barcodes_list.txt')
fw=open('unmapped_stats_with_16S_library_barcodes_list.txt','w')
count_table={}
total=0
barcode_list=[]
for line in f:
	barcode_list.append(line[:-1])

print len(barcode_list), "total num of 16S library barcodes"
fw.write(str(len(barcode_list)) + " total num of 16S library barcodes\n")
f.close()

from Bio.Seq import Seq

f1=open('unmapped_stats.txt')

i=0
j=0
f1.readline()
f1.readline()
f1.readline()
for total, line in enumerate(f1):
	total+=1
	barcode=line[:-1].split('\t')[0]
	reads=line[:-1].split('\t')[1]

	if barcode in barcode_list:
		print barcode, "barcode found", reads
		fw.write(barcode+ '\tbarcode found\t'+reads+'\n')
		i+=1

	seq=Seq(barcode)
	rev=seq.reverse_complement()
	if rev in barcode_list:
		print barcode, "reverse complementary barcode found", reads
		fw.write(barcode+ '\treverse complementary barcode found\t'+reads+'\n')
		j+=1

print total, "total num of unmapped barcodes"
print i, "total num of barcodes found"
print j, "total num of reverse complementary barcode found"

fw.write(str(total) + " total num of unmapped barcodes\n")
fw.write(str(i) + " total num of barcodes found\n")
fw.write(str(j) + " total num of reverse complementary barcode found\n")

f1.close()
fw.close()