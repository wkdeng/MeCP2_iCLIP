import numpy as np
barcodes={
'M_CTRL':'CGTGAT',
'M1':'ACATCG',
'M2':'GCCTAA',
'M3':'TGGTCA',
'H1':'CACTGT',
'H2':'ATTGGC',
'H3':'GATCTG'
}
barcode_position={}
for key in barcodes:
    barcode_position[key]=np.zeros(150)


undeterminded=0
line_count=0
total_read=0
for line in open('/home/dengw1/workspace/MeCP2_iCLIP/raw_fq/Undetermined_S0_L003_R1_001.fastq'):
    if line.startswith('@'):
        line_count=1
    else:
        line_count+=1
    if line_count==2:
        total_read+=1
        if not  total_read%10000000:
            print('Total reads: %s'%total_read)
        flag=False
        for key in barcodes:
            pos=line.find(barcodes[key])
            if pos>=0:
                barcode_position[key][pos]+=1
                flag=True
                break
        if not flag:
            undeterminded+=1


print(total_read)
print('Undetermined: %s' % undeterminded)
print(barcode_position)
