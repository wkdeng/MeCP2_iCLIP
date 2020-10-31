import sys
fin=sys.argv[1]
fo=sys.argv[2]

barcode_dict = {'CGTGAT': 'MC',
                'ACATCG':'M1',
                'GCCTAA':'M2',
                'TGGTCA': 'M3',
                'CACTGT': 'H1',
                'ATTGGC': 'H2',
                'GATCTG': 'H3'}
fo_writers={
    'MC':open(fo+'/MC.fq','w'),
    'M1': open(fo+'/M1.fq', 'w'),
    'M2': open(fo+'/M2.fq', 'w'),
    'M3': open(fo+'/M3.fq', 'w'),
    'H1': open(fo+'/H1.fq', 'w'),
    'H2': open(fo+'/H2.fq', 'w'),
    'H3': open(fo+'/H3.fq', 'w'),
    'Undetermined': open(fo+'/Undetermined.fq', 'w')
}

last_read_content=[]
last_read_group=''
line_count=0
read_count=0
for line in open(fin):
    # if read_count>=10000:
    #     break
    if line_count==3 and len(last_read_content) >0:
        last_read_content.append(line[15:])
        fo_writers[last_read_group].write('\n'.join(last_read_content))
        line_count=0
        last_read_group=''
        read_count+=1
        last_read_content=[]
    else:
        line_count+=1
        if line_count==2:
            barcode=line[5:11]
            umi=line[:15]
            last_read_content[0] = '@%s:%s'%(last_read_content[0][1:].split(' ')[0],umi)
            if barcode in barcode_dict:
                last_read_group=barcode_dict[barcode]
            else:
                last_read_group='Undetermined'
            line=line[15:]
        last_read_content.append(line.strip())


for key in fo_writers:
    fo_writers[key].close()
