import os, operator, sys

def SamSub(filename):
	ref = open('../../../IGV/Ref.fa').readlines()[1].strip()
	#ref = open('../WT/mWT.fasta').readlines()[1].strip()
	with open(filename+'.sorted.sam') as f, open(filename+'_sub.txt','w') as fw:
		cnt = 0
		pos_d = {}
		m = 0
		for line in f:
			if line[0] == '@':
				continue
			a = line.split('\t')
			st = int(a[3])
			md = a[12]
			mut_seq = a[9]
			md = md[md.find('Z:')+2:].strip()
			pos_mut = 0
			pos_ref = st - 1
			p = ''
			del_v = False
			cigar = a[5]
			ins_p = 0
			cut_v = False
			ins_l = []
			
			for i in cigar:
				p += i
				try: 
					int(p)
				except:
					p = int(p[:-1])
					if i in ['S'] and cut_v == False:
						mut_seq = mut_seq[p:]
					cut_v = True
					if i == 'M':
						ins_p += p
					if i == 'I':
						for x in range(p):
							ins_l.append(ins_p+x+1)
						ins_p += p
					p = ''
			m +=1
			for i in md:
				p += i
				if del_v ==True and i not in 'ATGC':
					del_v = False
				try:
					int(p)
				except:
					p = p[:-1]
					if p == '':
						p = 0
					else:
						p = int(p)
					if i == '^':
						pos_mut += p
						pos_ref += p
						del_v = True
					if i in 'ATGC':
						if del_v == True:
							pos_ref += 1
							if ref[pos_ref-1] != i:
								print(md)
								print(pos_ref-1)
								print(ref[pos_ref-1])
								print(i)
								print('!!!! Del')
								input()
						else:
							
							pos_mut += p + 1
							pos_ref += p + 1
							cnt+=1
							ins_con = 0
							v = 0
							for x in range(len(ins_l)):
								if pos_mut +ins_con <= ins_l[x]:
									if pos_mut + ins_con == ins_l[x] or v == 1:
										if x !=len(ins_l)-1:
											if ins_l[x+1] - ins_l[x] != 1:
												print(pos_ref)
												if v== 1:
													ins_con += 1
												break
											else:
												v = 1
												ins_con += 1
									else:
										break
								else:
									ins_con += 1
							pos_mut += ins_con
							ins_l = ins_l[ins_con:]
							if ref[pos_ref-1] != i:
								print(md)
								print(pos_ref-1)
								print(ref[pos_ref-1])
								print(i)
								print('!!!! Sub')
								input()
							if i == mut_seq[pos_mut-1]:
								print(cigar[:20])
								print(ins_con)
								print(ins_l)
								print(pos_ref)
								print(ref[pos_ref-20:pos_ref])
								print(mut_seq[pos_mut-20:pos_mut])
								input()
							if pos_ref not in pos_d.keys():
								pos_d[pos_ref] = ['', i+'to'+mut_seq[pos_mut-1]]
							else:
								pos_d[pos_ref] = ['_dup', i+'to'+mut_seq[pos_mut-1]]
					p = ''
		fw.write('Total coutn : '+str(cnt)+'\n')
		for p, l in sorted(pos_d.items(), key=operator.itemgetter(0)):
			fw.write('{0}{1}\t{2}\n'.format(p,l[0], l[1]))

#list = [['Ldh', 'NC_003450_sequence'], ['BOL_1', 'NC_003450_sequence'], ['mWT', 'NC_003450_sequence']]
list = ['BOL1','BOL2','Idh1','Idh2','WT1','WT2']
#list = ['Ldh_mWT', 'BOL-1_mWT']
for n in list:
	SamSub(n)

"""
for l in list:
	os.system('bwa mem ../NC_003450/NC_003450_sequence.fasta {0}.fasta > {0}_{1}.bam'.format(l[0], l[1]))
	os.system('samtools sort {0}_{1}.bam > {0}_{1}_sorted.bam'.format(l[0], l[1]))
	os.system('samtools view -h -o {0}_{1}_sorted.sam {0}_{1}_sorted.bam'.format(l[0], l[1]))
	print(l)
	SamSub('{0}_{1}'.format(l[0], l[1])) 
"""
