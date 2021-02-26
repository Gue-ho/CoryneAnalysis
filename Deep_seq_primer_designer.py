import os

ref_seq = open('../Ref/Ref.fasta').readlines()[1].strip()

pri_d = {}
res_d = {}

def rc(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

fw = open('input.txt','w')
fw.write('/mnt/d/Project/SKKU_WHM/WHS/Primer_design/Ref/\n'+'N'*20+'\n')

offinder_d = {}

out_p = []
for i in open('outlist.txt').readlines():
	out_p.append(int(i.strip().split('\t')[1]))

for i in open('list3.txt').readlines():
		
	isp = i.split('\t')
	refnt = isp[2][0]
	p = int(isp[1])
	if p in out_p:
		print(isp)
		print(p)
		print('!!!')
		input()

	if ref_seq[p-1] != refnt:
		print(i)
	
	res_d[isp[1]] = isp
	pri_d[isp[1]] = [[],[]]

	seq = ref_seq[p - 1 - 130: p + 130]

	for i in range(50):
	
		s = seq[i:i+20]
		if 8 < s.count('G') + s.count('C') < 12:
			pri_d[isp[1]][0].append(s)
			if s in offinder_d.keys():
				pass
			else:
				offinder_d[s] = [0,0,0]
				fw.write(s+' 2\n')

		s = rc(seq[len(seq)-i-20: len(seq)])
		if 8 < s.count('G') + s.count('C') < 12:
			pri_d[isp[1]][1].append(s)
			if s in offinder_d.keys():
				pass
			else:
				offinder_d[s] = [0,0,0]
				fw.write(s+' 2\n')

print(len(res_d))

fw.close()

os.system('cas-offinder input.txt G output.txt')

for i in open('output.txt').readlines():
	isp = i.split('\t')
	offinder_d[isp[0]][int(isp[5])] += 1


fw = open('result3.txt','w')

for p, isp in res_d.items():
	
	for_pri = ''
	for pri in pri_d[isp[1]][0]:
		if offinder_d[pri] == [1,0,0]:
			for_pri = pri
			break
	
	back_pri = ''
	for pri in pri_d[isp[1]][1]:
		if offinder_d[pri] == [1,0,0]:
			back_pri = pri
			break
	
	if for_pri == '':
		for pri in pri_d[isp[1]][0]:
			if offinder_d[pri][:2] == [1,0]:
				for_pri = pri
				break
	
	if back_pri == '':
		for pri in pri_d[isp[1]][1]:
			if offinder_d[pri][:2] == [1,0]:
				back_pri = pri
				break
		
	seq = pri_seq[p-1-130: p+130]
	
	isp += [for_pri, back_pri, seq[seq.find(for_pri):seq.rfind(back_pri)+20], str(len(seq[seq.find(for_pri):seq.rfind(back_pri)+20]))]
	fw.write('\t'.join(isp)+'\n')

fw.close()


