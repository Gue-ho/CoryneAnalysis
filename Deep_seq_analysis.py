import operator

sptoid = {}
stop = {}
fd = {}
order_l = []

for i in open('ana_list2.txt').readlines():
	isp = i.strip().split('\t')
	if i[0].find('Both') != -1:
		fd[isp[0]] = i.strip().split('\t') + ['','']
	else:
		fd[isp[0]] = i.strip().split('\t') + ['','','']
	sptoid[isp[5]+':'+isp[11]] = isp[0]
	order_l.append(isp[0])

	if isp[5] not in stop.keys():
		stop[isp[5]] = [isp[11]]
	else:
		stop[isp[5]].append(isp[11])

for i in ['23','84']:
	
	seq_count = {}
	print(str(i)+'...')
	with open('20200417/{0}.fastqjoin'.format(i)) as f:
		n = 0
		for line in f:
			n += 1
			if n%4 != 2:
				continue
			if line[:20] not in stop.keys():
				continue
			if line[:20] in seq_count.keys():
				seq_count[line[:20]].append(line)
			else:
				seq_count[line[:20]] = [line]
	for pri, seq_l in seq_count.items():
		for n in stop[pri]:
			wt_cnt = 0
			mut_cnt = 0
			all_cnt = 0
			sid = sptoid[pri+':'+n]
			m = fd[sid][2]
			ATGC = {'A':0, 'T':0, 'G':0, 'C':0}
			if len(seq_l) <500:
				continue
			for seq in seq_l:
				all_cnt += 1
				if len(seq) <= int(n):
					continue
				if seq[int(n)-1] == 'N':
					continue
				ATGC[seq[int(n)-1]] += 1
			nt = sorted(ATGC.items(), key=operator.itemgetter(1), reverse=True)[0][0]
			if nt == m[0]:
				nt = 'WT'
			elif nt == m[3]:
				nt = 'Mut'
			if sid.find('Both') != -1:
				if 13<=int(i)<=24:
					fd[sid][13] += str(nt)
				elif 37<=int(i)<=48:
					fd[sid][14] += str(nt)
				elif 85<=int(i)<=96:
					fd[sid][12] += str(nt)
			else:
				if nt == 'WT':
					fd[sid][12] += str(nt)
				elif nt == 'Mut':
					fd[sid][13] += str(nt)
fw=open('NGS_result2.txt','w')
for i in order_l:
	fw.write('\t'.join(fd[i])+'\n')
fw.close()
