

def a(fn):
	
	all_cnt = 3309401
	m_cnt = 0

	for i in open(fn).readlines():
		isp = i.strip().split('\t')
		if len(isp)< 2:
			continue
		if isp[2] != '0':
			m_cnt += 1
	
	print(all_cnt)
	print(m_cnt)
	print(round(m_cnt*100/all_cnt, 2))

for i in ['WT1','WT2']:
	a('{0}.depth'.format(i))

