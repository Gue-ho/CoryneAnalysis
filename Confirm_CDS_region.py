l = 1000

ID_db = {}
ID_tocg = {}

for i in open('ID_list.txt').readlines():
	isp = i.strip().split('\t')
	if isp[1] != '-':
		ID_db[isp[0]] = isp[1]
		ID_tocg[isp[1].upper()] = isp[0]
	elif isp[2] != '-':
		ID_db[isp[0]] = isp[2]
	else:
		ID_db[isp[0]] = isp[0]
	if isp[2] != '-':
		ID_tocg[isp[2].upper()] = isp[0]
	ID_tocg[isp[0].upper()] = isp[0]

TLS = []
TSS = {}

fl = open('NC_003450_TLS.txt').readlines()
fs = open('TLS_TSS_combine.txt').readlines()
ref = open('../../database/NC_003450_sequence.fasta').readlines()[1].strip()

amino_sym = {"TGG": "W", "GGG": "G", "GTT": "V", "CGG": "R", "ACT": "T", "TGA": "X", "CTA": "L", "TCC": "S","GAA": "E", "CCA": "P", "GAC": "D", "ACC": "T", "TTT": "F", "CTC": "L", "GCT": "A", "CCC": "P", "TCG": "S", "CAT": "H", "GTC": "V", "CGA": "R", "CAG": "Q", "ATA": "I", "AAG": "K", "CCG": "P", "GGA": "G", "AGC": "S", "TAT": "Y", "CTG": "L", "ACG": "T", "GAG": "E", "GCT": "A", "TGC": "C", "TGT": "C", "AGG": "R", "ATG": "M", "TTA": "L", "GCA": "A", "AAT": "N", "GTA": "V", "GGT": "G", "AGA": "R", "CGC": "R", "ATC": "I", "TAC": "Y", "TTG": "L", "ACA": "T", "GCG": "A", "CTT": "L", "ATT": "I", "CGT": "R", "CAC": "H", "TCA": "S", "CCT": "P", "TAA": "X", "GAT": "D", "GTG": "V", "AAA": "K", "AAC": "N", "GGC": "G", "TTC": "F", "CAA": "Q", "AGT": "S", "TAG": "X", "TCT": "S", "GCC": "A"}

def rc(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

prop = []
rprop = []
for i in range(len(ref)):
	prop.append(0)
	rprop.append(0)

remove_gene = []

for i in  open('operon.txt').readlines():
	isp = i.split()
	ope = isp[0].split('-')
	if len(ope) < 2:
		continue
	if isp[2] == '+':
		remove_gene += ope[1:]
	elif isp[2] == '-':
		remove_gene += ope[:-1]

for i in range(len(fl)):
	tls = fl[i].strip().split('\t')
	tss = fs[i].strip().split('\t')
	if tls[2] == '+':
		st = int(tls[1])
		ed = int(tls[3])
		for x in range(st - l, st):
			prop[x-1] += 1
		if tls[0].upper() not in ID_tocg.keys() or ID_tocg[tls[0].upper()] not in remove_gene:
			for x in range(st - l, st):
				rprop[x-1] += 1
	elif tls[2] == '-':
		ed = int(tls[1])
		st = int(tls[3])
		for x in range(ed, ed + l):
			try:
				prop[x-1] += 1
			except:
				pass
		if tls[0].upper() not in ID_tocg.keys() or ID_tocg[tls[0].upper()] not in remove_gene:
			for x in range(ed, ed + l):
				try:
					rprop[x-1] += 1
				except:
					pass
	
	TLS.append([st, ed, tls[0], tls[2]])
	sts = int(tss[1])

print('promoter_rate : {0}'.format(round((len(ref) - prop.count(0))*100/len(ref),2)))

print('operon_rate : {0}'.format(round((len(ref) - rprop.count(0))*100/len(ref),2)))

for i in open('coding_list.txt').readlines():
	isp = i.strip().split()
	if isp[0] in ID_db.keys():
		iid = ID_db[isp[0]]
	else:
		iid = isp[0]
	if len(isp) < 5: continue
	if isp[1] == '+':
		st = int(isp[2].replace(',',''))
		ed = int(isp[3].replace(',',''))
	elif isp[1] == '-':
		st = int(isp[3].replace(',',''))
		ed = int(isp[2].replace(',',''))
	if iid in TSS.keys():
		if int(isp[5]) > TSS[iid][3]:
			TSS[iid] = [isp[0], st, ed, int(isp[5])]
	else:
		TSS[iid] = [isp[0], st, ed, int(isp[5])]

fw=open('result_20200525.txt','w')
for i in open('pos.txt').readlines():
	isp = i.strip().split('\t')
	ip = int(isp[0])
	ref_nt = isp[1][0]
	mut_nt = isp[1][3]
	gene_name = ''
	for x, y in TSS.items():
		if y[1] <= ip <= y[2]:
			gene_name += 'UTR_' + x+';'
	gene_name += '\t'
	pro_name = ''
	ope_name = ''
	mut_type = ''
	aa_mut = ''
	pos_percent = ''
	for x in TLS:
		if (x[3] == '+' and x[0] - l <= ip <x[0]) or (x[3] == '-' and x[1] < ip <= x[1] + l):
			pro_name += 'Pro_' + x[2] + ';'
			if x[2].upper() not in ID_tocg.keys() or ID_tocg[x[2].upper()] not in remove_gene:
				ope_name += 'Ope_' + x[2] + ';'
		if x[0] <= ip <= x[1]:
			gene_name += 'CDS_' + x[2] + ';'
			gene_seq = ref[x[0] - 1: x[1]]
			gp = ip - x[0]
			if x[3] == '+':
				pos_percent = round((ip - x[0]) * 100 / (x[1] - x[0]), 2)
			else:
				pos_percent = round((x[1] - ip) * 100 / (x[1] - x[0]), 2)
			if x[3] == '-':
				gene_seq = rc(gene_seq)
				gp = x[1] - ip
				ref_nt = rc(ref_nt)
				mut_nt = rc(mut_nt)
			mut_seq = gene_seq[:gp] + mut_nt + gene_seq[gp+1:]
			ap = (gp-1) // 3
			if gene_seq[gp] != ref_nt:
				print(gp)
				print(x)
				print(i)
				input()
			ref_aa = amino_sym[gene_seq[3*ap: 3*ap+3]]
			mut_aa = amino_sym[mut_seq[3*ap: 3*ap+3]]
			if ref_aa != 'X' and mut_aa == 'X':
				mut_type = 'Nonsense'
			elif ref_aa == mut_aa:
				mut_type = 'Silent'
			elif ref_aa != mut_aa:
				mut_type = 'Missense'
			aa_mut = ref_aa + 'to' + mut_aa
	fw.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n'.format(i.strip(), gene_name, mut_type, aa_mut, pos_percent, ref[ip-1 -5: ip-1 + 6], rc(ref[ip-1-5:ip-1+6]), pro_name, ope_name))
fw.close()
