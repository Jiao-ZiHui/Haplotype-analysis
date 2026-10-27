gpn = ["ruf", "niv"]
with open("new_sample_ID", 'r') as R:
    grp = {line.split('\t')[0]: line.split('\t')[1].rstrip('\n') for line in R}
with open('sample.vcf', 'r') as R, open('sample.fas', 'w+') as W:
    headline = ""
    tsam = []
    sam = []
    pdt = []
    for line in R:
        if line.startswith("#CHROM"):
            headline = line
        if line.startswith("#"):
            continue
        pdt.append(line.strip())
    pgp = [[] for _ in range(len(gpn))]
    #pgn用于接住样号名
    sam_pos = {}
    headline = headline.strip()
    tm = headline.split('\t')
    tsam = tm[9:]
    for i, sample in enumerate(tsam):
        sam_pos[sample] = i + 9
    # 获得样号的位置索引
    for tsm in tsam:
        if grp.get(tsm) in gpn:
            sam.append(tsm)
    # 如果样号属于目标群体，则加入sam列表
    for j in sam:
        for i, item in enumerate(gpn):
            if grp[j] == item:
                pgp[i].append(j)
    #对样号分组
    snp_pos = [] #snp位点
    hap_phs = {} #单倍型
    for sample in sam:
        hap_phs[sample + '.1'] = ""
        hap_phs[sample + '.2'] = ""
        #给每个样号命名两条链
    for i, var in enumerate(pdt):
        tm = var.split('\t')
        snp_pos.append(tm[1])
        snpref = tm[3]
        snpalt = tm[4].split(',')
        # 拆分0|1，判断单双位点
        if len(snpalt) == 1:
            for sample in sam:
                tn = tm[sam_pos[sample]].split('|')
                if tn[0] == '0' and tn[1] == '0':
                    hap_phs[sample + '.1'] += snpref
                    hap_phs[sample + '.2'] += snpref
                elif tn[0] == '0' and tn[1] == '1':
                    hap_phs[sample + '.1'] += snpref
                    hap_phs[sample + '.2'] += snpalt[0]
                elif tn[0] == '1' and tn[1] == '0':
                    hap_phs[sample + '.1'] += snpalt[0]
                    hap_phs[sample + '.2'] += snpref
                elif tn[0] == '1' and tn[1] == '1':
                    hap_phs[sample + '.1'] += snpalt[0]
                    hap_phs[sample + '.2'] += snpalt[0]
        else:
            for sample in sam:
                tn = tm[sam_pos[sample]].split('|')
                if tn[0] == '0':
                    hap_phs[sample + '.1'] += snpref
                elif tn[0] == '.':
                    hap_phs[sample + '.1'] += 'N'
                else:
                    hap_phs[sample + '.1'] += snpalt[int(tn[0])-1]
                if tn[1] == '0':
                    hap_phs[sample + '.2'] += snpref
                elif tn[1] == '.':
                    hap_phs[sample + '.2'] += 'N'
                else:
                    hap_phs[sample + '.2'] += snpalt[int(tn[1])-1]
    for key in hap_phs:
        W.write('>'+key+'\n'+hap_phs[key]+'\n')

