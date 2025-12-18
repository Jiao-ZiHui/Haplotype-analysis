#!/usr/bin/env python3
gpn = ["ruf", "niv"]
with open("sample_ID", 'r') as R:
    grp = {line.split('\t')[0]: line.split('\t')[1].rstrip('\n') for line in R}
with open('sample.hap.vcf', 'r') as R, open('sample.out', 'w') as W:
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
    sam_pos = {}
    headline = headline.strip()
    tm = headline.split('\t')
    tsam = tm[9:]
    for i, sample in enumerate(tsam):
        sam_pos[sample] = i + 9
    for tsm in tsam:
        if grp.get(tsm) in gpn:
            sam.append(tsm)
    for j in sam:
        for i, item in enumerate(gpn):
            if grp[j] == item:
                pgp[i].append(j)
    snp_pos = [] 
    hap_phs = {}
    for sample in sam:
        hap_phs[sample + '.1'] = ""
        hap_phs[sample + '.2'] = ""
    for i, var in enumerate(pdt):
        tm = var.split('\t')
        snp_pos.append(tm[1])
        snpref = tm[3]
        snpalt = tm[4].split(',')
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
                else:
                    hap_phs[sample + '.1'] += snpalt[int(tn[0])-1]
                if tn[1] == '0':
                    hap_phs[sample + '.2'] += snpref
                else:
                    hap_phs[sample + '.2'] += snpalt[int(tn[1])-1]
    hap_dat = [hap_phs[sample + '.1'] for sample in sam] + [hap_phs[sample + '.2'] for sample in sam]
    hap = [hap_dat[0]]
    for dat in hap_dat:
        if dat not in hap:
            hap.append(dat)



    hap_total = len(hap_dat)
    hap_num = {key: 0 for key in hap}
    for key in hap_dat:
        hap_num[key] += 1
    nhp = [item for item in hap if hap_num[item] / hap_total > 0.01]
    nhp_no = {item: f'H{i + 1}' for i, item in enumerate(nhp)}
    W.write("Variation of sites\n")
    for i, var in enumerate(pdt):
        W.write(f"{i + 1}\t{var}\n")
    W.write("\n")

    W.write(f"A total of {len(nhp)} gene haplotypes in rice\n")
    for i, item in enumerate(nhp):
        W.write(f"H{i + 1}\t{item}\n")
    W.write("\n")

    W.write(f"The distribution of gene haplotype in rice groups\n")
    W.write("Haplotype\tAll")
    for group in gpn:
        W.write(f"\t{group}")
    W.write("\n")
    for i, np in enumerate(nhp):
        W.write(f"H{i + 1}")
        htn = 0
        for sample in sam:
            if hap_phs[sample + '.1'] in nhp_no and hap_phs[sample + '.2'] in nhp_no:
                if hap_phs[sample + '.1'] == np:
                    htn += 1
                if hap_phs[sample + '.2'] == np:
                    htn += 1
        W.write(f"\t{htn}")
        for j, group in enumerate(gpn):
            htn = 0
            for sample in pgp[j]:
                if hap_phs[sample + '.1'] in nhp_no and hap_phs[sample + '.2'] in nhp_no:
                    if hap_phs[sample + '.1'] == np:
                        htn += 1
                    if hap_phs[sample + '.2'] == np:
                        htn += 1
            W.write(f"\t{htn}")
        W.write("\n")
    W.write("\n")

    W.write("The hyplotype of each sample in each group\n")
    for group in gpn:
        W.write(f"{group}\thaplotype\n")
        for sample in pgp[gpn.index(group)]:
            if hap_phs[sample + '.1'] in nhp_no and hap_phs[sample + '.2'] in nhp_no:
                W.write(f"{sample}_1\t{nhp_no[hap_phs[sample + '.1']]}\n")
                W.write(f"{sample}_2\t{nhp_no[hap_phs[sample + '.2']]}\n")