for i in range(1,23):
    j = str(i)
    print "grep -P 'chr"+j+"\\t snp137.txt > chr"+j+"_snp137.txt"

print "grep -P 'chrX\\t snp137.txt > chrX_snp137.txt"
print "grep -P 'chrY\\t snp137.txt > chrY_snp137.txt"