

# Installation instructions
# 1) Download hg38.fa file (938M)
# ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/

# 2) Install bedtools：
  # Way 1:
  # sudo apt-get update
  # apt install bedtools
  # Way 2:
  # apt-get https://github.com/arq5x/bedtools2/archive/bedtools-v2.30.0.tar.gz
  # tar xzvf bedtools-v2.30.0 
  # cd bedtools2
  # make
  # # sudo apt install build-essential # if missing gcc
  # # sudo apt-get install zlib1g-dev # if missing Zlib.h
  # cd bin/
  # export PATH=$PWD:$PATH

# 3) DNA sequences were obtained using bedtools
# sudo -s
bedtools getfasta -fi hg38.fa -bed allBED.bed > results.fa

  # 3) samtool tool is more convenient and the results are consistent with bedtools
  # samtools faidx hg38.fa chr4:122548477-122548979

# 4) Obtain TFs list from PROMO (a web tool)
# http://alggen.lsi.upc.es/cgi-bin/promo_v3/promo/promoinit.cgi?dirDB=TF_8.3

# 5) Determine the gene name
# https://www.genecards.org/
