#!/bin/bash

# RAxML bootstraps on Apodemia GBS SNP data
# Jeff Oliver
# jcoliver@email.arizona.edu
# 2016-09-30

# WARNING: Many hard-coded paths which likely do not exist
cd ~/Documents/Other\ Work/Apodemia
# 10 bootstraps
~/Executables/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -f a -p 201605121632 -x 201605151632 -s Data/Apodemia-2016-05-12.phyx -n Apo-RAxML-2016-05-12 -m GTRGAMMA -#10 -T 6

# 2 bootstraps, Dockweiler specimens as outgroup
~/Executables/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -f a -p 201605121632 -x 201605151632 -s Data/Apodemia-2016-05-12.phyx -n Apo-RAxML-2016-05-12 -m GTRGAMMA -#2 -T 6 -o Dockweiler_10501,Dockweiler_10502,Dockweiler_10503,Dockweiler_10504,Dockweiler_10505

# 10 bootstraps
cd ~/Documents/Other\ Work/Apodemia/Tree\ Inference/10\ BS
~/Executables/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -f a -p 201603021034 -x 201603021034 -s filtered_Apo.phyx -n filtered_Apo -m GTRGAMMA -#10 -T 6

# 100 bootstraps
cd ~/Documents/Other\ Work/Apodemia/Tree\ Inference/100\ BS
~/Executables/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -f a -p 201603021034 -x 201603021034 -s filtered_Apo.phyx -n filtered_Apo -m GTRGAMMA -#100 -T 6
