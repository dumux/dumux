rsync -zauP ~/Stupid lauser@login2.iws.uni-stuttgart.de:~/

ssh -p 1234 -l lauser localhost "FILTER='.*' /home/lauser/sb/run $@"
