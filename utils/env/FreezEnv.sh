

tstmp=$(date +%s)
dateRow=`date`

(pip3 freeze; echo ${dateRow} ) > ~/utils/env/pip_timestamp.${tstmp}.txt


Rscript ~/utils/env/rbackup.R
mv ~/utils/env/Renv.txt ~/utils/env/renv_timestamp.${tstmp}.txt
echo ${dateRow} >> ~/utils/env/renv_timestamp.${tstmp}.txt


