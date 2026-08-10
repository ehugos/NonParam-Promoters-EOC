#!/usr/bin/bash

# fullname="USER INPUT"
read -p "Enter path for downloading raw data: " fullname

wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133556/suppl/GSE133556_RAW.tar' -O fullname"/GSE133556/

#wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155760/suppl/GSE155760_RAW.tar' -O fullname"/GSE155760/

#wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185008/suppl/GSE185008_RAW.tar'
# wget 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226872/suppl/GSE226872_RAW.tar'