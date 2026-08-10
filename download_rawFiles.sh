#!/usr/bin/bash

# Enter directory in which to create and download data
read -p "Enter path for downloading raw data: " fullname
echo $fullname

# Create directories in the designated path and download the files into these
mkdir "${fullname}/GSE133556/"
wget -P "${fullname}/GSE133556/" 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE133nnn/GSE133556/suppl/GSE133556_RAW.tar'

mkdir "${fullname}/GSE155760/"
wget -P "${fullname}/GSE155760/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE155nnn/GSE155760/suppl/GSE155760_RAW.tar'

mkdir "${fullname}/GSE185008/"
wget -P "${fullname}/GSE185008/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185008/suppl/GSE185008_RAW.tar'

mkdir "${fullname}/GSE211686/"
wget -P "${fullname}/GSE211686/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE211nnn/GSE211686/suppl/GSE211686_RAW.tar'

mkdir "${fullname}/GSE263434/"
wget -P "${fullname}/GSE263434/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263434/suppl/GSE263434_RAW.tar'

mkdir "${fullname}/GSE226872/"
wget -P "${fullname}/GSE226872/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226872/suppl/GSE226872_RAW.tar'

mkdir "${fullname}/GSE267068/"
wget -P "${fullname}/GSE267068/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267068/suppl/GSE267068_RAW.tar'

mkdir "${fullname}/GSE51820/"
wget -P "${fullname}/GSE51820/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51820/suppl/GSE51820_non_normalized.txt.gz'

mkdir  "${fullname}/GSE226823/"
wget -P "${fullname}/GSE226823/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226823/suppl/GSE226823_RAW.tar'

mkdir "${fullname}/GSE263434/"
wget -P "${fullname}/GSE263434/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263434/suppl/GSE263434_RAW.tar'

mkdir "${fullname}/Train/"
wget -P "${fullname}/Train/" 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE326nnn/GSE326000/suppl/GSE326000_RAW.tar'

# Untar all the raw data files
cd $fullname
find . -type f -iname "*.tar" -print0 -execdir tar xf {} \; -delete