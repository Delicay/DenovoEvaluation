# DenovoEvaluation

## Overview
Using high-density genetic anchor (4.4M) to evaluate de novo genome assembly of maize. If genetic anchors are available, it also can be used in other species.

## Snapshot

Nearly 5,000 genetic anchors are used to validate the quality of a 10Mb scaffold from de novo genome assembly 

![validation](https://c1.staticflickr.com/5/4204/34901070415_5fd5902466_k.jpg)

## Prerequisites

Java 8

http://www.oracle.com/technetwork/java/javase/overview/java8-2100321.html

R

https://cran.r-project.org/

## Release

Currently at v1.0.1

## Usage

Please keep the lib folder in the working directory. From Linux command line, type:

java -Xms2g -Xmx16g -jar DenovoEvaluation.jar parameter.txt > log.txt &

The user-custom parameter file is in the released java package. More details of usage are in the parameter file. Results will be generated in the result directory.

## Contributor

Fei Lu

flu@genetics.ac.cn; dr.lufei@gmail.com

https://sites.google.com/site/feilu0000/

https://plantgeneticslab.weebly.com

## Citation

[Fei Lu, Maria C Romay, Jeff C Glaubitz, Peter J Bradbury, Rob J Elshire, Tianyu Wang, Yu Li, Yongxiang Li, Kassa Semagn, Xuecai Zhang, Alvaro G. Hernandez, Mark A. Mikel, Ilya Soifer, Omer Barad, Edward S Buckler. (2015) High-resolution genetic mapping of maize pan-genome sequence anchors. Nature Communications. DOI:10.1038/ncomms7914!](https://www.nature.com/articles/ncomms7914)
