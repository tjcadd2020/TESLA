TESLA初步注释代码:
import os
from os import listdir
import sys,getopt
from ete3 import Tree

#需要两个输入数据，第一个为待处理数据，第二个为该进程所用核数
a=getopt.getopt(sys.argv[1:],'-i-n', ['input','cores'])[1][0]
b=getopt.getopt(sys.argv[1:],'-i-n', ['input','cores'])[1][2]

path=os.path.abspath(sys.argv[0]).rsplit('/',1)[0]+'/'
name=a.rsplit('.',1)[0]
classified=[]
unclassified=[]
notax={}
out_dir=name+'_temp/'

#建立中间文件所用的文件夹
def makedir():
    if os.path.exists(out_dir) is False:
       print(os.system('mkdir '+out_dir))
    if os.path.exists(out_dir+'/classified') is False:
       print(os.system('mkdir '+out_dir+'classified'))
    if os.path.exists(out_dir+'/query') is False:
       print(os.system('mkdir '+out_dir+'query'))
    if os.path.exists(out_dir+'/seq') is False:
       print(os.system('mkdir '+out_dir+'seq'))
    if os.path.exists(out_dir+'/tree') is False:
       print(os.system('mkdir '+out_dir+'tree'))
    if os.path.exists(out_dir+'/refer') is False:
       print(os.system('mkdir '+out_dir+'refer'))
    if os.path.exists(out_dir+'/mash') is False:
       print(os.system('mkdir '+out_dir+'mash'))

makedir()

#对待处理数据进行剪切，提取其V4区
def trim(a,b):
    trim_out='trim_'+a
    trim_pipline=os.system('~/software/anaconda3/envs/qiime2-2021.4/bin/qiime cutadapt trim-paired --i-demultiplexed-sequences '+a+' --p-front-f GTGCCAGCMGCCGCGGTAA --p-front-r GGACTACHVGGGTWTCTAAT --o-trimmed-sequences '+out_dir+trim_out+' --p-cores '+b+' --verbose')
    return(trim_pipline)

#使用qiime2中的dada2模块降噪形成asv序列
def denoise(a,b):
    trim_out='trim_'+a
    denoise_out='dada2_'+a
    table_out='table_'+a
    stats_out='stat_'+a
    denoise_pipeline=os.system('~/software/anaconda3/envs/qiime2-2021.4/bin/qiime dada2 denoise-paired --i-demultiplexed-seqs '+out_dir+trim_out+' --p-trim-left-f 0 --p-trunc-len-f 200 --p-trim-left-r 0 --p-trunc-len-r 200 --p-n-threads '+b+' --o-table '+out_dir+table_out+' --o-representative-sequences '+out_dir+denoise_out+' --o-denoising-stats '+out_dir+stats_out)
    return(denoise_pipeline)

#使用课题构建的TESLA初步注释分类器进行分类
def classify(a):
    classify_in='dada2_'+a
    classify_out='suptax_'+a
    classify_pipeline=os.system('~/software/anaconda3/envs/qiime2-2021.4/bin/qiime feature-classifier classify-sklearn --i-classifier '+path+ 'tools/sup_classifier.qza '+'--i-reads '+out_dir+classify_in+' --o-classification '+out_dir+classify_out+' --p-confidence 0.5')
    return(classify_pipeline)

if os.path.exists(out_dir+'trim_'+a) is False:
   trim(a,b)
if os.path.exists(out_dir+'dada2_'+a) is False or os.path.exists(out_dir+'stat_'+a) is False: 
   denoise(a,b)
if os.path.exists(out_dir+'suptax_'+a) is False:
   classify(a)
