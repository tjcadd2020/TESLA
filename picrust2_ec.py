import os
from os import listdir
import sys,getopt
from ete3 import Tree
#输入两个参数，第一个为待处理数据，第二个为该进程所用核数
a=getopt.getopt(sys.argv[1:],'-i-n', ['input','cores'])[1][0]
b=getopt.getopt(sys.argv[1:],'-i-n', ['input','cores'])[1][2]
path=os.path.abspath(sys.argv[0]).rsplit('/',1)[0]+'/'
name=a.rsplit('.',1)[0]
out_dir=name+'_temp/'
#建立中间文件所用的文件夹
def makedir():
    if os.path.exists(out_dir+'/ec') is False:
       print(os.system('mkdir '+out_dir+'ec'))
    if os.path.exists(out_dir+'ec/ref_id') is False:
       print(os.system('mkdir '+out_dir+'ec/ref_id'))
    if os.path.exists(out_dir+'ec/ref_16s') is False:
       print(os.system('mkdir '+out_dir+'ec/ref_16s'))
    if os.path.exists(out_dir+'ec/ref_ec') is False:
       print(os.system('mkdir '+out_dir+'ec/ref_ec'))
    if os.path.exists(out_dir+'ec/ref_table') is False:
       print(os.system('mkdir '+out_dir+'ec/ref_table'))
    if os.path.exists(out_dir+'ec/marker') is False:
       print(os.system('mkdir '+out_dir+'ec/marker'))
    if os.path.exists(out_dir+'ec/ec') is False:
       print(os.system('mkdir '+out_dir+'ec/ec'))
    if os.path.exists(out_dir+'ec/ec') is False:
       print(os.system('mkdir '+out_dir+'ec/ec'))
    if os.path.exists(out_dir+'ec/ref_ec') is False:
       print(os.system('mkdir '+out_dir+'ec/ref_ec'))
    if os.path.exists(out_dir+'ec/ref_table') is False:
       print(os.system('mkdir '+out_dir+'ec/ref_table'))
    if os.path.exists(out_dir+'ec/meta_ec') is False:
       print(os.system('mkdir '+out_dir+'ec/meta_ec'))
makedir()

query_id=[]

tree_dic={}

classify_in=out_dir+'query/'+listdir(out_dir+'query')[0]+'/data/'+listdir(out_dir+'query/'+listdir(out_dir+'query')[0]+'/data')[0]

query_in=[]

for line in open(classify_in):
    query_in.append(line.strip('\n').strip('>'))

def tax(sheet_in):
    sheet=open(sheet_in).readlines()
    for i in range(1,len(sheet)):
        tax=sheet[i].split('\t')[1].strip('\n')
        ids=sheet[i].split('\t')[0].strip('>')
        if 's__'  in tax and ';' in tax:
           gen=sheet[i].split('\t')[1].split(';')[5]
           tree_dic.setdefault(ids,gen)
           with open(out_dir+'seq/'+ids,'w') as f:
                f.write('>'+ids+'\n')
                f.write(query_in[query_in.index(ids)+1])
    return None
#读取TESLA初步注释得到的分类文件，将其分为注释到种以及未注释到种两部分
sheet_loc=out_dir+'classified/'+listdir(out_dir+'classified')[0]+'/data/'+listdir(out_dir+'classified/'+listdir(out_dir+'classified')[0]+'/data')[0]

tax(sheet_loc)

#注释到种的被鉴定序列由于TESLA二次注释未能进行系统发育分析，Picrust2功能步骤中需将其插入到对应属水平进化树中

def tree_place(tax):
   global tree_dic
   global out_dir
   for key in tree_dic.keys():
      if key+'.tre' not in listdir(out_dir+'tree'):
         place_tree=os.system('place_seqs.py -s '+out_dir+'seq/'+key+' -o '+out_dir+'tree/'+key+'.tre'+' --ref_dir '+path+'data/genus/'+tree_dic[key]+' -p '+b)
         print(place_tree)

tree_place(tax)

for i in range(len(listdir(out_dir+'tree'))):
    if '.tre' in listdir(out_dir+'tree')[i]:
       query_id.append(listdir(out_dir+'tree')[i].replace('.tre',''))

#读取参考序列信息
def tree_id(query_id):
    for i in range(len(query_id)):
        id_show=os.system('python '+path+'tools/node.py '+out_dir+'tree/'+query_id[i]+'.tre '+out_dir+'ec/ref_id')
        print(id_show)

tree_id(query_id)   

def search_query(txt,suf1,suf2,loc1,loc2):
    global query_id
    txt_id=[]
    for j in range(len(txt)):
        txt_id.append(txt[j].split('\t')[0])
    for i in range(len(query_id)):
        with open(loc2+query_id[i]+suf2,'w') as f:
             f.write(txt[0])
             for line in open(loc1+query_id[i]+suf1):
                 line=line.strip('\n')
                 if line in txt_id:             
                    f.write(txt[txt_id.index(line)])
#导入参考序列相关数据                
s16=open(path+'tools/16S.txt','r').readlines()
ec_assembly=open(path+'tools/ec_db.txt','r').readlines()
search_query(s16,'_id','_16s',out_dir+'ec/ref_id/',out_dir+'ec/ref_16s/')
search_query(ec_assembly,'_id','_ec',out_dir+'ec/ref_id/',out_dir+'ec/ref_ec/')

def table_make(query):
    global out_dir
    for i in range(len(query)):
        with open(out_dir+'ec/ref_table/'+query[i]+'_table.txt','w') as f:
             f.write('#OTU_num'+'\t'+'c1'+'\n')
             f.write(query[i]+'\t'+'1')

table_make(query_id)
#隐藏状态预测，得到基因家族丰度
def hsp1(query_id):
    for i in range(len(query_id)):
        hsp_1=os.system('hsp.py --tree '+out_dir+'tree/'+query_id[i]+'.tre --output '+out_dir+'ec/marker/'+query_id[i]+'_marker.tsv.gz --observed_trait_table '+out_dir+'ec/ref_16s/'+query_id[i]+'_16s')
        print(hsp_1)

def hsp2(query_id):
    for i in range(len(query_id)):
        hsp_2=os.system('hsp.py --tree '+out_dir+'tree/'+query_id[i]+'.tre --output '+out_dir+'ec/ec/EC_'+query_id[i]+'.tsv.gz --observed_trait_table '+out_dir+'ec/ref_ec/'+query_id[i]+'_ec')
        print(hsp_2)
#确定每个样本的基因家族丰度
def metapipe(query_id):
    for i in range(len(query_id)):
        pipeline=os.system('metagenome_pipeline.py -i '+out_dir+'ec/ref_table/'+query_id[i]+'_table.txt -m '+out_dir+'ec/marker/'+query_id[i]+'_marker.tsv.gz -f '+out_dir+'ec/ec/EC_'+query_id[i]+'.tsv.gz -o '+out_dir+'ec/meta_ec/'+query_id[i]+'_meta_out')
        print(pipeline)
#根据通路映射文件推断通路丰度
def pathpipe(query_id):
    for i in range(len(query_id)):
        pipeline=os.system('pathway_pipeline.py -i '+out_dir+'ec/meta_ec/'+query_id[i]+'_meta_out/pred_metagenome_unstrat.tsv.gz -o '+out_dir+query_id[i]+'_pathways_out --map '+path+'tools/metacyc_path2rxn_struc_filt_pro.txt')
        print(pipeline)

hsp1(query_id)
hsp2(query_id)

metapipe(query_id)
pathpipe(query_id)
