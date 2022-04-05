import os
from os import listdir
import sys,getopt
from ete3 import Tree

#输入两个参数，第一个为待处理数据，第二个为该进程所用核数
a=getopt.getopt(sys.argv[1:],'-i-n', ['input','cores'])[1][0]
b=getopt.getopt(sys.argv[1:],'-i-n', ['input','cores'])[1][2]
path=os.path.abspath(sys.argv[0]).rsplit('/',1)[0]+'/'
name=a.rsplit('.',1)[0]
classified=[]
unclassified=[]

notax={}
havespe={}
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

#解压TESLA初步注释所得的分类文件，之后读取其分类情况等相关信息
def unzip(a):
    zip_in='suptax_'+a
    zip_qu='dada2_'+a
    name=a.rsplit('.',1)[0]
    unzip_pipeline=os.system('y|unzip '+out_dir+zip_in+' -d '+out_dir+'classified')
    unzip_pipeline_qu=os.system('y|unzip '+out_dir+zip_qu+' -d '+out_dir+'query')
    return(unzip_pipeline,unzip_pipeline_qu)

if os.path.exists(out_dir+'suptax_'+a) is True and os.path.exists(out_dir+'dada2_'+a) is True and len(listdir(out_dir+'classified'))==0 and len(listdir(out_dir+'query'))==0:
   unzip(a)

def chose_tax(a):
    while a.rsplit(';',1)[1][1:]=='__':
          a=a.rsplit(';',1)[0]
    return(a.rsplit(';',1)[1])
#获取初步注释分类结果中注释到种水平的序列和未注释到种水平序列的信息
def tax(sheet_in):
    sheet=open(sheet_in).readlines()
    global classified
    global unclassified
    global notax
    global havespe
    for i in range(1,len(sheet)):
        tax=sheet[i].split('\t')[1].strip('\n')
        ids=sheet[i].split('\t')[0].strip('>')
        if 's__' not in tax and ';' in tax:
            notax[ids]=chose_tax(tax)
            unclassified.append(tax)
        if 's__' in tax:
            classified.append(tax)
            havespe[ids]=sheet[i].split('\t')[1].strip('\n')
    return None

sheet_loc=out_dir+'classified/'+listdir(out_dir+'classified')[0]+'/data/'+listdir(out_dir+'classified/'+listdir(out_dir+'classified')[0]+'/data')[0]

tax(sheet_loc)
#输出未注释到种水平的分类
if os.path.exists(name+'_unclassified.txt') is False:
   with open(name+'_unclassified.txt','w') as f:
        for i in range(len(unclassified)):
            f.write(unclassified[i]+'\n')
#输出注释到种水平的分类
if os.path.exists(name+'_classified.txt') is False:
   with open(name+'_classified.txt','w') as f:
        for i in range(len(classified)):
            f.write(classified[i]+'\n')

classify_in=out_dir+'query/'+listdir(out_dir+'query')[0]+'/data/'+listdir(out_dir+'query/'+listdir(out_dir+'query')[0]+'/data')[0]

query_in=[]

for line in open(classify_in):
    query_in.append(line.strip('\n').strip('>'))

for key in notax.keys():
    with open(out_dir+'seq/'+key,'w') as f:
         f.write('>'+key+'\n')
         f.write(query_in[query_in.index(key)+1])

#确定未注释到种水平序列的最近分类，之后将未注释到种水平序列插入到对应进化树中
def tree_place(notax):
   global out_dir
   for key in notax.keys():
      if notax[key][0]=='p':
         loc='phylum'
      if notax[key][0]=='c':
         loc='class'
      if notax[key][0]=='o':
         loc='order'
      if notax[key][0]=='f':
         loc='family'
      if notax[key][0]=='g':
         loc='genus'
      place_tree=os.system('place_seqs.py -s '+out_dir+'seq/'+key+' -o '+out_dir+'tree/'+key+'.tre'+' --ref_dir '+path+'data/'+loc+'/'+notax[key]+' -p '+b)
      print(place_tree)

if len(listdir(out_dir+'tree'))==0:
   tree_place(notax)

tree_in=listdir(out_dir+'tree')

#计算未注释到种的被鉴定序列在进化树中的邻近序列
def tree_distance(tree):
    query=tree.rsplit('/',1)[1].split('.')[0]
    t=Tree(tree)
    dis=[]
    no=[]
    for node in t:
        if node.name != query:
           no.append(node.name)
    for i in range(len(no)):
        dis.append(t.get_distance(query,no[i]))
    return(no[dis.index(min(dis))])

query={}

for i in range(len(tree_in)):
    tree_id=tree_in[i].split('.')[0]
    query.setdefault(tree_id,[])
    query[tree_id].append(tree_distance(out_dir+'tree/'+tree_id+'.tre'))
    query[tree_id].append(notax[tree_id])

reference_seq=[]
for line in open(path+'tools/reference.fna'):
    reference_seq.append(line.strip('\n').strip('>'))

def ref_seq(dic):
    global reference_seq
    for key in dic.keys():
        with open(out_dir+'refer/'+dic[key][0]+'.fna','w') as f:
             f.write('>'+dic[key][0]+'\n')
             f.write(reference_seq[reference_seq.index(dic[key][0])+1]+'\n')

ref_seq(query)
#计算被鉴定序列与进化树邻近序列的mash距离
def mashdis(query):
   for key in query.keys():
      query_in=out_dir+'seq/'+key
      refer=out_dir+'refer/'+query[key][0]+'.fna'
      mash=os.system(path+'tools/mash dist '+query_in+' '+refer+'>>'+out_dir+'mash_result')
      print(mash)
   return

if os.path.exists(out_dir+'/mash_result') is False:
   mashdis(query)

if os.path.exists(out_dir+'/mash_result') is True:
   for line in open(out_dir+'/mash_result'):
       query[line.split('\t')[0].rsplit('/',1)[1]].append(float(line.split('\t')[2]))

for line in open(sheet_loc):
    if line.split('\t')[0] in query.keys():
       query[line.split('\t')[0]].append(line.split('\t')[1].strip('\n'))

tax_id=[]
tax_ref=[]

for line in open(path+'tools/taxonomy.txt'):
    tax_id.append(line.split('\t')[0])
    tax_ref.append(line.split('\t')[1].strip('\n')) 

for key in query.keys():
    query[key].append(tax_ref[tax_id.index(query[key][0])]) 
#根据被鉴定序列与其邻近序列的mash距离以及邻近序列的分类为其赋予新命名
def add_name(tax,n,name):
    tax_list=[]
    if tax.count(';')==6:
       tax_l=tax.split(';')
       for i in range(len(tax_l)):
           tax_list.append(tax_l[i][3:])
       new_tax=[]
       for i in range(0,n):
           new_tax.append(tax_list[i])
       for i in range(n,len(tax_list)):
           new_tax.append(name+tax_list[i])
       seq='k__'+new_tax[0]+';'+'p__'+new_tax[1]+';'+'c__'+new_tax[2]+';'+'o__'+new_tax[3]+';'+'f__'+new_tax[4]+';'+'g__'+new_tax[5]+';'+'s__'+new_tax[6]
       return seq

def addtax(query):
    for key in query.keys():
        dis=query[key][2]
        if dis >=0 and dis <0.002:
           pre='iso_'
           query[key].append(1)
        if dis >=0.002 and dis <0.005:
           pre='simili_'
           query[key].append(2) 
        if dis >=0.005 and dis <0.035:
           pre='iuxta_'
           query[key].append(3)
        if dis >=0.035 and dis <0.06:
           pre='cognati_'
           query[key].append(4)
        if dis >=0.06 and dis <0.08:
           pre='anomeo_'
           query[key].append(5)
        if dis >=0.08:
           pre='tele_'
           query[key].append(6)

    for key in query.keys():
        if query[key][1][0]=='g' and query[key][5]==2:
           query[key].append(add_name(query[key][4],6,pre))
        if query[key][1][0]=='g' and query[key][5]<2:
           query[key].append(add_name(query[key][4],6,'extra-'+pre))
        if query[key][1][0]=='g' and query[key][5]>2:
           query[key].append(add_name(query[key][4],6,'meta-'+pre))
#gen
        if query[key][1][0]=='g' and int(query[key][5])==2:
           query[key].append(add_name(query[key][4],6,pre))
        if query[key][1][0]=='g' and int(query[key][5])<2:
           query[key].append(add_name(query[key][4],6,'extra-'+pre))
        if query[key][1][0]=='g' and int(query[key][5])>2:
           query[key].append(add_name(query[key][4],6,'meta-'+pre))
#fam
        if query[key][1][0]=='f' and query[key][5]==3:
           query[key].append(add_name(query[key][4],5,pre))
        if query[key][1][0]=='f' and query[key][5]<3:
           query[key].append(add_name(query[key][4],5,'extra-'+pre))
        if query[key][1][0]=='f' and query[key][5]>3:
           query[key].append(add_name(query[key][4],5,'meta-'+pre))
#ord
        if query[key][1][0]=='o' and query[key][5]==4:
           query[key].append(add_name(query[key][4],4,pre))
        if query[key][1][0]=='o' and query[key][5]<4:
           query[key].append(add_name(query[key][4],4,'extra-'+pre))
        if query[key][1][0]=='o' and query[key][5]>4:
           query[key].append(add_name(query[key][4],4,'meta-'+pre))
#cla
        if query[key][1][0]=='c' and query[key][5]==5:
           query[key].append(add_name(query[key][4],3,pre))
        if query[key][1][0]=='c' and query[key][5]<5:
           query[key].append(add_name(query[key][4],3,'extra-'+pre))
        if query[key][1][0]=='c' and query[key][5]>5:
           query[key].append(add_name(query[key][4],3,'meta-'+pre))
#phy
        if query[key][1][0]=='p' and query[key][5]==6:
           query[key].append(add_name(query[key][4],2,pre))
        if query[key][1][0]=='p' and query[key][5]<6:
           query[key].append(add_name(query[key][4],2,'extra-'+pre))
        if query[key][1][0]=='p' and query[key][5]>6:
           query[key].append(add_name(query[key][4],2,'meta-'+pre))
    return(query)

addtax(query)
#输出所有被鉴定序列的id以及其分类，为了方便相关流程，其中初步注释未注释到种的序列分类为其进化树邻近序列分类
if os.path.exists(name+'_allclassified.txt') is False:
   with open(name+'_allclassified.txt','w') as f:
        for key in havespe.keys():
            f.write(key+'\t'+havespe[key]+'\n')
        for key in query.keys():
            f.write(key+'\t'+query[key][4]+'\n')

#输出初步注释未能注释种水平，经过二次注释赋予新命名的序列分类
if os.path.exists(name+'_twiceclassified.txt') is False:
   with open(name+'_twiceclassified.txt','w') as f:
        for key in query.keys():
            f.write(query[key][6]+'\n')
