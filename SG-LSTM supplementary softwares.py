
import sys
from PyQt5 import QtWidgets,uic
import pandas as pd
import time
import numpy as np
import datetime
from PyQt5.QtWidgets import QWidget, QApplication, QPushButton, QColorDialog, QFontDialog, QTextEdit, QFileDialog
#from sklearn.externals import joblib#引入joblib用
import joblib
import pandas as pd
qtCreatorFile = "mainwindow.ui" # Enter file here.
 
Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
#获取数字模板
#  try:
#          index_number = int(self.index.toPlainText())
#        except:
#          index_number=1
class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        global string_history
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        string_history=''
        self.setWindowTitle("SG-LSTM supplement sofware,version 1.0.1")
        self.load_mirna.clicked.connect(self.Load_mirna)
        self.load_gene.clicked.connect(self.Load_gene)
        #self.generate_s_predict.clicked.connect(self.Generate_S)
        #self.generate_P.clicked.connect(self.Generate_P)
        self.load_real_data.clicked.connect(self.Load_real_data)
        #self.case_study_button.clicked.connect(self.Case_study)
        self.G_pos.clicked.connect(self.gene_posi)
        self.G_nega.clicked.connect(self.gene_manu_nega)
        self.gene_csv.clicked.connect(self.Gene_csv)
        self.load_posi.clicked.connect(self.Load_posi)
        self.nega_file.clicked.connect(self.load_nega_file)
        self.load_nega.clicked.connect(self.Load_nega)
        #self.test_R.clicked.connect(self.test_retained_test)
        #self.Feature1.clicked.connect(self.Feature_load1)
        #self.Feature2.clicked.connect(self.Feature_load2)
        #self.test_mode.clicked.connect(self.Test_mode)
        #self.newdataset.clicked.connect(self.new_dataset)
        self.calculate_distance.clicked.connect(self.CDA)
    
    def CDA(self):
        global mirna_file
        global mirna_label
        global gene_file
        global gene_label
        
        def calculate_euclidean(a,b):
              from scipy.spatial import distance
              dst = distance.euclidean(a,b)
              return dst
        def cos_sim(vector_a, vector_b):
             vector_a = np.mat(vector_a)
             vector_b = np.mat(vector_b)
             num = float(vector_a * vector_b.T)
             denom = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
             cos = num / denom
             sim = 0.5 * cos
             return sim
        global string_history
        try:
          print(len(mirna_file))
          print(len(gene_file))
        except:
          string_history=string_history+'gene file/ mirna file unavailable'+'\n'
          self.state_show.setText(string_history)
          return self
        try:
          global positive_mean
          print(positive_mean)
          global positive_list
          print('rows of PS',len(postive_list))
        except:
          string_history=string_history+'please load postive samples file first'+'\n'
          self.state_show.setText(string_history)
          return self
        th1=len(mirna_file)
        th2=len(gene_file)
        cos_list=[]
        eu_list=[]
        for i in range(0,th1):
            print('calculating {} mirna'.format(i))
            for j in range(0,th2):
                key={'mirna':mirna_label[i],'gene':gene_label[j]}
                if key not in postive_list:
                    feature1=mirna_file[i:i+1]
                    feature1=feature1.reset_index(drop=True)
                    
                    feature2=gene_file[j:j+1]
                    feature2=feature2.reset_index(drop=True)
                    temp=pd.concat([feature1,feature2],axis=1,join='outer')
                    #print((temp))
                    temp_cos=cos_sim(positive_mean,temp)
                    temp_dis=calculate_euclidean(positive_mean,temp)
                    if(j==th2-1):
                      print(temp_cos,temp_dis)
                    cos_list.append(temp_cos)
                    eu_list.append(temp_dis)
            np.save('cos_list',cos_list)
            np.save('eu_list',eu_list)
            self.state_show.setText('calculated distance successfully')
    def new_dataset(self):
        
        def check_out(index_name,index_gene):
             global mirna_label
             global gene_label
             global list_temp
             key={'miRNA':mirna_label[index_name],'gene':gene_label[index_gene]}
             if key in list_temp:
                 return 1
             else:
                 return 0
        def find_mirna(mirna):#找到miRNA返回其index
          global mirna_label
          for i in range(0,len(mirna_label)):
            if mirna_label[i]==mirna:
                return i
        def find_gene(gene):
          global gene_label
          for i in range(0,len(gene_label)):
            if gene_label[i]==gene:
                return i
        #try:
        def ab():
          global string_history
          fname = QFileDialog.getOpenFileName(self, '打开文件','./', "Text Files (*.csv)")
          file=pd.read_csv('{}'.format(fname[0]))
          S_mirna=file['mirna']
          S_gene=file['gene']
          string_history=string_history+'文件读取成功{}'.format(fname[0])+'\n'
          self.state_show.setText(string_history)#结果显示
          skip=0
          count=0
          list_label=[]
          print(len(S_mirna))
          print(len(S_gene))
          for i in range(0,3):
             print('{}个miRNA新数据集生成完毕'.format(i))
             for j in range(0,len(S_gene)):
              
              try:   
               index_name=find_mirna(S_mirna[i])
               feature1=mirna_file[index_name:index_name+1]
               feature1=feature1.reset_index(drop=True)
               index_gene=find_gene(S_gene[j])
               feature2=gene_file[index_gene:index_gene+1]
               feature2=feature2.reset_index(drop=True)
               if check_out(index_name,index_gene)==1:
                 list_label.append({'0_mirna':S_mirna[i],'1_gene':S_gene[j],'label':1})
               else:
                 list_label.append({'0_mirna':S_mirna[i],'1_gene':S_gene[j],'label':0})
               if count==0:
                  temp0=pd.concat([feature1,feature2],axis=1,join='outer')
                  count+=1
               elif count==1:                  
                  temp1=pd.concat([feature1,feature2],axis=1,join='outer')
                  feature=pd.concat([temp0,temp1],axis=0)
                  count+=1
               else:
                  temp1=pd.concat([feature1,feature2],axis=1,join='outer')
                  feature=pd.concat([feature,temp1],axis=0)
               # 说明二者是有关系的 label应该是1 否则二者没有关系 标签为0
              except:
                  skip+=1
          fea_label=pd.DataFrame(list_label)
          feature=feature.reset_index(drop=True)
         #positive_mean=feature.mean()
          feature=pd.concat([feature,fea_label],axis=1)#这个feature 就包含了所有的抽取的正样本
         #print(len(feature),len(fea_label),feature,fea_label)#mirna+gene的嵌入合在一起的
         #string_result=string_result+'正样本生成完毕'
         #self.result_show.setText(string_result)#结果显示
          string_history=string_history+'读取文件结束,新数据集new_dataset.csv生成完毕'+'\n'
          self.state_show.setText(string_history)#结果显示
          feature.to_csv('new_dataset.csv',index=None)
        ab()
        #except:
          #result_string='输入截取miRNA和gene 从而生成新的样本集，文件应该有两列 分别是miRNA以及基因'
          #self.result_show.setText(result_string)
                  
          
 
    def load_nega_file(self):
        global feature_negative
        global negative_mean
        global string_history
        result_string=''
        try:
          fname = QFileDialog.getOpenFileName(self, '打开文件','./', "Text Files (*.csv)")
          file=pd.read_csv('{}'.format(fname[0]))
          feature_negative=file
          file=file.drop('0_mirna',axis=1)
          file=file.drop('1_gene',axis=1)
          file=file.drop('label',axis=1)
          negative_mean=file.mean()
          string_history=string_history+'NS file{}loaded'.format(fname[0])+'\n'
          self.state_show.setText(string_history)#结果显示
        except:
          result_string='pls load already made NS file,columns should be /1_gene,label/respectively'
          self.result_show.setText(result_string)
    def Load_nega(self):#使用文件生成负样本
         global string_history
         global gene_file
         global mirna_file
         global feature#这里是正样本
         global positive_mean
         global feature_negative
         global negative_mean
         def cos_sim(vector_a, vector_b):
             vector_a = np.mat(vector_a)
             vector_b = np.mat(vector_b)
             num = float(vector_a * vector_b.T)
             denom = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
             cos = num / denom
             sim = 0.5 * cos
             return sim
               
         def calculate_cos(feature):
             global positive_mean #这里需要一个正样本的mean_值
             return cos_sim(positive_mean,feature)
         def calculate_euclidean(feature):
              from scipy.spatial import distance
              dst = distance.euclidean(positive_mean,feature)
              return dst
             #看看欧式距离看看    
         def find_mirna(mirna):
          global mirna_label
          for i in range(0,len(mirna_label)):
            if mirna_label[i]==mirna:
                return i
         def find_gene(gene):
          global gene_label
          for i in range(0,len(gene_label)):
            if gene_label[i]==gene:
                return i
            
         fname = QFileDialog.getOpenFileName(self, '打开文件','./', "Text Files (*.csv)")
         file=pd.read_csv('{}'.format(fname[0]))
         
         S_mirna=file['miRNA']
         S_gene=file['Target Gene']
         print('read{}'.format(fname[0]))
         
         list_label=[]
         count=0
         skip=0
         string_result=''
         list_distance=[]
         for i in range(0,len(file)):
            if i%1000==0:
                print(i)
                string_result=string_result+'{}finished'.format(i)
                self.result_show.setText(string_result)#结果显示
            try: #如果找到某一对才做
             
             index_name=find_mirna(S_mirna[i])
             feature1=mirna_file[index_name:index_name+1]
             feature1=feature1.reset_index(drop=True)
             index_gene=find_gene(S_gene[i])
             feature2=gene_file[index_gene:index_gene+1]
             feature2=feature2.reset_index(drop=True)
             temp=pd.concat([feature1,feature2],axis=1,join='outer')
             cor=calculate_cos(temp)
             eu=calculate_euclidean(temp)
             list_label.append({'0_mirna':S_mirna[i],'1_gene':S_gene[i],'label':0})#'cor':cor,'grount_truth':1})
             #print('余弦相似度',calculate_cos(temp))
             #print('和正样本的欧几里得距离',eu)
             var=calculate_euclidean(temp)
             list_distance.append(var)
             
             if count==0:
                temp0=pd.concat([feature1,feature2],axis=1,join='outer')
                count+=1
             elif count==1:                  
                temp1=pd.concat([feature1,feature2],axis=1,join='outer')
                feature_negative=pd.concat([temp0,temp1],axis=0)
                count+=1
             else:
                temp1=pd.concat([feature1,feature2],axis=1,join='outer')
                feature_negative=pd.concat([feature_negative,temp1],axis=0)
                #print(index_name,index_gene,feature)#mirna+gene的嵌入合在一起的
            
            except:
                skip+=1
                print('skip','mirna',S_mirna[i],'gene',S_gene[i])
         fea_label=pd.DataFrame(list_label)
         feature_negative=feature_negative.reset_index(drop=True)
         negative_mean=feature_negative.mean()
         feature_negative=pd.concat([feature_negative,fea_label],axis=1)#这个feature 就包含了所有的抽取的正样本
        ## print(len(feature_negative),len(fea_label),feature_negative,fea_label)#mirna+gene的嵌入合在一起的
         #print('和正样本的平均欧几里得距离',(sum(list_distance)/len(list_distance)))#4.684385237024212
         #string_result=string_result+'正样本生成完毕'
         #self.result_show.setText(string_result)#结果显示
         string_history=string_history+'file read{}finished,negative_sample generation done'.format(fname[0])+'\n'
         self.state_show.setText(string_history)#结果显示
         feature_negative.to_csv('negative_sample.csv',index=None)
        
    def Load_posi(self):
        global feature
        global positive_mean
        global string_history
        global postive_list
        result_string=''
        try:
          fname = QFileDialog.getOpenFileName(self, '打开文件','./', "Text Files (*.csv)")
          file=pd.read_csv('{}'.format(fname[0]))
          feature=file
          postive_list=[]
          mirna=file['0_mirna']
          gene=file['1_gene']
          for i in range(0,len(file)):
                 postive_list.append({'mirna':mirna[i],'gene':gene[i]})
          print('postive list generation complete')
          file=file.drop('0_mirna',axis=1)
          file=file.drop('1_gene',axis=1)
          file=file.drop('label',axis=1)
          positive_mean=file.mean()
          string_history=string_history+'positive file{}loaded，rows{},finished'.format(fname[0],len(file))+'\n'
          self.state_show.setText(string_history)#结果显示
        except:
          result_string='pls load already made postive samples ，columns should be /0_mirna,1_gene,label/respectively'
          self.result_show.setText(result_string)
    
    def Gene_csv(self):
        global feature_negative
        global feature
        global mirna_name
        def alter(mirna_name):
            i=mirna_name.find('.csv')
            return mirna_name[0:i]
        mirna_name=alter(mirna_name)
        
        feature=feature.reset_index(drop=True)
        feature_negative=feature_negative.reset_index(drop=True)
        #print(feature)
        #print(feature.columns)
        #print(feature_negative.columns)
        #print(feature_negative)
        sum_up=pd.concat([feature,feature_negative],axis=0)#等于0是加到后面
        #print(sum_up)
        sum_up=sum_up.reset_index(drop=True)
        
        try:
            plain_name=str(self.name.toPlainText())
            if len(plain_name)==0:
               name='{}_default.csv'.format(mirna_name)
            else:
               name='{}.csv'.format(plain_name) 
        except:
            name='{}_default.csv'.format(mirna_name)
        sum_up.to_csv('{}'.format(name),index=None)
        global string_history
        string_history=string_history+'file // {}  // generated'.format(name)+'\n'
        self.state_show.setText(string_history)#结果显示
    def gene_manu_nega(self):#生成人工的数据集 默认数目200条
         global string_history
         global gene_file
         global mirna_file
         global mirna_label
         global gene_label
         global feature_ne
         global string_result
         global feature_negative #需要通过这个负的特征 来最终得到所有的特征
         #global negative_mean
         global positive_mean
         count=0;
         list_label=[]
         early_stop=0
         try:
          array_number = int(self.NS_number.toPlainText())
         except:
          array_number=29725
         try:
          th_cosin=float(self.th_cos.value())
          print('cos threshold received',th_cosin)
         except:
          th_cosin=0.2
         try:
          th_distance=float(self.th_dis.value())
          print('distance threshold received',th_distance)
         except:
          th_distance=5
         def check_out(index_name,index_gene):#检查是否是在一开始的东西中就有的
             global mirna_label
             global gene_label
             global list_temp
             key={'miRNA':mirna_label[index_name],'gene':gene_label[index_gene]}
             if key in list_temp:
                 return 1
             else:
                 return 0
         def calculate_euclidean(a,b):
              from scipy.spatial import distance
              dst = distance.euclidean(a,b)
              return dst
         def cos_sim(vector_a, vector_b):
             vector_a = np.mat(vector_a)
             vector_b = np.mat(vector_b)
             num = float(vector_a * vector_b.T)
             denom = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
             cos = num / denom
             sim = 0.5 * cos
             return sim

         def calculate_cos(feature):
             global positive_mean #这里需要一个正样本的mean_值
             return cos_sim(positive_mean,feature)
         #def calculate_cos2(feature):
            # global negative_mean #这里需要一个正样本的mean_值
            # return cos_sim(negative_mean,feature)
         
          #如果自动生成的样本 是存在于list_temp 里的 就要重新选择
         string_result=''
         #zf=cos_sim(positive_mean,negative_mean)
         while count<array_number or early_stop>=10000 :
          if count%1000==0 and count!=0:
                print('{}negaive samples generated,please wait until all of the samples are generated'.format(count)+'\n')
                string_result=string_result+'{}negaive samples generated,please wait until all of the samples are generated'.format(count)+'\n'
                self.result_show.setText(string_result)#结果显示
          
          y=(np.random.choice(len(gene_label), 1))
          x=(np.random.choice(len(mirna_label), 1))
          index_name=x[0]
          index_gene=y[0]
          early_stop+=1
          if check_out(index_name,index_gene)==0:
             early_stop=0#如果有 就
             feature1=mirna_file[index_name:index_name+1]
             feature1=feature1.reset_index(drop=True)
             feature2=gene_file[index_gene:index_gene+1]
             feature2=feature2.reset_index(drop=True)
             temp=pd.concat([feature1,feature2],axis=1,join='outer')
             th1=calculate_cos(temp)#正样本的相似度小于0，3
             #th2=calculate_cos2(temp)#
             eu_po=calculate_euclidean(temp,positive_mean)
             #eu_ne=calculate_euclidean(temp,negative_mean)
             #var=eu_ne-eu_po #(和负的距离近，和正的距离远 距离大于064即可 负样本大于这个4.684385237024212
             #print(zf)和负样本的应该金 和正样本的应该远
             #print('该样本和负样本均值间的欧几里得距离',eu_ne)
             #print('该样本和正样本均值间的欧几里得距离',eu_po)
             #print('该样本的var value 为',var)
             #print('该样本和正样本的相似度',th1,'该样本和负样本的相似度',th2)
             #if th1<th_negative and th2>th_positive:#余弦相似度小于0,3才能称之为负样本
             if th1<=th_cosin and eu_po>th_distance:
             #if eu_po>4.6843 and var<0:
             #print(count)
              #print('该样本可行')
              if count%1000==0:
                  print(count)
              list_label.append({'0_mirna':mirna_label[index_name],'1_gene':gene_label[index_gene],'label':0})
              if count==0:
                temp0=pd.concat([feature1,feature2],axis=1,join='outer')
                count+=1
              elif count==1:                  
                temp1=pd.concat([feature1,feature2],axis=1,join='outer')
                feature_ne=pd.concat([temp0,temp1],axis=0)
                count+=1
              else:
                temp1=pd.concat([feature1,feature2],axis=1,join='outer')
                feature_ne=pd.concat([feature_ne,temp1],axis=0)
                count+=1
         #print(feature_ne)
         #print(list_label)
         fea_label=pd.DataFrame(list_label)
         feature_ne=feature_ne.reset_index(drop=True)
         positive_mean=feature.mean()
         feature_ne=pd.concat([feature_ne,fea_label],axis=1)#这个feature 就包含了所有的抽取的正样本
         #print(len(feature_ne),len(fea_label),feature_ne,fea_label)#mirna+gene的嵌入合在一起的
         #global string_history
         generate_number=len(feature_ne)
         feature_negative=feature_ne
         #feature_negative=pd.concat([feature_negative,feature_ne],axis=0)
         #加到后面
         string_history=string_history+'generation of negative samples {} completed,generated rows number{}'.format(len(feature_ne),generate_number)+'\n'
         self.state_show.setText(string_history)#结果显示
    def gene_posi(self):
         global string_history
         global gene_file
         global mirna_file
         global feature#这里是正样本
         global positive_mean
         def find_mirna(mirna):
          global mirna_label
          for i in range(0,len(mirna_label)):
            if mirna_label[i]==mirna:
                return i
         def find_gene(gene):
          global gene_label
          for i in range(0,len(gene_label)):
            if gene_label[i]==gene:
                return i
         try:   
           fname = QFileDialog.getOpenFileName(self, '打开文件','./', "Text Files (*.csv)")
           file=pd.read_csv('{}'.format(fname[0]))
           S_mirna=file['mirna']
           S_gene=file['gene']
           print('read {}，begin to generate'.format(fname[0]))
         except:
           result_string="pls select a file to generate postive file ，//mirna  gene//should be the label of columns respectively"+'\n'
           self.result_show.setText(result_string)  
         list_label=[]
         count=0
         skip=0
         string_result=''
         try:
          array_number = int(self.PS_number.toPlainText())
         except:
          array_number=len(file)
         for i in range(0,array_number):
            if i%100==0:
                print(i)
                string_result=string_result+'{}rows finished'.format(i)
                #self.result_show.setText(string_result)#结果显示
            try: #如果找到某一对才做
             
             index_name=find_mirna(S_mirna[i])
             feature1=mirna_file[index_name:index_name+1]
             feature1=feature1.reset_index(drop=True)
             index_gene=find_gene(S_gene[i])
             feature2=gene_file[index_gene:index_gene+1]
             feature2=feature2.reset_index(drop=True)
             list_label.append({'0_mirna':S_mirna[i],'1_gene':S_gene[i],'label':1})
             if count==0:
                temp0=pd.concat([feature1,feature2],axis=1,join='outer')
                count+=1
             elif count==1:                  
                temp1=pd.concat([feature1,feature2],axis=1,join='outer')
                feature=pd.concat([temp0,temp1],axis=0)
                count+=1
             else:
                temp1=pd.concat([feature1,feature2],axis=1,join='outer')
                feature=pd.concat([feature,temp1],axis=0)
                #print(index_name,index_gene,feature)#mirna+gene的嵌入合在一起的
            
            except:
                skip+=1
                #print('skip','mirna',S_mirna[i],'gene',S_gene[i])
         fea_label=pd.DataFrame(list_label)
         feature=feature.reset_index(drop=True)
         positive_mean=feature.mean()
         feature=pd.concat([feature,fea_label],axis=1)#这个feature 就包含了所有的抽取的正样本
         print(len(feature),len(fea_label),feature,fea_label)#mirna+gene的嵌入合在一起的
         #string_result=string_result+'正样本生成完毕'
         #self.result_show.setText(string_result)#结果显示
         string_history=string_history+'read {} finished,PS{} finished'.format(fname[0],array_number)+'\n'
         self.state_show.setText(string_history)#结果显示
         feature.to_csv('positive_sample.csv',index=None)
         #positive_mean.to_csv('positive_mean.csv',index=None)
         #print(positive_mean)
            
    def Load_real_data(self):
       global string_history
       global list_temp
       result_string=''
       try:
        fname = QFileDialog.getOpenFileName(self, '打开文件','./', "Text Files (*.csv)")
        #print(fname)
        df=pd.read_csv('{}'.format(fname[0]),encoding='ANSI')
        list_temp=[]
        hsa_mirna=df['miRNA']
        hsa_gene=df['Target Gene']
        for i in range(0,len(hsa_mirna)):
            list_temp.append({'miRNA':hsa_mirna[i],'gene':hsa_gene[i]})
        string_history=string_history+'hsa_MTI.csv authentic dataset loaded successfully'+'\n'
        self.state_show.setText(string_history)
       except:
         result_string="please choose the authentic dataset，the //miRNA,Target Gene//should be the columns' name respectively"+'\n'
         self.result_show.setText(result_string)
   
         
    def Load_mirna(self):
        global string_history
        global mirna_file
        global mirna_label
        global mirna_name
        #mirna_name
        string_result=''
        try:
          fname = QFileDialog.getOpenFileName(self, 'open file','./', "Text Files (*.csv)")
        #print(fname)
          mirna_name=fname[0]
          mirna_file=pd.read_csv('{}'.format(fname[0]))
          mirna_label=mirna_file['label']
          mirna_file=mirna_file.drop(['label'],axis=1)
          mirna_file=mirna_file.add_suffix('_m')
          string_history=string_history+"open mirna file{},file load succesfully".format(fname[0])+'\n'
          self.state_show.setText(string_history)
        except:
         string_result="please chose a miRNA representation file，the columns of the Mirna should be label"+'\n'
         self.result_show.setText(string_result)
    def Load_gene(self):
        global string_history
        global gene_file
        global gene_label
        string_result=''
        try:
         fname = QFileDialog.getOpenFileName(self, 'open file','./', "Text Files (*.csv)")
        #print(fname)
         gene_file=pd.read_csv('{}'.format(fname[0]))
         gene_label=gene_file['label']
         gene_file=gene_file.drop(['label'],axis=1)
         gene_file=gene_file.add_suffix('_g')
         string_history=string_history+"open gene file{},file load succesfully".format(fname[0])+'\n'
         self.state_show.setText(string_history)
        except:
         string_result="please chose a gene representation file，the columns of the gene should be label"+'\n'
         self.result_show.setText(string_result)
    def Generate_S(self):
        def find_mirna(mirna):
          global mirna_label
          for i in range(0,len(mirna_label)):
            if mirna_label[i]==mirna:
                return i
        def find_gene(gene):
          global gene_label
          for i in range(0,len(gene_label)):
            if gene_label[i]==gene:
                return i
        global string_history
        global gene_file
        global mirna_file
        global new_model
        try:
         mirna_input=str(self.Mirna_input.toPlainText())
         index_name=find_mirna(mirna_input)
         feature1=mirna_file[index_name:index_name+1]
         feature1=feature1.reset_index(drop=True)
         string_history=string_history+mirna_input+'/mirna获取正常'
        except:
          index_name=1
          feature1=mirna_file[index_name:index_name+1]
          string_history=string_history+'mirna获取异常,自动设置为1'
        try:
          gene_input=str(self.Gene_input.toPlainText())
          index_gene=find_gene(gene_input)
          feature2=gene_file[index_gene:index_gene+1]
          feature2=feature2.reset_index(drop=True)
          string_history=string_history+'/gene获取正常'+'\n'
        except:
          index_gene=1
          feature2=gene_file[index_gene:index_gene+1]
          string_history=string_history+gene_input+'/gene获取异常,自动设置为1'+'\n'
        feature=pd.concat([feature1,feature2],axis=1,join='outer')
        y_pred=new_model.predict(feature)
        y_pred2=new_model.decision_function(feature)
        string_result='mirna:{}/gene:{}  label:{}/score:{}'.format(index_name,index_gene,y_pred,y_pred2)
        #print(mirna_input,index_name)
        #print(type(feature1))
        #print(feature1)
        #print(gene_input,index_gene)
        #print(feature2)
        #print(feature)
        #print(type(feature))
        self.result_show.setText(string_result)
        self.state_show.setText(string_history)
    def Generate_P(self):
        global gene_file
        global string_history
        global mirna_file
        global mirna_label
        global gene_label
        global new_model
        
        try:
           th= int(self.continue_n.toPlainText())
        except:
           th=0 
        for i in range(th,len(mirna_file)):#len(mirna_file)
            string_state=''
            feature1=mirna_file[i:i+1]
            feature1=feature1.reset_index(drop=True)
            mirna_name=mirna_label[i]
            list_result=[]
            string_history=string_history+'mirna'+str(mirna_name)+'分数计算中'+'\n'
            print(string_history)
            #self.state_show.setText(string_history)
            for j in range(0,len(gene_label)):#len(gene_label)
                feature2=gene_file[j:j+1]
                feature2=feature2.reset_index(drop=True)
                gene_name=gene_label[j]
                feature=pd.concat([feature1,feature2],axis=1,join='outer')
                y_pred=new_model.predict(feature)
                y_pred2=new_model.decision_function(feature)
                list_result.append({'3_label':y_pred[0],'1_mirna':mirna_name,'2_gene':gene_name,'4_score':y_pred2[0]})
                if j%1000==0:
                    string_state='{}/{}'.format(j,len(gene_label))+'\n'
                    print(string_state)
                    #self.result_show.setText(string_state)
            dff=pd.DataFrame(list_result)
            dff.to_csv('{}_{}scores.csv'.format(i,mirna_name),index=None)
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    window = MyApp()
    window.show()
    sys.exit(app.exec_())
