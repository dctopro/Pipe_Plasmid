"""##!/usr/bin/env python"""
#-------------------------Librerias-------------#
try:
    from datetime import datetime
    import time
except ImportError:
    print("libreria datetime no instalado")
try:
  from sklearn.metrics import precision_score
except ImportError:
    print("libreria sklearn no instalado")
try:
    import os
except ImportError:
    print("libreria os no instalado")
try:
    import getopt
except ImportError:
    print("libreria getopt no instalado")
try:
    import csv
except ImportError:
    print("libreria csv no instalado")
try:
    import sys
except ImportError:
    print("libreria sys no instalado")
try:
    import copy
except ImportError:
    print("libreria copy no instalado")
try:
    import difflib
except ImportError:
    print("libreria difflib no instalado")#no
try:
    import glob
except ImportError:
    print("libreria glob no instalado")
try:
    import math
except ImportError:
    print("libreria math no instalado")#no
try:
    import pprint
except ImportError:
    print("libreria pprint no instalado")#no
try:
    import shutil
except ImportError:
    print("libreria shutil no instalado")#no
try:
    import re
except ImportError:
    print("libreria re no instalado")#no
try:
    from Bio import SeqIO
except ImportError:
    print("libreria SeqIO no instalado")
try:
    import numpy as np
except ImportError:
    print("libreria np de numpy no instalado")
try:
    import statistics
except ImportError:
    print("libreria statistics no instalado")
try:
    from typing import Dict, Any
except ImportError:
  print("libreria typing no instalado")
try:
    from email.mime.multipart import MIMEMultipart
except ImportError:
    print("libreria email.mime.multipart no instalado")
try:
    from email.mime.text import MIMEText
except ImportError:
    print("libreria email.mime.text no instalado")
try:
    import smtplib
except ImportError:
    print("libreria smtplib no instalado")
try:
    import gzip
except ImportError:
    print("libreria gzip no instalado")
#------------------------------------Reporte de Profundidad incluyendo script de calculo de tamano de reads
def reportQC(ruta,lista):#explorar y hacer el histograma
  #-- contar y crear histogramas---
  if str(type(lista))[-3:-2]=='t':
    for li in lista:
      arch= ruta+'/'+li+'/Reads_Clean/'
      if os.path.isfile(arch+li+'.counts.txt')==False:
        reads_f=glob.glob(arch+'*fastq.gz')
        countseqs = 0
        readscounts={}
        print('\t\nLeyendo archivos de reads para calculo de tamano ponderado ')
        for read in reads_f:
          with gzip.open(read, "rt") as handle:
            for record in SeqIO.parse(handle,"fastq"):
              countseqs += 1
              if readscounts.get(len(record.seq), None) == None:
                readscounts[len(record.seq)] = 1
              else:    
                readscounts[len(record.seq)] += 1
        j1=0
        fileout =  arch + li+'.counts.txt'
        largos=[]
        for vals in readscounts:
          largos.append(int(vals))
        largos_ord=np.sort(largos)
        with open(fileout, 'w') as fasta_file:
           for valor in largos_ord:
               j = str(valor)+ ',' + str(readscounts[valor])+','+ str('{:.3}'.format(float(readscounts[valor])*100/float(countseqs)))
               j1 += valor * readscounts[valor]
               fasta_file.write(j+'\n')
           fasta_file.write('N_seqs='+str(countseqs)+'\n')
           fasta_file.write('Prom='+str(j1/countseqs)+'\n')
      else:
        print('\t\nArchivo de conteo de reads encontrado')
  else:
    arch= ruta+'/'+lista+'/Reads_Clean/'
    if os.path.isfile(arch+lista+'.counts.txt')==False:
      arch= ruta+'/'+lista+'/Reads_Clean/'
      reads_f=glob.glob(arch+'*fastq.gz')
      countseqs = 0
      readscounts={}
      print('\t\nLeyendo archivos de reads para calculo de tamano ponderado ')
      for read in reads_f:
        with gzip.open(read, "rt") as handle:
          for record in SeqIO.parse(handle,"fastq"):
            countseqs += 1
            if readscounts.get(len(record.seq), None) == None:
              readscounts[len(record.seq)] = 1
            else:    
              readscounts[len(record.seq)] += 1
      j1=0
      fileout =  arch + lista+'.counts.txt'
      largos=[]
      for vals in readscounts:
        largos.append(int(vals))
      largos_ord=np.sort(largos)
      with open(fileout, 'w') as fasta_file:
         for valor in largos_ord:
             j = str(valor)+ ',' + str(readscounts[valor])+','+ str('{:.3}'.format(float(readscounts[valor])*100/float(countseqs)))
             j1 += valor * readscounts[valor]
             fasta_file.write(j+'\n')
         fasta_file.write('N_seqs='+str(countseqs)+'\n')
         fasta_file.write('Prom='+str(j1/countseqs)+'\n')
    else:
        print('\t\nArchivo de conteo de reads encontrado')
  #------------ crear vector de tamanos por muestra----#
  reads_len_d={}
  if str(type(lista))[-3:-2]=='t':
    for lis in lista:
      archi=ruta+'/'+lis+'/Reads_Clean/'+lis+'.counts.txt'
      if os.path.isfile(archi):
        with open(archi, 'r') as histo:
          for line in histo:
            if line[:5]=='Prom=':
              tam=float((line[5:]).strip())
              reads_len_d[lis]=int(tam)
      else:
        print('no se encontro histograma')
        reads_len_d[lis]=250        
  else:
    archi=ruta+'/'+lista+'/Reads_Clean/'+lista+'.counts.txt'
    if os.path.isfile(archi):
      with open(archi, 'r') as histo:
        for line in histo:
          if line[:5]=='Prom=':
            tam=float((line[5:]).strip())
            reads_len_d[lista]=int(tam)
    else:
      print('no se encontro histograma')
      reads_len_d[lista]=250
    print('a')
  print('\t\tTamano calculado de la distribucion de reads '+str(reads_len_d))
  #-------------------
  muestras='Muestra\tProfundidad_Ponderada\t'
  tablap=[]
  if str(type(lista))[-3:-2]=='t':
    #poner indice a la tabla
    for m in lista:
      #print('muestra'+str(m))
      spades_file=ruta+'/'+m+'/Assembly/Spades/contigs.filter1000.fasta'
      spade = {} 
      suma_tam=0
      prof_pon=0
      with open(spades_file, 'r') as input_fasta:
        for contig in SeqIO.parse(input_fasta, 'fasta'):
          line = contig.id.strip()
          fields = line.split('_')
          #campos 1= ,5= ,3= 
          spade[str(fields[1])] = float(((float(fields[5])*reads_len_d[m])/(reads_len_d[m]-127+1))*float(fields[3]))
          suma_tam+=int(fields[3])
      #print('dict'+str(spade))
      #print('suma genoma'+str(suma_tam))
      for contig in spade:
        prof_pon+=spade[contig]/suma_tam
      tablap.append(prof_pon)
  else:
    spades_file=ruta+'/'+lista+'/Assembly/Spades/contigs.filter1000.fasta'
    spade = {} 
    suma_tam=0
    prof_pon=0
    with open(spades_file, 'r') as input_fasta:
      for contig in SeqIO.parse(input_fasta, 'fasta'):
        line = contig.id.strip()
        fields = line.split('_')
        spade[str(fields[1])] = float(((float(fields[5])*reads_len_d[lista])/(reads_len_d[lista]-127+1))*float(fields[3]))
        suma_tam+=int(fields[3])
    print('dict'+str(spade))
    print('suma genoma'+str(suma_tam))
    for contig in spade:
      prof_pon+=spade[contig]/suma_tam
    tablap.append(prof_pon)
  if os.path.isdir(ruta+'/QC/'):
    print('Carpeta: '+ruta+'/QC/'+'  Encontrado!!!!')
  else:
    crearruta(ruta+'/QC/')
  tab=open(ruta+'/QC/deep_report.txt','w')
  tab.write(muestras+'\n')
  if str(type(lista))[-3:-2]=='t':
    conteo=0
    for fila in tablap:
      linea=str(lista[conteo]+'\t'+str(fila))
      conteo+=1
      tab.write(linea+'\n')
    tab.close()
  else:
    linea=str(lista+'\t'+str(tablap[0]))
    tab.write(linea+'\n')
    tab.close()
#----------------------tabla general de calificacion------------------------#
def tabladiv(ruta,lista):
  if str(type(lista))[-3:-2]=='t':
    encab='#\tContig\tPlasmidspades\tPlasflow\tProfundidad\tBD_plasmidos\tBD_No_cromosoma\tReal\n'
    tabla=ruta+'/test/mmodelo/tabla_modelo_dbnoS_crom9595.txt'
    No_ctgs=0
    with open(tabla, mode='w') as tabla_m:
      tabla_m.write(encab)
      for m in lista:
        contigs={}
        with open(ruta+'/'+m+'/Assembly/Spades/contigs.filter1000.fasta', mode='r') as cont:
          for contig in SeqIO.parse(cont, 'fasta'):
            contigs[contig.id]=[0,0,0,0,0,0]#creacion del dictado por muestra 
        with open(ruta+'/test/True/'+m+'/verdaderos.txt', mode='r') as verd:
          for v in verd:
            if contigs.get(v.strip(),None)!=None:
              contigs[v.strip()][5]=1
        #plasmidspades
        if os.path.isfile(ruta+'/'+m+'/Assembly/Plasmitest/Plasmidspades/nodos_plasmidspades.txt'):
          with open(ruta+'/'+m+'/Assembly/Plasmitest/Plasmidspades/nodos_plasmidspades.txt', mode='r') as plasmidspades:
            for pred in plasmidspades:
              if contigs.get(pred.strip(),None)!=None:
                contigs[pred.strip()][0]=1
        #plasflow
        if os.path.isfile(ruta+'/'+m+'/Assembly/Plasmitest/Plasflow/Plasflow.txt'):
          with open(ruta+'/'+m+'/Assembly/Plasmitest/Plasflow/Plasflow.txt', mode='r') as plasflow:
            for pred in plasflow:
              if contigs.get(pred.strip(),None)!=None:
                contigs[pred.strip()][1]=1
        #profundidad
        if os.path.isfile(ruta+'/'+m+'/Assembly/Plasmitest/Profundidad/nodos_cov1.5.txt'):
          with open(ruta+'/'+m+'/Assembly/Plasmitest/Profundidad/nodos_cov1.5.txt', mode='r') as profundidad:
            for pred in profundidad:
              if contigs.get(pred.strip(),None)!=None:
                contigs[pred.strip()][2]=1        
        #homologia plasmidos
        if os.path.isfile(ruta+'/'+m+'/Assembly/Plasmitest/Homologia/PlasmidsDB_nodos_noS.txt'):
          with open(ruta+'/'+m+'/Assembly/Plasmitest/Homologia/PlasmidsDB_nodos_noS.txt', mode='r') as homologia:
            for pred in homologia:
              if contigs.get(pred.strip(),None)!=None:
                contigs[pred.strip()][3]=1
        #homologia con cromosoma
        if os.path.isfile(ruta+'/'+m+'/Assembly/Plasmitest/Homologia/No_chromosomeDB_nodos_noS_9595.txt'):
          with open(ruta+'/'+m+'/Assembly/Plasmitest/Homologia/No_chromosomeDB_nodos_noS_9595.txt', mode='r') as homologiac:
            for pred in homologiac:
              if contigs.get(pred.strip(),None)!=None:
                contigs[pred.strip()][4]=1
        for c in contigs:
          linea_c=str(No_ctgs)+'\t'+str(c)+'\t'+str(contigs[c][0])+'\t'+str(contigs[c][1])+'\t'+str(contigs[c][2])+'\t'+str(contigs[c][3])+'\t'+str(contigs[c][4])+'\t'+str(contigs[c][5])+'\n'
          No_ctgs+=1
          tabla_m.write(linea_c)
      print('lo quiero y quierase  ') 
  else:
    tabla=ruta+'/'+lista+'/tabla_modelo.txt'
    with open(tabla, mode='w') as tabla_m:
      contigs={}
      with open(ruta+'/'+lista+'/Assembly/Spades/contigs.filter1000.fasta', mode='r') as cont:
        for contig in SeqIO.parse(cont, 'fasta'):
          contigs[contig.id]=[0,0,0,0,0,0]#creacion del dictado por muestra 
      with open(ruta+'/test/True/'+lista+'/verdaderos.txt', mode='r') as verd:
        for v in verd:
          if contigs.get(v.strip(),None)!=None:
            contigs[v.strip()][5]=1
      #plasmidspades
      if os.path.isfile(ruta+'/'+lista+'/Assembly/Plasmitest/Plasmidspades/nodos_plasmidspades.txt'):
        with open(ruta+'/'+lista+'/Assembly/Plasmitest/Plasmidspades/nodos_plasmidspades.txt', mode='r') as plasmidspades:
          for pred in plasmidspades:
            if contigs.get(pred.strip(),None)!=None:
              contigs[pred.strip()][0]=1
      #plasflow
      if os.path.isfile(ruta+'/'+lista+'/Assembly/Plasmitest/Plasflow/Plasflow.txt'):
        with open(ruta+'/'+lista+'/Assembly/Plasmitest/Plasflow/Plasflow.txt', mode='r') as plasflow:
          for pred in plasflow:
            if contigs.get(pred.strip(),None)!=None:
              contigs[pred.strip()][1]=1
      #profundidad
      if os.path.isfile(ruta+'/'+lista+'/Assembly/Plasmitest/Profundidad/nodos_cov1.5.txt'):
        with open(ruta+'/'+lista+'/Assembly/Plasmitest/Profundidad/nodos_cov1.5.txt', mode='r') as profundidad:
          for pred in profundidad:
            if contigs.get(pred.strip(),None)!=None:
              contigs[pred.strip()][2]=1        
      #homologia plasmidos
      if os.path.isfile(ruta+'/'+lista+'/Assembly/Plasmitest/Homologia/PlasmidsDB_nodos_noS.txt'):
        with open(ruta+'/'+lista+'/Assembly/Plasmitest/Homologia/PlasmidsDB_nodos.txt', mode='r') as homologia:
          for pred in homologia:
            if contigs.get(pred.strip(),None)!=None:
              contigs[pred.strip()][3]=1
      for c in contigs:
        linea_c=str(c)+'\t'+str(contigs[c][0])+'\t'+str(contigs[c][1])+'\t'+str(contigs[c][2])+'\t'+str(contigs[c][3])+'\t'+str(contigs[c][4])+'\t'+str(contigs[c][5])+'\n'
        tabla_m.write(linea_c) 
    print('bien camilin, dele y vera que acaba. :*')
#-------------- tabla de recall para todas las muestras---------------------#Experimento homologia cromosomas
def recallhomocro(ruta,lista):
  print(lista)
  muestras='ID\tCOV\t'
  if str(type(lista))[-3:-2]=='t':
    tablap=np.zeros((2500,len(lista)+2))
    tablar=np.zeros((2500,len(lista)+2))
    tablaf=np.zeros((2500,len(lista)+2))
    vp=np.zeros(len(lista))
    me=np.zeros(len(lista))
    sc=np.zeros(len(lista))
    pc=np.zeros(len(lista))
    tp=np.zeros(len(lista))
    pre=np.zeros(len(lista))
    rec=np.zeros(len(lista))
    fsc=np.zeros(len(lista))
    for m in lista[:-1]:
      muestras+=str(m)+'\t'
    muestras+=str(lista[-1:])[2:-2]
  else:
    muestras+=str(lista)
    tablap=np.zeros((2500,3))
    tablar=np.zeros((2500,3))
    tablaf=np.zeros((2500,3))
    vp=0
    me=0
    sc=0
    pc=0
    tp=0
    pre=0
    rec=0
    fsc=0
  itam=0
  jtam=0
  fila=0
  for i in range(50,100):
    jtam=0 
    for j in range(50,100):
      tablap[fila][0]=i#indicar combinacion de id y cov 
      tablar[fila][0]=i
      tablaf[fila][0]=i
      tablap[fila][1]=j
      tablar[fila][1]=j
      tablaf[fila][1]=j
      neg_pr=np.zeros(len(lista))
      verdad=0
      conteo_m=0
      conteo_n=0
      #sys.exit()
      if str(type(lista))[-3:-2]=='t':
        vp=np.zeros(len(lista))
        tp=np.zeros(len(lista))
        #poner indice a la tabla
        for m in lista:
          #print('muestra'+str(m))
          blast_csv=ruta+'/'+m+'/Assembly/Plasmitest/Homologia/No_chromosomeDB_blastn3_out.csv'
          #blast_csv=ruta+'/'+m+'/Assembly/Plasmitest/Homologia/No_chromosomeDB_blastn3_olddb_out.csv'#prueba old DB#
          #blast_csv=ruta+'/'+m+'/Assembly/Plasmitest/Homologia/No_chromosomeDB_blastn_out.csv'#prueba old DB
          nodos={}
          contact=''
          contpas=''
          cvacum=0
          with open(blast_csv, mode='r') as blast_file:
            for line in blast_file:
              line = line.strip()
              fields = line.split(',')
              contact=fields[0]
              if contact!=contpas:
                cvs= (float(fields[3])/float(fields[9]))*100
                cv = float(fields[10])
                ide = float(fields[2])
                if ide >= i or cv >= j :#fue mejor usar solo cobertura de query
                  if nodos.get(fields[0],None)==None:
                    conteo_n+=1
                    nnod=str(fields[0]).split('_')
                    nodos[nnod[1]]=1#son las etiquetas de los contigs de spades extraidas del alineamiento BD y spades
                  else:
                    print('repetido')
              contpas=contact  
          #------------------------seccion para mirar los nodos excluidos de los cromosomas
          conteo_n=0
          cantidad_ex=0
          nodos_ex={}
          nodos_spades={}
          spades_fi=ruta+'/'+m+'/Assembly/Spades/contigs.filter1000.fasta'
          with open(spades_fi, mode='r') as multifasta:
            for contig in SeqIO.parse(multifasta, 'fasta'):
              if nodos_spades.get(((str(contig.id).split('_'))[1]),None)==None:
                nodos_spades[(str(contig.id).split('_'))[1]]=1
          cantidad_ex=0
          for n in nodos_spades:
            if nodos.get(n,None)==None:
              nodos_ex[n]=1
              conteo_n+=1
            else:
              cantidad_ex+=1
          #print('cantidad de nodos de cromosoma excluidos '+str(cantidad_ex))
          #------------------------------------------------------fin seccion excluidos         
          me[conteo_m]=conteo_n#Cantidad de nodos predichos por homologia que pasan el filtro
          conteo_n=0
          #print('tamanio de nodos que pasaron el filtrado '+str(len(nodos)))
          true=ruta+'/test/True/'+m+'/verdaderos.txt'
          with open(true, mode='r') as verdaderos:
            for line in verdaderos:
              field = (line.strip()).split('_')
              tp[conteo_m]+=1
              #if nodos.get(str(field[1]),None)!=None:
              if nodos_ex.get(str(field[1]),None)!=None:
                vp[conteo_m]+=1        
          pre[conteo_m]=(vp[conteo_m]/me[conteo_m])
          rec[conteo_m]=(vp[conteo_m]/(tp[conteo_m]))
          fsc[conteo_m]=2*(((pre[conteo_m])*(rec[conteo_m]))/((pre[conteo_m])+(rec[conteo_m])))
          fsc[np.isnan(fsc)]=0#eliminar los nan
          tablap[fila][conteo_m+2]=round(pre[conteo_m],4)#precision
          tablar[fila][conteo_m+2]=round(rec[conteo_m],4)#recall
          tablaf[fila][conteo_m+2]=round(fsc[conteo_m],4)#fscore
          conteo_m+=1
      else:
        vp=0
        tp=0
        blast_csv=ruta+'/'+lista+'/Assembly/Plasmitest/Homologia/No_chromosomeDB_blastn3.csv'
        nodos={}
        contact=''
        contpas=''
        print(lista)
        with open(blast_csv, mode='r') as blast_file:
          for line in blast_file:
            line = line.strip()
            fields = line.split(',')
            contact=fields[0]
            if contact!=contpas:
              cvs= (float(fields[3])/float(fields[9]))*100
              cv = float(fields[10])
              ide = float(fields[2])
              #print(ide,cv,cvs)
              if ide >= i or cv >= j:#prueba tomando buenos hits y luego excluyendo
                if nodos.get(fields[0],None)==None:
                  conteo_n+=1
                  nnod=str(fields[0]).split('_')
                  nodos[nnod[1]]=1#son los numeros de los contigs
                else:
                  print('repetido')
            contpas=contact
        #------------------------seccion para mirar los nodos excluidos de los cromosomas
        conteo_n=0
        cantidad_ex=0
        nodos_ex={}
        nodos_spades={}
        spades_fi=ruta+'/'+lista+'/Assembly/Spades/contigs.filter1000.fasta'
        with open(spades_fi, mode='r') as multifasta:
          for contig in SeqIO.parse(multifasta, 'fasta'):
            if nodos_spades.get(((str(contig.id).split('_'))[1]),None)==None:
              nodos_spades[(str(contig.id).split('_'))[1]]=1
        cantidad_ex=0
        for n in nodos_spades:
          if nodos.get(n,None)==None:
            nodos_ex[n]=1
            conteo_n+=1
          else:
            cantidad_ex+=1
        #print('cantidad de nodos de cromosoma excluidos '+str(cantidad_ex))
        #------------------------------------------------------fin seccion excluidos 
        me=conteo_n#Cantidad de nodos homologos que pasan el filtro del alineamiento
        #print(me)
        conteo_n=0
        #print(nodos)
        true=ruta+'/test/True/'+lista+'/verdaderos.txt'#Los verdaderos contigs son seleccionados 
        #de un aliniamiento blast con metricas 90 70
        with open(true, mode='r') as verdaderos:
          for line in verdaderos:
            field = (line.strip()).split('_')
            tp+=1
            #if nodos.get(str(field[1]),None)!=None:#prueba excluyendo nodos
            if nodos_ex.get(str(field[1]),None)!=None: 
              vp+=1
        try:
          rec=(vp/me)
        except:
          rec=0
          print('indeterminado en rec')
        try:
          pre=(vp/(tp))
        except:
          pre=0
          print('indeterminado en pre')
        try:
          fsc=2*(((vp/me)*(vp/(tp)))/((vp/me)+(vp/(tp))))
        except:
          print('operacion indeterminada x/0')
        tablap[fila][2]=round(pre,4)#precision
        tablar[fila][2]=round(rec,4)#recall
        tablaf[fila][2]=round(fsc,4)#fscore
      jtam+=1
      fila+=1
    itam+=1
  print('presicion '+str(tablap))
  print('recall '+str(tablar))
  print('f-score '+str(tablaf))
  tab=open(ruta+'/test/homologia/tabla_croNO_fscore_nos_or_3hit.txt','w')
  tab.write(muestras+'\n')
  for fila in tablaf:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
#-------------- tabla de recall para todas las muestras---------------------#Experimento homologia plasmidos
def recallhomo(ruta,lista):
  print(lista)
  muestras='ID\tCOV\t'
  if str(type(lista))[-3:-2]=='t':
    tablap=np.zeros((1600,len(lista)+2))
    tablar=np.zeros((1600,len(lista)+2))
    tablaf=np.zeros((1600,len(lista)+2))
    vp=np.zeros(len(lista))
    me=np.zeros(len(lista))
    sc=np.zeros(len(lista))
    pc=np.zeros(len(lista))
    tp=np.zeros(len(lista))
    pre=np.zeros(len(lista))
    rec=np.zeros(len(lista))
    fsc=np.zeros(len(lista))
    for m in lista[:-1]:
      muestras+=str(m)+'\t'
    muestras+=str(lista[-1:])
  else:
    muestras+=str(lista)
    tablap=np.zeros((1600,3))
    tablar=np.zeros((1600,3))
    tablaf=np.zeros((1600,3))
    vp=0
    me=0
    sc=0
    pc=0
    tp=0
    pre=0
    rec=0
    fsc=0
  itam=0
  jtam=0
  fila=0
  for i in range(80,100):
    jtam=0 
    for j in range(20,100):
      tablap[fila][0]=i#indicar combinacion de id y cov 
      tablar[fila][0]=i
      tablaf[fila][0]=i
      tablap[fila][1]=j
      tablar[fila][1]=j
      tablaf[fila][1]=j
      neg_pr=np.zeros(len(lista))
      verdad=0
      conteo_m=0
      conteo_n=0
      #sys.exit()
      if str(type(lista))[-3:-2]=='t':
        vp=np.zeros(len(lista))
        tp=np.zeros(len(lista))
        #poner indice a la tabla
        for m in lista:
          #print('muestra'+str(m))
          #blast_csv=ruta+'/'+m+'/Assembly/Plasmitest/Homologia/PlasmidsDB_blastn3_out.csv'#prueba old DB
          blast_csv=ruta+'/'+m+'/Assembly/Plasmitest/Homologia/PlasmidsDB_blastn3_olddb_out.csv'
          nodos={}
          contact=''
          contpas=''
          cvacum=0
          with open(blast_csv, mode='r') as blast_file:
            for line in blast_file:
              line = line.strip()
              fields = line.split(',')
              contact=fields[0]
              if contact!=contpas:
                cvs= (float(fields[3])/float(fields[9]))*100
                cv = float(fields[10])
                ide = float(fields[2])
                #if ide >= i and (cv >= j or cvs >= j):
                if ide >= i and cv >=j:
                  if nodos.get(fields[0],None)==None:
                    conteo_n+=1
                    nnod=str(fields[0]).split('_')
                    nodos[nnod[1]]=1#son las etiquetas de los contigs de spades extraidas del alineamiento BD y spades
                  else:
                    print('repetido')
              contpas=contact           
          me[conteo_m]=conteo_n#Cantidad de nodos predichos por homologia que pasan el filtro
          conteo_n=0
          #print('tamanio de nodos que pasaron el filtrado '+str(len(nodos)))
          true=ruta+'/test/True/'+m+'/verdaderos.txt'
          with open(true, mode='r') as verdaderos:
            for line in verdaderos:
              field = (line.strip()).split('_')
              tp[conteo_m]+=1
              if nodos.get(str(field[1]),None)!=None:
                vp[conteo_m]+=1        
          pre[conteo_m]=(vp[conteo_m]/me[conteo_m])
          rec[conteo_m]=(vp[conteo_m]/(tp[conteo_m]))
          fsc[conteo_m]=2*(((pre[conteo_m])*(rec[conteo_m]))/((pre[conteo_m])+(rec[conteo_m])))
          fsc[np.isnan(fsc)]=0#eliminar los nan
          tablap[fila][conteo_m+2]=round(pre[conteo_m],4)#precision
          tablar[fila][conteo_m+2]=round(rec[conteo_m],4)#recall
          tablaf[fila][conteo_m+2]=round(fsc[conteo_m],4)#fscore
          conteo_m+=1
      else:
        vp=0
        tp=0
        blast_csv=ruta+'/'+lista+'/Assembly/Plasmitest/Homologia/PlasmidsDB_blastn3_out.csv'
        nodos={}
        contact=''
        contpas=''
        with open(blast_csv, mode='r') as blast_file:
          for line in blast_file:
            line = line.strip()
            fields = line.split(',')
            contact=fields[0]
            if contact!=contpas:
              cvs= (float(fields[3])/float(fields[9]))*100
              cv = float(fields[10])
              ide = float(fields[2])
              if ide >= i and (cv >= j or cvs >= j):
                if nodos.get(fields[0],None)==None:
                  conteo_n+=1
                  nnod=str(fields[0]).split('_')
                  nodos[nnod[1]]=1#son los numeros de los contigs
                else:
                  print('repetido')
            contpas=contact
        me=conteo_n#Cantidad de nodos homologos que pasan el filtro del alineamiento
        #print(me)
        conteo_n=0
        #print(nodos)
        true=ruta+'/test/True/'+lista+'/verdaderos.txt'#Los verdaderos contigs son seleccionados 
        #de un aliniamiento blast con metricas 90 70
        with open(true, mode='r') as verdaderos:
          for line in verdaderos:
            field = (line.strip()).split('_')
            tp+=1
            if nodos.get(str(field[1]),None)!=None:
              vp+=1
        rec=(vp/me)
        pre=(vp/(tp))
        try:
          fsc=2*(((vp/me)*(vp/(tp)))/((vp/me)+(vp/(tp))))
        except:
          print('operacion indeterminada x/0')
        tablap[fila][2]=round(pre,4)#precision
        tablar[fila][2]=round(rec,4)#recall
        tablaf[fila][2]=round(fsc,4)#fscore
      jtam+=1
      fila+=1
    itam+=1
  print('presicion '+str(tablap))
  print('recall '+str(tablar))
  print('f-score '+str(tablaf))
  tab=open(ruta+'/test/homologia/tabla_fscore_olddb.txt','w')
  tab.write(muestras+'\n')
  for fila in tablaf:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
##
#-----------------------------------------Experimento
def explorprof(ruta,lista):#explorar y hacer el histograma
  print(lista)
  #muestras='Porcentaje\t'
  read_len={'SRR2244244':'301','SRR4025861':'251','SRR5167853':'300','SRR5168393':'300','SRR6348595':'300','SRR6514350':'300','SRR8175017':'151','SRR8607467':'300','SRR9042857':'251','ERR2929690':'151','SRR3465532':'300','SRR4245476':'251','SRR5168236':'300','SRR5168488':'300','SRR6519357':'150','SRR8607449':'300','SRR8668707':'251','ERR2929692':'151','SRR3465557':'300','SRR5146463':'300','SRR5168385':'300','SRR5714002':'300','SRR6514141':'300','SRR6675860':'300','SRR8607459':'300','SRR8778550':'300'}
  muestras=''
  print(type(lista))
  if str(type(lista))[-3:-2]=='t':
    tablap=np.zeros((100,len(lista)+1))
    tablar=np.zeros((100,len(lista)+1))
    tablaf=np.zeros((100,len(lista)+1))
    tablahis=np.zeros((200,len(lista)*3))
    vp=np.zeros(len(lista))
    me=np.zeros(len(lista))
    sc=np.zeros(len(lista))
    pc=np.zeros(len(lista))
    tp=np.zeros(len(lista))
    pre=np.zeros(len(lista))
    rec=np.zeros(len(lista))
    fsc=np.zeros(len(lista))
    for m in lista[:-1]:
      muestras+=str(m)+'_ID\tNucle_cov\tTipo_sec\t'
    muestras+=str(lista[-1:])[2:-2]+'_ID\tKmer_Cov\tTipo_sec'
  else:
    muestras+=str(lista)
    tablap=np.zeros((100,2))
    tablar=np.zeros((100,2))
    tablaf=np.zeros((100,2))
    vp=0
    me=0
    sc=0
    pc=0
    tp=0
    pre=0
    rec=0
    fsc=0
    itam=0
    jtam=0
    fila=0
  if str(type(lista))[-3:-2]=='t':
    vp=np.zeros(len(lista))
    tp=np.zeros(len(lista))
    conteo_m=0
    #poner indice a la tabla
    for m in lista:
      #print('muestra'+str(m))
      spades_file=ruta+'/'+m+'/Assembly/Spades/contigs.filter1000.fasta'
      spade = {}  # type: Dict[Any, float]
      #cov_nodes = {}
      #quan=[]
      with open(spades_file, 'r') as input_fasta:
        for contig in SeqIO.parse(input_fasta, 'fasta'):
            line = contig.id.strip()
            fields = line.split('_')
            if read_len.get(m,None)!=None:
              spade[str(fields[1])] = float(float(fields[5])*int(read_len[m])/(int(read_len[m])-127+1))
      true=ruta+'/test/True/'+m+'/verdaderos.txt'
      verdad={}
      with open(true, mode='r') as verdaderos:
        for line in verdaderos:
          num=str(((line.strip()).split('_'))[1])
          verdad[str(num)]=1
      fila=0
      print(len(verdad))
      for contig in spade:
        if len(verdad)>0:
          if verdad.get(contig,None)!=None:
            tablahis[fila][conteo_m]=int(contig)
            tablahis[fila][conteo_m+1]=int(spade[contig])
            tablahis[fila][conteo_m+2]=1
          else:
            tablahis[fila][conteo_m]=int(contig)
            tablahis[fila][conteo_m+1]=int(spade[contig])
            tablahis[fila][conteo_m+2]=0
        else:
          print('ingresa al except ')
          tablahis[fila][conteo_m]=int(contig)
          tablahis[fila][conteo_m+1]=int(spade[contig])
          tablahis[fila][conteo_m+2]=0
        fila+=1
      conteo_m+=3
  else:
    spades_file=ruta+'/'+lista+'/Assembly/Spades/contigs.filter1000.fasta'
    spade = {}  # type: Dict[Any, float]
    #cov_nodes = {}
    #quan=[]
    with open(spades_file, 'r') as input_fasta:
          for contig in SeqIO.parse(input_fasta, 'fasta'):
              line = contig.id.strip()
              fields = line.split('_')
              spade[str(fields[1])] = float(fields[5])
    true=ruta+'/test/True/'+m+'/verdaderos.txt'
    verdad={}
    with open(true, mode='r') as verdaderos:
      for line in verdaderos:
        num=str(((line.strip()).split('_'))[1])
        verdad[str(num)]=1
    histo=ruta+'/test/profundidad/histograma.txt'       
    with open(histo, 'w') as txt_sal:
      for contig in spade:
        if verdad.get(contig,None)!=None:
          txt_sal.write(str(contig)+'\t'+str(spade[contig])+'\t1\n')
        else:
          txt_sal.write(str(contig)+'\t'+str(spade[contig])+'\t0\n')
  if str(type(lista))[-3:-2]=='t':
    tab=open(ruta+'/test/profundidad/histograma_nucl.txt','a')
    tab.write(muestras+'\n')
    for fila in tablahis:
      linea=''
      for col in fila[:-1]:
        linea+=str(col)+'\t'
      linea+=str(fila[-1])
      tab.write(linea+'\n')
    tab.close()
    """quan.append(float(fields[5]))  #all contig coverages
    q1 = np.percentile(quan, 25)
    q3 = np.percentile(quan, 75)
    limite = q3 + ((q3-q1)*3) #external outlier limit
    cov_file=salida+'/Profundidad/nodos_cov.txt'
    print ("\n coverage " + str(limite) +"\n")
    with open(cov_file, 'w') as txt_sal:
      for contig in spade:
          if float(spade[contig]) > limite:
             cov_nodes[contig]=1
             txt_sal.write(str(contig)+'\n')
    return cov_nodes"""
#------------------------recall coverage------------------------------------#Experimento
def recallprof(ruta,lista):
  #print(lista)
  read_len={'SRR2244244':301,'SRR4025861':251,'SRR5167853':300,'SRR5168393':300,'SRR6348595':300,'SRR6514350':300,'SRR8175017':151,'SRR8607467':300,'SRR9042857':251,'ERR2929690':151,'SRR3465532':300,'SRR4245476':251,'SRR5168236':300,'SRR5168488':300,'SRR6519357':150,'SRR8607449':300,'SRR8668707':251,'ERR2929692':151,'SRR3465557':300,'SRR5146463':300,'SRR5168385':300,'SRR5714002':300,'SRR6514141':300,'SRR6675860':300,'SRR8607459':300,'SRR8778550':300}
  muestras='Porcentaje\t'
  #print(type(lista))
  if str(type(lista))[-3:-2]=='t':
    tablap=np.zeros((100,len(lista)+1))
    tablar=np.zeros((100,len(lista)+1))
    tablaf=np.zeros((100,len(lista)+1))
    vp=np.zeros(len(lista))
    fp=np.zeros(len(lista))
    vn=np.zeros(len(lista))
    fn=np.zeros(len(lista))
    pre=np.zeros(len(lista))
    rec=np.zeros(len(lista))
    fsc=np.zeros(len(lista))
    for m in lista[:-1]:
      muestras+=str(m)+'\t'
    muestras+=str(lista[-1:])[2:-2]
  else:
    muestras+=str(lista)
    tablap=np.zeros((100,2))
    tablar=np.zeros((100,2))
    tablaf=np.zeros((100,2))
    vp=0
    fp=0
    vn=0
    fn=0
    pre=0
    rec=0
    fsc=0
    fila=0
  for i in range(0,100):
    tablap[i][0]=(1.5+(i*1.5/100))#indicar porcentaje de exclucion
    tablar[i][0]=(1.5+(i*1.5/100))
    tablaf[i][0]=(1.5+(i*1.5/100))
    if str(type(lista))[-3:-2]=='t':
      conteo_m=0
      #poner indice a la tabla
      for m in lista:
        spades_file=ruta+'/'+m+'/Assembly/Spades/contigs.filter1000.fasta'
        spade = {}  # type: Dict[Any, float]
        cov=[]
        tam=[]
        #print('archivo de lecturas de la muestra '+str(m)+' : '+str(os.path.getsize(spades_file)))
        with open(spades_file, 'r') as input_fasta:
          for contig in SeqIO.parse(input_fasta, 'fasta'):
            line = contig.id.strip()
            fields = line.split('_')
            spade[str(fields[1])] = float(fields[5])
            cov.append(float(fields[5]))
            tam.append(float(fields[3]))
        #print('muestra '+str(lista)+' cantidad de contigs '+str(len(cov)))
        q1 = np.percentile(cov, 25)
        q3 = np.percentile(cov, 75)
        limite = q3 + ((q3-q1)*(1.5+(i*1.5/100))) #external outlier limit
        """media=np.mean(cov)
        mediapon=np.average(cov,weights=tam)
        dev_st=np.std(cov)
        limite=mediapon+((dev_st)/200)*i"""
        true=ruta+'/test/True/'+m+'/verdaderos.txt'
        verdad={}
        #print('archivo de nodos de referencia de la muestra '+str(m)+' : '+str(os.path.getsize(true)))
        with open(true, mode='r') as verdaderos:
          for line in verdaderos:
            num=str(((line.strip()).split('_'))[1])
            verdad[str(num)]=1
        
        for contig in spade:
          if float(spade[contig])>limite:
            if verdad.get(contig):
              vp[conteo_m]+=1
            else:
              fp[conteo_m]+=1
          else:
            if verdad.get(contig):
              fn[conteo_m]+=1
            else:
              vn[conteo_m]+=1
        pre=(vp[conteo_m])/(vp[conteo_m]+fp[conteo_m])
        try:
          rec=(vp[conteo_m])/(vp[conteo_m]+fn[conteo_m])
        except:
          rec=0
        tablap[i][conteo_m+1]=pre
        tablar[i][conteo_m+1]=rec
        try:
          tablaf[i][conteo_m+1]=2*((pre*rec)/(pre+rec))
        except:
          tablaf[i][conteo_m+1]=0
        conteo_m+=1
    else:
      spades_file=ruta+'/'+lista+'/Assembly/Spades/contigs.filter1000.fasta'
      spade = {}  # type: Dict[Any, float]
      cov=[]
      tam=[]
      with open(spades_file, 'r') as input_fasta:
        for contig in SeqIO.parse(input_fasta, 'fasta'):
          line = contig.id.strip()
          fields = line.split('_')
          spade[str(fields[1])] = float(fields[5])
          if read_len.get(lista,None)!=None:
            cov.append(float(float(fields[5])*read_len[lista]/(read_len[lista]-127+1)))
          tam.append(float(fields[3]))
      #print('muestra '+str(lista)+' cantidad de contigs '+str(len(cov)))
      media=np.mean(cov)
      mediapon=np.average(cov,weights=tam)
      dev_st=np.std(cov)
      true=ruta+'/test/True/'+lista+'/verdaderos.txt'
      verdad={}
      with open(true, mode='r') as verdaderos:
        for line in verdaderos:
          num=str(((line.strip()).split('_'))[1])
          verdad[str(num)]=1
      limite=mediapon+((dev_st*3)/100)*i
      for contig in spade:
        if float(spade[contig])>limite:
          if verdad.get(contig):
            vp+=1
          else:
            fp+=1
        else:
          if verdad.get(contig):
            fn+=1
          else:
            vn+=1
      pre=(vp)/(vp+fp)
      rec=(vp)/(vp+fn)
      tablap[i][1]=pre
      tablar[i][1]=rec
      try:
        tablaf[i][1]=2*((pre*rec)/(pre+rec))
      except:
        tablaf[i][1]=0
  print(tablap)
  print(tablar)
  print(tablaf)
  tab=open(ruta+'/test/profundidad/tablafscore_at_eat.txt','a')
  tab.write(muestras+'\n')
  for fila in tablaf:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
#-------------- tabla de recall para todas las muestras plasmidspades---------------------#Experimento
def recallplasmid(ruta,lista):
  print(lista)
  muestras='ID\tCOV\t'
  if str(type(lista))[-3:-2]=='t':
    tablap=np.zeros((1600,len(lista)+2))
    tablar=np.zeros((1600,len(lista)+2))
    tablaf=np.zeros((1600,len(lista)+2))
    vp=np.zeros(len(lista))
    me=np.zeros(len(lista))
    sc=np.zeros(len(lista))
    pc=np.zeros(len(lista))
    tp=np.zeros(len(lista))
    pre=np.zeros(len(lista))
    rec=np.zeros(len(lista))
    fsc=np.zeros(len(lista))
    for m in lista[:-1]:
      muestras+=str(m)+'\t'
    muestras+=str(lista[-1:])
  else:
    muestras+=str(lista)
    tablap=np.zeros((1600,3))
    tablar=np.zeros((1600,3))
    tablaf=np.zeros((1600,3))
    vp=0
    me=0
    sc=0
    pc=0
    tp=0
    pre=0
    rec=0
    fsc=0
  itam=0
  jtam=0
  fila=0
  for i in range(80,100):
    jtam=0 
    for j in range(20,100):
      tablap[fila][0]=i#indicar combinacion de id y cov 
      tablar[fila][0]=i
      tablaf[fila][0]=i
      tablap[fila][1]=j
      tablar[fila][1]=j
      tablaf[fila][1]=j
      neg_pr=np.zeros(len(lista))
      verdad=0
      conteo_m=0
      conteo_n=0
      #sys.exit()
      if str(type(lista))[-3:-2]=='t':
        vp=np.zeros(len(lista))
        tp=np.zeros(len(lista))
        #poner indice a la tabla
        for m in lista:
          #print('muestra'+str(m))
          spades=ruta+'/'+m+'/Assembly/Spades/contigs.filter1000.fasta'
          renglones = os.popen("grep -c '>' "+spades).read()
          sc[conteo_m] = int(renglones.strip())#cabtidad de nodos spades
          plspades=ruta+'/'+m+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.fasta'
          renglones = os.popen("grep -c '>' "+plspades).read()
          #print(renglones)
          pc[conteo_m] = int(renglones.strip())#cabtidad de nodos Plasmidspades
          blast_csv=ruta+'/'+m+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.blastn3_csv'
          nodos={}
          contact=''
          contpas=''
          cvacum=0
          with open(blast_csv, mode='r') as blast_file:
            for line in blast_file:
              line = line.strip()
              fields = line.split(',')
              contact=fields[0]
              if contact!=contpas:
                cvs= (float(fields[4])/float(fields[5]))*100
                cv = float(fields[3])
                ide = float(fields[2])
                if ide >= i and (cv >= j or cvs >= j):
                  if nodos.get(fields[0],None)==None:
                    conteo_n+=1
                    nnod=str(fields[0]).split('_')
                    nodos[nnod[1]]=1#son las etiquetas de los contigs de spades extraidas del alineamiento entre plasmid y spades
                  else:
                    print('repetido')
              contpas=contact           
          me[conteo_m]=conteo_n#Cantidad de nodos predichos por plasmidspades que pasan el filtro
          conteo_n=0
          #print('tamanio de nodos que pasaron el filtrado '+str(len(nodos)))
          true=ruta+'/test/True/'+m+'/verd_mau.txt'
          with open(true, mode='r') as verdaderos:
            for line in verdaderos:
              fields = line.strip()
              tp[conteo_m]+=1
              if nodos.get(str(fields),None)!=None:
                vp[conteo_m]+=1        
          pre[conteo_m]=(vp[conteo_m]/me[conteo_m])
          #print(rec[ii])
          rec[conteo_m]=(vp[conteo_m]/(tp[conteo_m]))
          #print(pre[ii])
          fsc[conteo_m]=2*(((pre[conteo_m])*(rec[conteo_m]))/((pre[conteo_m])+(rec[conteo_m])))
          #print(fsc[ii])
          fsc[np.isnan(fsc)]=0#eliminar los nan
          tablap[fila][conteo_m+2]=round(pre[conteo_m],4)#precision
          tablar[fila][conteo_m+2]=round(rec[conteo_m],4)#recall
          tablaf[fila][conteo_m+2]=round(fsc[conteo_m],4)#fscore
          conteo_m+=1
      else:
        vp=0
        tp=0
        """pre=0
        rec=0
        fsc=0"""
        #print('muestra'+str(lista))
        spades=ruta+'/'+lista+'/Assembly/Spades/contigs.filter1000.fasta'
        renglones = os.popen("grep -c '>' "+spades).read()
        sc = renglones.strip()#cabtidad de nodos spades
        plspades=ruta+'/'+lista+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.fasta'
        renglones = os.popen("grep -c '>' "+plspades).read()
        pc = int(renglones.strip())#cabtidad de nodos Plasmidspades
        #print('nodos de plasmidspades'+str(pc))
        blast_csv=ruta+'/'+lista+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.blastn3_csv'
        nodos={}
        contact=''
        contpas=''
        with open(blast_csv, mode='r') as blast_file:
          for line in blast_file:
            line = line.strip()
            fields = line.split(',')
            contact=fields[0]
            if contact!=contpas:
              cvs= (float(fields[4])/float(fields[5]))*100
              cv = float(fields[3])
              ide = float(fields[2])
              if ide >= i and (cv >= j or cvs >= j):
                if nodos.get(fields[0],None)==None:
                  conteo_n+=1
                  nnod=str(fields[0]).split('_')
                  nodos[nnod[1]]=1#son los numeros de los contigs
                else:
                  print('repetido')
            contpas=contact
        me=conteo_n#Cantidad de nodos predichos por plasmid que pasan el filtro del alineamiento
        #print(me)
        conteo_n=0
        #print(nodos)
        #sys.exit()
        true=ruta+'/test/True/'+lista+'/verd_mau.txt'#Los verdaderos contigs son seleccionados 
        #de un aliniamiento grafico con mauve entre spades y plasmidspades
        with open(true, mode='r') as verdaderos:
          for line in verdaderos:
            fields = line.strip()
            tp+=1
            if nodos.get(str(fields),None)!=None:
              vp+=1
        rec=(vp/me)
        #print(pre)
        pre=(vp/(tp))
        #print(rec)
        #print('nodos que pasan el filtro'+str(me))
        #print(vp)
        try:
          fsc=2*(((vp/me)*(vp/(tp)))/((vp/me)+(vp/(tp))))
        except:
          print('operacion indeterminada x/0')
        tablap[fila][2]=round(pre,4)#precision
        tablar[fila][2]=round(rec,4)#recall
        tablaf[fila][2]=round(fsc,4)#fscore
                #sacar el nombre del nodo ************************************************
      #print('cantidad plasmidos reales '+str(tp))
      #print('plasmidos predichos y reales vp '+str(vp))
      #print('nodos predichos '+str(me))
      #print('nodos spades '+str(sc))
      #tablap[itam][jtam]=round(np.mean(pre),4)#precision--------------------------------
      #tablar[itam][jtam]=round(np.mean(rec),4)#recall
      #tablaf[itam][jtam]=round(np.mean(fsc),4)#fscore
      jtam+=1
      fila+=1
    itam+=1
  print('presicion '+str(tablap))
  print('recall '+str(tablar))
  print('f-score '+str(tablaf))
  """tab=open(ruta+'/test/plasmidspades/tabla_precision557.txt','a')
  tab.write(muestras+'\n')
  for fila in tablap:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
  tab=open(ruta+'/test/plasmidspades/tabla_recall557.txt','a')
  tab.write(muestras+'\n')
  for fila in tablar:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()"""
  tab=open(ruta+'/test/plasmidspades/tabla_fscore853.txt','a')
  tab.write(muestras+'\n')
  for fila in tablaf:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
  #por cada ID,COV sacar los contigs de cada muestra
  #comparar con los True
  #evaluar el conjunto en las metricas"""
#-----------------------------Prokka------------------------------------------#
def prokka(fasta_file, out_dir,sample, prefix):
  fasta_lines = []
  """with open(fasta_file, mode='r') as fasta:
    for line in fasta:
      fasta_lines.append(line)
  with open(fasta_file, mode='w') as fasta:
    for line in fasta_lines:
      line = line.replace('_component_', '_c_')
      fasta.write(line)"""#funcion para quitar el component de los archivos multifasta resultradfo de plasmidspades
  prok = ('source /vault2/soft/miniconda2/bin/activate prokka;prokka --outdir ' + out_dir + '/Annotation/Prokka --force --cpus 0 --prefix ' + prefix +
          ' --addgenes --locustag IBUN_' + sample + ' --centre IBUN_INAS --genus Providencia ' +
          '--species rettgeri --strain ' + sample + ' --kingdom Bacteria --gcode 11 --usegenus ' +
          '--proteins /home/ebarreto/Downloads/CardDB.fasta --evalue 1e-9 --rfam ' + fasta_file)
  os.system(prok)
#-----------------------------Roary------------------------------------------#
# Run Roary
def roary(gff_list,ofolder):
  command = 'roary -e --mafft -p 8 -f ' + ofolder+ '/roary '
  for gff in gff_list:
    command = command + ' ' + gff
  os.system(command)
#--------------------------creacion archivo consolidado----------------------#
def consolidar(spades_filtrado,rutal,result_nodes):
  prob= 3  #number of test to be plasmid contig
  fasta_out = rutal + '/consolidated_plasmids.fasta'
  with open(fasta_out, 'w') as fasta_file:
      with open(spades_filtrado, 'r') as spades_fasta:
          def filtered_contigs_generator():
              for contig in SeqIO.parse(spades_fasta, 'fasta'):
                  if result_nodes.get(contig.id, None) != None and result_nodes[contig.id] >= prob:
                      print('\t' + contig.id)
                      yield contig
          SeqIO.write(filtered_contigs_generator(), fasta_file, 'fasta')
  fasta_out1 = rutal + '/consolidated_chrom.fasta'
  with open(fasta_out1, 'w') as fasta_file:
      with open(spades_filtrado, 'r') as spades_fasta:
          def filtered_contigs_generator():
              for contig in SeqIO.parse(spades_fasta, 'fasta'):
                  if result_nodes.get(contig.id, None) != None and result_nodes[contig.id] >= prob:
                      print('\t' + contig.id)
                  else:
                      yield contig
          SeqIO.write(filtered_contigs_generator(), fasta_file, 'fasta')
  return fasta_out, fasta_out1
#------------------------consolidar nodos {}---------------------------------#
def consolidate_fasta(result_nodes,nodes):
    for contig_id in nodes:
        if result_nodes.get(contig_id, None) != None:
            result_nodes[contig_id] += 1
        else:
            result_nodes[contig_id] = 1
#------------Recall para herramientas de clasificacion sobre contigs--------#
def recall(ruta,lista):
  #calificacion=ruta+'/test/calificacion_plasflow.txt'
  #calificacion=ruta+'/test/plasmidspades/calificacion_no_dust.txt'
  calificacion=ruta+'/test/calificacion_profundidad1.5.txt'
  with open(calificacion, mode ='w')as salida:
    salida.write('Muestra\tPrecision\tRecall\tF-score\n')
    for m in lista:
      if m != 'SRR9042857':
        escribir=str(m)+'\t'  
        #herramienta=ruta+'/'+m+'/Assembly/Plasmitest/Plasflow/Plasflow.txt'
        #herramienta=ruta+'/'+m+'/Assembly/Plasmitest/Plasmidspades/nodos_plasmidspades.txt'#-------------------------
        verdaderos=ruta+'/test/True/'+m+'/verdaderos.txt'
        herramienta=ruta+'/'+m+'/Assembly/Plasmitest/Profundidad/nodos_cov1.5.txt'
        #verdaderos=ruta+'/test/True/'+m+'/verd_mau.txt'
        nodos_reales={}
        predichos=0
        verd_pos=0
        real=0
        with open(verdaderos, mode ='r')as reales:
         for nodo in reales:
           num=str(((nodo.strip()).split('_'))[1])
           #num=str(nodo.strip())# mau es solo el numero
           if nodos_reales.get(num,None)==None:
             nodos_reales[num]=1
             real+=1
        print('nodos reales ')
        print(nodos_reales)
        with open(herramienta, mode='r') as clasificados:
         print('\tleyendo archivo de nodos clasificados, muestra: '+str(m))
         for clas in clasificados:
           predichos+=1
           num_clas=str(((clas.strip()).split('_'))[1])
           print(num_clas)
           if nodos_reales.get(num_clas,None)!=None:
             verd_pos+=1
        print('nodos predichos ')
        print(predichos)
        pre=verd_pos/predichos
        rec=verd_pos/real
        print(pre,rec)
        try:
          fsc=2*((pre*rec)/(pre+rec))
        except:
          fsc=0
        escribir+=str(pre)+'\t'+str(rec)+'\t'+str(fsc)+'\n'
        salida.write(escribir)
#------------------------------Filtrado homologos------------------------------#
def process_blast_csv(csv_file,subjects,tipo,salida):
  contact=''
  nodes={}
  contpas=''
  retorno={}
  with open(salida, mode='w') as salid:
    with open(csv_file, mode='r') as blast_file:
      for line in blast_file:
        line = line.strip()
        fields = line.split(',')
        node = fields[0]
        subject= fields[1]
        contact=str(fields[0])
        if contact != contpas:
          id = float(fields[2])
          cv = float(fields[10])
          #print(cv)
          cvs= (float(fields[3])/float(fields[9]))*100
          if tipo == 'p':#----------------------------
            #if id >= 90 and (cv >= 60 or cvs >= 60):
            if id >= 95 and cv >= 42.96 :
              #print('verificando contigs posibles plasmidos')
              if nodes.get(node, None) == None:
                nodes[node] = fields[1:]
                retorno[node]=1
                subjects.append(subject)
              salid.write(str(node)+'\n')
          elif tipo =='c':
            if id <= 72 and cv <= 50 :#(cv <= 60 or cvs <= 60)
              #print('verificando contigs que no pertenecen a cromosomas')
              if nodes.get(node, None) == None:
                nodes[node] = fields[1:]
                retorno[node]=1
                subjects.append(subject)
              salid.write(str(node)+'\n')
        contpas=contact
  #print(nodes, subjects)
  return retorno
#---------------------------------parametros.log extraer-----------------------#
def parmext_bd(entrada):
  print('*_________________________--------*---------______________________________*')
  print('\n\tDefinir Base de Datos a implementar, en clasificacion por homologia')
  param=open(entrada,'r')
  lineas=param.readlines()
  param.close()
  can=1
  #bp='/vault2/sgig/pipelinedbs/Plasmidos/Blastdb/Plasmidos'#prueba con base antigua de plasmidos --------------------------
  bp='/vault2/purdue/plasmid/PlasmidDBs/NCBI_rrefseq_full_plasmids/blastdb/base_plasmidos'
  bc=''
  #bc='/vault2/purdue/plasmid/PlasmidDBs/Kpneumoniae/blastdb/Kpneumoniae'#++++++++++++++++++++++++++++++++++++++++++++++-----
  #return bp,bc
  if len(lineas)>0:
    print('\t...... Ultimas Bacterias usadas : \n\n')
    for l in lineas:
      print('\t<'+str(can)+'> '+l)
      can+=1
    print('\t<'+str(can)+'> Si quiere usar otra bacteria, esta opcion correra el escript Create_Ubdate_DB.py')
    r=input('\n\tEscriba el codigo de la Bacteria que desea usar \n\t:')
    try:
      rr=int(r)
      if rr>0 and rr<can:
        bc=('/vault2/sgig/pipelinedbs/Cromosomas/'+str(lineas[rr-1])[:-1]+'/Blastdb/Cromosomas')
        return bp,bc
      elif rr==can:
        print('******************************************************************************************')
        print('********* SCRIPT CREATE UPDATE DATA BASE OF PLASMID AND CROMOSOLMAL BACTERIAL*************')
        print('******************************************************************************************')
        os.system('python /vault2/plasmid_pipeline/Scripts/Base_de_datos/Create_Update_BD.py')
      else:
        print('\n\tEl numero ingresado no esta dentro del rango')
        print('\t ----Adios----')
        sys.exit()
    except:
      print('\n\tEl dato ingresado no es un numero ..')
      print('\t ----Adios----')
      sys.exit()
  else:
    print('\tEl archivo de Bacterias esta vacio, debe correr el escript Create_Ubdate_DB.py en edirectorio /vault2/plasmid_pipeline/Scripts/Base_de_datos/')
    sys.exit()
#----------------------seleccion de contigs por homologia-------------------#
def blastn(fasta_file, out_dir,BDp,BDc,tipo):
    nodos={}
    subjects = []
    subjects1 = []
    if tipo =='c':
      #blast_csv = out_dir+'No_chromosomeDB_blastn3_olddb_out.csv'
      blast_csv = out_dir+'No_chromosomeDB_blastn3qqq.csv'
      #blast_csv = out_dir+'No_chromosomeDB_blastn_out.csv'
      #print(os.path.getsize(blast_csv))
      dir_nodos=blast_csv[:16]+'nodos_7850.txt'
      if os.path.isfile(blast_csv)==False or int(os.path.getsize(blast_csv))==0:
        print('Base de datos usada para homologia de secuencias en cromosoma: ')
        print(BDc)
        os.system('blastn -query '+fasta_file+' -db '+BDc+' -task blastn -max_hsps 3 -max_target_seqs 1 -num_threads 20 -dust no -outfmt "10 qseqid sseqid pident length qstart qend sstart send qlen slen qcovus qcovs qcovhps" -out ' + blast_csv)
      else:
        print('\talineamiento blast prebio encontrado para cromosomas')  
      if os.path.isfile(dir_nodos)==False or os.path.getsize(dir_nodos)==0:#repetir
        print('no lo haria  melo')
        #nodos=process_blast_csv(blast_csv, subjects,'c',dir_nodos)
    elif tipo =='p':
      blast_csv =out_dir+'PlasmidsDB_blastn3.csv'#prueba con base de datos antigua para comparar eficiencia y FP
      dir_nodos=blast_csv[:-11]+'nodos.txt'
      if os.path.isfile(blast_csv)==False or os.path.getsize(blast_csv)==0:
        print('Base de datos usada para homologia de secuencias en plasmido: ')
        print(BDp)
        os.system('blastn -query '+fasta_file+' -db '+BDp+' -task blastn -max_hsps 3 -max_target_seqs 1 -num_threads 20 -dust no -outfmt "10 qseqid sseqid pident length qstart qend sstart send qlen slen qcovus qcovs qcovhps" -out ' + blast_csv)
      else:
        print('\talineamiento blast prebio encontrado para plasmidos')
      if os.path.isfile(dir_nodos)==False or os.path.getsize(dir_nodos)==0:
        #nodos=process_blast_csv(blast_csv, subjects1,'p',dir_nodos)
        print('saltando plasmidos, buscando '+str(blast_csv[:-15]+'nodos.txt'))
    return nodos
#--------------------------Lectura de dictados desde dir--------------------#
def leer_dict(direc):
  nodos={}
  with open(direc, 'r') as arch:
    for linea in arch:
      contig=linea.strip()
      nodos[contig]=1
  return nodos
#--------------------Clasificacion de nodos por profundidad-----------------#
def deep_cov(spades_file,salida):
  if os.path.isdir(salida+'/Profundidad'):
      print('directory: ' + salida+'/Profundidad  EXIT!!!!')
  else:
      os.mkdir(salida+'/Profundidad')
  spade = {}  # type: Dict[Any, float]
  cov_nodes = {}
  quan=[]
  print("\n Contigs deep Coverage results:\n")
  with open(spades_file, 'rU') as input_fasta:
        for contig in SeqIO.parse(input_fasta, 'fasta'):
            line = contig.id.strip()
            fields = line.split('_')
            spade[contig.id] = float(fields[5])
            quan.append(float(fields[5]))  #all contig coverages
  q1 = np.percentile(quan, 25)
  q3 = np.percentile(quan, 75)
  limite = q3 + ((q3-q1)*1.5) #external outlier limit
  cov_file=salida+'/Profundidad/nodos_cov1.5.txt'
  print ("\n coverage " + str(limite) +"\n")
  with open(cov_file, 'w') as txt_sal:
    for contig in spade:
        if float(spade[contig]) > limite:
           cov_nodes[contig]=1
           txt_sal.write(str(contig)+'\n')
  return cov_nodes
#-----------------------------Plasflow--------------------------------------#
def plasflow(fasta_file, out_dir):
  print('activando ambiente plasflow ')
  result_nodespf = {}
  # Checck output directory
  if os.path.isdir(out_dir+'/Plasflow'):
      print('directory: ' + out_dir + '/Plasflow  EXIT!!!!')
  else:
      os.mkdir(out_dir+'/Plasflow')
  plasflow_tsv = out_dir+ '/Plasflow/Plasflow.tsv'
  plasflow_txt= out_dir+ '/Plasflow/Plasflow.txt'
  print(os.path.isfile(plasflow_txt))
  if os.path.isfile(plasflow_txt)==False:
    print('\tInicio Plasflow')
    os.system('source /vault2/soft/miniconda2/bin/activate plasflow;PlasFlow.py --input ' + fasta_file + ' --output ' + plasflow_tsv)
    #os.system('PlasFlow.py --input ' + fasta_file + ' --output ' + plasflow_tsv)
    os.system('rm -v ' + fasta_file + '*_kmer_*')
  with open(plasflow_tsv, mode='r') as plasflow_file:
    with open(plasflow_txt, mode='a') as plasflow:
      plasflow_reader = csv.reader(plasflow_file, delimiter='\t')
      print("\n Plasflow results:\n")
      for row in plasflow_reader:
        if plasflow_reader.line_num == 1: # Ignore header
          contig_id_idx = row.index('contig_name')
          label_idx = row.index('label')
          continue
        if row[label_idx].startswith('plasmid'):
          contig_id = row[contig_id_idx]
          print('\t' + contig_id)
          result_nodespf[contig_id] = 1
          plasflow.write(str(contig_id)+'\n')
  return result_nodespf
#-----------------------recall tabla por metodo------------------------------#
def table_recall(spdes_cont, t_nodes,nodes,out_dir, file_test, sample):
    j=0
    for contig_id in t_nodes:
        if nodes.get(contig_id, None) != None:
            j += 1
    recall = sample + ',' + str(j) + ','+ str(len(nodes)) + ','+ str(len(t_nodes)) + ','+ spdes_cont
    if os.path.isdir(out_dir):
        print(out_dir + '\n')
    else:
        os.mkdir(out_dir)
    os.system('echo ' + recall + '>>'+ out_dir + '/'+   file_test) 
#----------------Extraer la identidad de los contigs de spades---------------#temporal experimento
def true_nodes(spades_out,out_dir,lis):
  if os.path.isdir(out_dir)==False:
    os.mkdir(out_dir)
    print('carpeta creada con exito')
  blast_csv = out_dir+'/true_nodes_blastn3_out_spe.csv'
  if os.path.isfile(blast_csv)==False:
    blast_db = '/vault2/homehpc/dtalero/Kpneumoniae/BD_experimento/plasmidos/BD/'+lis+'/'+lis
    os.system('source activate sgig;blastn -query %s -db %s -max_target_seqs 1 -max_hsps 3 -task blastn -dust no -outfmt "10 qseqid sseqid pident qcovus length slen" -out %s'
              % (spades_out, blast_db, blast_csv))
  result_nodestn = {}
  print('extrayendo nodos verdaderos de : '+str(lis))#------------------------------------------------------
  contact=''
  contpas=''  
  with open(blast_csv, mode='r') as blast_file:
      for line in blast_file:
          line = line.strip()
          fields = line.split(',')
          id = float(fields[2])
          cv = float(fields[3])
          cvs= (float(fields[4])/float(fields[5]))*100
          contact=str(fields[0])
          if contact != contpas:
            #print('nodo '+str(fields[0]))
            #print(id)
            #print(cv)
            #print(cvs)
            if id >= 98 and (cv >= 70 or cvs >= 70):
                result_nodestn[fields[0]] = 1
                #print(fields[0])
          contpas=contact
  if os.path.isfile(out_dir+'/verdaderos.txt')==False:
    tabla=open(out_dir+'/verdaderos.txt','w')
    print('escribiendo nodos verdaderos')
    for n in result_nodestn:
      tabla.write(str(n)+'\n')
    tabla.close()
  return result_nodestn
#-----------------------Buscar nodos de spades desde plasmidspades-----------#
def blastplasmid(plasm_filtrado,spades_filtrado,out_dir):
  blast_csv=plasm_filtrado.replace('fasta','blastn3_pos_csv')
  if os.path.isfile(blast_csv)==False:
    print('alineamineto blast no encontrado ..')
    os.system("makeblastdb -in %s -dbtype 'nucl' -out %s"%(plasm_filtrado,out_dir+'/db_spades_out'))
    os.system('blastn -query %s -db %s/db_spades_out -max_target_seqs 1 -max_hsps 3 -task blastn -dust no -outfmt "10 qseqid sseqid pident qcovus length slen qlen qstart qend sstart send" -out %s'%(spades_filtrado,out_dir, blast_csv))
  result_nodesbp = {}
  contact=''
  contpas=''
  nodos_file=out_dir+'/nodos_plasmidspades.txt'
  #_no_bord2000.txt' 
  with open(nodos_file, 'w') as salida:
    with open(blast_csv, mode='r') as blast_file:
        for line in blast_file:
            line = line.strip()
            fields = line.split(',')
            id = float(fields[2])
            cv = float(fields[3])
            cvs= (float(fields[4])/float(fields[5]))*100
            qtam=float(fields[6])
            qal_ini=float(fields[7])
            qal_fin=float(fields[8])
            contact=str(fields[0])
            if contact != contpas:
              if id >= 95 and (cv >= 63 or cvs >= 63):
                #if qal_ini < 2000 or (qtam-qal_fin) < 2000:
                result_nodesbp[fields[0]] = 1
                salida.write(str(fields[0])+'\n')
            contpas=contact
  return result_nodesbp
#--------------------------------correr PlasmidSpades------------------------#
def plasmidspades(fastq1,fastq2,fastq1u,fastq2u,ks,out_dir):
  if os.path.isdir(out_dir):
    print('directory: '+out_dir+'  Encontrado!!!!')
  else:
    crearruta(out_dir)
  print('\n\tComienzo PlasmidSpades----')
  if len(fastq1u)==0 and len(fastq2u)==0:
    os.system('spades.py --plasmid --careful -m 70 -t 20 -k '+ks+' -1 '+fastq1+' -2 '+fastq2+' -o '+out_dir)
  else:
    os.system('spades.py --plasmid --careful -m 70 -t 20 -k '+ks+' -1 '+fastq1+' -2 '+fastq2+' --s1 '+fastq1u+' --s2 '+fastq2u+' -o '+out_dir)
  return (out_dir+'/contigs.fasta')   
#--------------------------------actualizar archivo Log-----------------------#
def log(ruta,texto):
  #print('\n\tActualizando archivo de ejecucion')
  dire=ruta+'/Anotation_Bacteria.log'
  if os.path.isfile(dire)==False:
    ar=open(dire,'w')
    ar.close()
  arch=open(dire,'a')
  arch.write(str(texto))
  arch.close()   
#--------------------------------filtrar y copiar salida de Spades---------------------#
def filter_spades(_file):
  min_length = 1000
  result_file = _file.replace('fa', 'filter{}.fa'.format(min_length))
  with open(result_file, 'w') as filtered_fasta:
    with open(_file, 'rU') as input_fasta:
      def filtered_contigs_generator(min):
        for contig in SeqIO.parse(input_fasta, 'fasta'):
          if len(contig) >= min:
            yield contig
          else:
            print('\tDiscarding %s with %d bases'%(contig.id, len(contig)))
      SeqIO.write(filtered_contigs_generator(int(min_length)), filtered_fasta, 'fasta')
  return result_file  
#--------------------------------crear(ruta)----------------------------------#
def crearruta(r):
  try:
    os.mkdir(r)
    print('\tCarpeta creada con exito ...\n')
  except:
    print('\tNo ha sido posible crear la ruta '+str(r))
#------------------------correr Spades3.15------------------------------------#
def spades(fastq1,fastq2,fastq1u,fastq2u,ks,out_dir):
  if os.path.isdir(out_dir):
    print('directory: '+out_dir+'  Encontrado!!!!')
  else:
    crearruta(out_dir)
  if len(fastq1u)==0 and len(fastq2u)==0:
    os.system('spades.py --careful -m 70 -t 20 -k '+ks+' -1 '+fastq1+' -2 '+fastq2+' -o '+out_dir)
  else:
    os.system('spades.py --careful -m 70 -t 20 -k '+ks+' -1 '+fastq1+' -2 '+fastq2+' --s1 '+fastq1u+' --s2 '+fastq2u+' -o '+out_dir)
  #return (out_dir+ '/spades_contigs.fasta')
  return (out_dir+'/contigs.fasta')
#------------------------Procesar Muestra-------------------------------------#
def process_file(ruta,lis,fastqs,kmer,BDc,BDp,repetir):
  #print(ruta,fastqs, kmer)
  plasm_filtrado_blast_dir='None'
  plasflow_nodes_dir='None'
  cov_nodes_dir='None'
  no_chromosome_nodes_dir='None'
  plasmid_nodes_dir='None'
  time1=datetime.now()
  log(ruta,'\n\t\n\t--------------------------------------------------------\n\t\n\t')
  log(ruta,'\n\tInicio:\n\t'+str(time1)+'\n\t')
  log(ruta,'\nIniciando proceso con la muestra:\n\t'+str(lis)+'\n\t')
  rutal=ruta+'/'+lis
  proceso=0#salto de proceso para la ejecucion de experimento
  if kmer ==0:
    kmer='21,33,55,77,99,127'
  if fastqs!=0:
    fast1 = str(fastqs[0])
    if len(str(fastqs[1]))>0:
      fast1u = str(fastqs[1])
    else:
      fast1u=''
    fast2 = str(fastqs[2])
    if len(str(fastqs[3]))>0:
      fast2u = str(fastqs[3])
    else:
      fast2u=''
  if proceso==0:
    if (os.path.isfile(rutal+'/Assembly/Spades/contigs.fasta') and os.path.getsize(rutal+'/Assembly/Spades/contigs.fasta')>0) and repetir == 0:
      print('\n\tEnsamblaje de novo previo, encontrado..')
      #log(ruta,'\tEnsamble previo encontrado..')
      spades_out=rutal+'/Assembly/Spades/contigs.fasta'
      proceso=1
    else:
      if os.path.isdir(rutal+'/Assembly'):
        print('Carpeta: '+rutal+'/Assembly'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Assembly')
      log(ruta,'\tInicio de ensamblaje de novo herramienta Spades:\n\t'+str(datetime.now()))
      spades_out=spades(fast1,fast2,fast1u,fast2u, kmer, rutal+'/Assembly/Spades')
      if os.path.isfile(rutal+'/Assembly/Spades/contigs.fasta'):
        proceso=1
    print('\n\tHa Fin proceso Spades')
    log(ruta,'\n\tFin proceso Spades')
  if proceso==1:
    if os.path.isfile(rutal+'/Assembly/Spades/contigs.filter1000.fasta') and os.path.getsize(rutal+'/Assembly/Spades/contigs.filter1000.fasta')>0 and repetir == 0:
      print('\n\tFiltrado de secuencias por tamano de ensamble de novo, encontrado..')
      log(ruta,'\tFiltrado de secuencias por tamano de ensamble de novo, encontrado..')
      spades_filtrado=rutal+'/Assembly/Spades/contigs.filter1000.fasta'
      spdes_cont = os.popen("grep -c '>' " + spades_filtrado).read()
      spdes_cont = spdes_cont.strip()
      proceso=2
    else:
      log(ruta,'\tInicio de Filtrado de secuencias:\n\t'+str(datetime.now()))
      spades_filtrado=filter_spades(spades_out)
      os.system('cp -v '+spades_filtrado+' '+rutal+'/spades_contigs.fasta ')
      spdes_cont = os.popen("grep -c '>' " + spades_filtrado).read()
      spdes_cont = spdes_cont.strip()
      log(ruta,'\tFiltrado secuencias menores a 1000 ejecutado con exito')
      if os.path.isfile(rutal+'/Assembly/Spades/contigs.filter1000.fasta'):
        proceso=2
  #Quast Saltar
  #proceso=3
  if proceso==2:
    if os.path.isfile(rutal+'/Assembly/Spades/QC/report.txt') and repetir == 0:
      print('\n\tCalculo de Metricas de ensamble de novo Quast tool, encontrado..')
      log(ruta,'\tCalculo de Metricas de ensamble de novo Quast tool, encontrado..')
      proceso=3
    else:
      if os.path.isdir(rutal+'/Assembly/Spades/QC/'):
        print('Carpeta: '+rutal+'/Assembly/Spades/QC/'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Assembly/Spades/QC/')
      log(ruta,'\tInicio de Herramienta Quast para analisis de Ensamble:\n\t'+str(datetime.now()))
      os.system('source /vault2/soft/miniconda2/bin/activate quast;quast -t 8 -o '+rutal+'/Assembly/Spades/QC/ -r /vault2/homehpc/dtalero/grupo/mishelle/pao1/sequence.fasta '+spades_filtrado)
      log(ruta,'\tFiltrado secuencias menores a 1000 ejecutado con exito')
      if os.path.isfile(rutal+'/Assembly/Spades/QC/report.txt'):
        proceso=3
    print('\n\tHa Fin proceso Metrics')
    log(ruta,'\n\tFin proceso Metrics')
  if lis == 'SRR9042857':
    print('muestra 857 no ejecutara pasos de plasmidspades ')
    proceso=6
    #Warning plasmidspades
  if proceso==3:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/warnings.log'):
      print('\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, Ha culminado sin encontrar secuencias ..')
      log(ruta,'\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, Ha culminado sin encontrar secuencias ..')
      plasm_filtrado_blast=False#-----------
      proceso=7#111111111111111111111111111111111111Salto por warning plasmidspades
    else:
      proceso=4
    print('\n\tHa Fin proceso warnings')
    log(ruta,'\n\tFin proceso warnings')
  if proceso==4:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.fasta')and os.path.getsize(rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.fasta')>0 and repetir == 0:
      print('\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, encontrado..')
      log(ruta,'\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, encontrado..')
      plasmispades_nodes=rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.fasta'
      proceso=5
    else:
      if os.path.isdir(rutal+'/Assembly/Plasmitest'):
        print('Carpeta: '+rutal+'/Assembly/Plasmitest'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Assembly/Plasmitest')
      if proceso==4:
        print('\n\tInicio de ensamble de novo de plasmidos herramienta PlasmidSpades...')
        log(ruta,'\n\tInicio de ensamble de novo de plasmidos herramienta PlasmidSpades:\n\t'+str(datetime.now()))
        plasmispades_nodes=plasmidspades(fast1,fast2,fast1u,fast2u,kmer,rutal+'/Assembly/Plasmitest/Plasmidspades')
        #Warning plasmidspades
        if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/warnings.log'):
          print('\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, Ha culminado sin encontrar secuencias ..')
          log(ruta,'\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, Ha culminado sin encontrar secuencias ..')
          plasm_filtrado_blast=False#-----------
          proceso=7#111111111111111111111111111111111111Salto por warning plasmidspades
          print('\n\tHa Fin proceso warnings')
          log(ruta,'\n\tFin proceso warnings')
        else:
          proceso=5
    print('\n\tFin proceso PlasmidSpades')
    log(ruta,'\n\tFin proceso PlasmidSpades')
  if proceso==5:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.fasta')and os.path.getsize(rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.fasta')>0 and repetir == 0:
      print('\n\tFiltrado por tamano de secuencias de plasmidos previo, encontrado..')
      log(ruta,'\n\tFiltrado por tamano de secuencias de plasmidos previo, encontrado..')
      plasm_filtrado=rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.fasta'
      proceso=6
    else:
      print('\n\tIniciando filtrado de Contigs ensamblados por PlasmidSpades')
      log(ruta,'\n\tIniciando filtrado de Contigs ensamblados por PlasmidSpades:\n\t'+str(datetime.now()))
      plasm_filtrado=filter_spades(plasmispades_nodes)
      os.system('cp -v '+plasm_filtrado+' '+rutal+'/Assembly/Plasmitest/PlasmidSpades_contigs.fasta ')
      if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.fasta'):
        proceso=6
    print('\n\tHa Fin proceso profundidad_Plasmidspades')
    log(ruta,'\n\tFin proceso profundidad_Plasmidspades')
  if proceso ==6:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/nodos_plasmidspades.txt') and os.path.getsize(rutal+'/Assembly/Plasmitest/Plasmidspades/nodos_plasmidspades.txt')>0 and repetir == 0:
      print('\n\tTabla de nodos Spades identificados por plasmidspades encontrado')
      log(ruta,'\n\tTabla de nodos Spades identificados por plasmidspades encontrado..')
      plasm_filtrado_blast_dir=rutal+'/Assembly/Plasmitest/Plasmidspades/nodos_plasmidspades.txt'
      proceso=7
    else:
      print('\n\tIniciando busqueda de contigs ensamblados por PlasmidSpades en nodos Spades:\n\t'+str(datetime.now()))
      log(ruta,'\n\tIniciando busqueda de contigs ensamblados por PlasmidSpades en nodos Spades')
      plasm_filtrado_blast=blastplasmid(plasm_filtrado,spades_filtrado,rutal+'/Assembly/Plasmitest/Plasmidspades')
      print(plasm_filtrado_blast)
      proceso=7
    print('\n\tHa Fin proceso nodos_plasmidspades')
    log(ruta,'\n\tFin proceso nodos_plasmidspades')
    log(ruta,'Ha ocurrido un error con el filtrado del archivo de plasmidSpades')
  if proceso==7:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasflow/Plasflow.txt')and os.path.getsize(rutal+'/Assembly/Plasmitest/Plasflow/Plasflow.txt')>0:
      print('\n\tTabla de contigs identificados por Plasflow encontrado')
      log(ruta,'\n\tTabla de nodos identificados por Plasflow encontrado..')
      plasflow_nodes_dir=rutal+'/Assembly/Plasmitest/Plasflow/Plasflow.txt'
      proceso=8
    else:
      print('\n\tIniciando busqueda de contigs de Plasmidos con herramienta Plasflow')
      log(ruta,'\n\tIniciando busqueda de contigs de Plasmidos con herramienta Plasflow:\n\t'+str(datetime.now()))
      plasflow_nodes = plasflow(spades_filtrado,rutal+'/Assembly/Plasmitest/')
      proceso=8
    print('\n\tHa Fin proceso Plasflow')
    log(ruta,'\n\tFin proceso Plasflow')
  if proceso==8:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Profundidad/nodos_cov1.5.txt')and os.path.getsize(rutal+'/Assembly/Plasmitest/Profundidad/nodos_cov1.5.txt')>0:#se hace evaluacion antigua 3 y con Atipico 1.5
      print('\n\tTabla de nodos identificados por profundidad, encontrado')
      log(ruta,'\n\tArchivo de contigs identificados por profundidad, encontrado..')
      cov_nodes_dir=rutal+'/Assembly/Plasmitest/Profundidad/nodos_cov.txt'#**
      proceso=9
    else:
      print('\n\tIniciando proceso de extraccion de contigs por profundidad')
      log(ruta,'\n\tIniciando proceso de extraccion de contigs por profundidad:\n\t'+str(datetime.now()))
      cov_nodes = deep_cov(spades_filtrado,rutal+'/Assembly/Plasmitest/')
      proceso=9
    print('\n\tHa Fin proceso profundidad')
    log(ruta,'\n\tFin proceso profundidad')
  #proceso=11
  #proceso 9 subdividido para crear subdirectorios de cromosoma y de plasmidos
  if proceso==9:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/No_chromosomeDB_nodos.txt') and os.path.getsize(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/No_chromosomeDB_nodos.txt')>0:
      print('\n\tTabla de nodos identificados por Homologia, encontrado')
      log(ruta,'\n\tTabla de nodos identificados porHomologia, encontrado..')
      no_chromosome_nodes_dir=rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/No_chromosomeDB_nodos.txt'
      proceso=10
    else:
      if os.path.isdir(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/'):
        print('Carpeta: '+rutal+'/Assembly/Plasmitest/Homologia'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/')
      print('\n\tIniciando proceso de extraccion de contigs por Homologia de cormosoma:')
      log(ruta,'\n\tIniciando proceso de extraccion de contigs por Homologia de cormosoma:\n\t'+str(datetime.now()))
      no_chromosome_nodes={}
      no_chromosome_nodes=blastn(spades_filtrado,rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/',BDp,BDc,'c')
      proceso=10
    print('\n\tHa Finalizado proceso Homologia Cromosoma')
    log(ruta,'\n\t Finalizado proceso Homologia Cromosoma')
  proceso=11#8888888888888888888888para evaluar solo cromosoma y luego si recorrerlo todo con sun nuevas carpetas y tral
  if proceso==10:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/PlasmidsDB_nodos.txtt')and os.path.getsize(rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/PlasmidsDB_nodos.txt')>0:#--------t--------corrida optimizado plasm
      print('\n\tTabla de nodos identificados por Homologia Plasmidos, encontrado')
      log(ruta,'\n\tTabla de nodos identificados porHomologia, encontrado..')
      plasmid_nodes_dir=rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/PlasmidsDB_nodos.txt'
      proceso=11
    else:
      if os.path.isdir(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/Plasmidos/'):
        print('Carpeta: '+rutal+'/Assembly/Plasmitest/Homologia'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/')
      print('\n\tIniciando proceso de extraccion de contigs por Homologia')
      log(ruta,'\n\tIniciando proceso de extraccion de contigs por Homologia:\n\t'+str(datetime.now()))
      plasmid_nodes={}
      plasmid_nodes=blastn(spades_filtrado,rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/',BDp,BDc,'p')
      proceso=11
    print('\n\tHa Fin proceso Homologia Plasmidos/')
    log(ruta,'\n\tFin proceso Homologia Plasmidos/')
  """if plasm_filtrado_blast_dir!='None':
    plasm_filtrado_blast=leer_dict(plasm_filtrado_blast_dir)
  if plasm_filtrado_blast==False:
    bandera_plasmidsp=1
  else:
    bandera_plasmidsp=0
  #print('nodos palsmidspades')
  #print(plasm_filtrado_blast)
  if plasflow_nodes_dir!='None':
    plasflow_nodes=leer_dict(plasflow_nodes_dir)
  #print('nodos plasflow')
  #print(plasflow_nodes)
  if cov_nodes_dir!='None':
    cov_nodes=leer_dict(cov_nodes_dir)
  #print('nodos covertura')
  #print(cov_nodes)
  if no_chromosome_nodes_dir!='None':
    no_chromosome_nodes=leer_dict(no_chromosome_nodes_dir)
  #print('nodos cromo homologos')
  #print(no_chromosome_nodes)
  if plasmid_nodes_dir!='None':
    plasmid_nodes=leer_dict(plasmid_nodes_dir)  
  #print('nodos plasmidos homologos')
  #print(plasmid_nodes) """
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
  """result_nodes = {}
  for i in (plasm_filtrado_blast,plasflow_nodes,cov_nodes,plasmid_nodes,no_chromosome_nodes):
    #print(i)
    consolidate_fasta(result_nodes, i)
  if proceso==11:
    if os.path.isfile(rutal+'/consolidated_plasmids.fasta') and os.path.isfile(rutal+'/consolidated_chrom.fasta')and os.path.getsize(rutal+'/consolidated_plasmids.fasta')>0 and os.path.getsize(rutal+'/consolidated_chrom.fasta')>0:
      print('\n\tTabla de contigs consolidados, encontrada')
      log(ruta,'\n\tTabla de contigs consolidados, encontrada..')
      tablaP_dir=rutal+'/consolidated_plasmids.fasta'
      tablaC_dir=rutal+'/consolidated_chrom.fasta'
      proceso=12
    else:
      print('\n\tIniciando proceso de generacion tabla de consolidados')
      log(ruta,'\n\tIniciando proceso de generacion tabla de consolidados:\n\t'+str(datetime.now()))
      tablaP_dir, tablaC_dir=consolidar(spades_filtrado,rutal,result_nodes)
      proceso=12
  else:
    print('\n\tHa ocurrido un error con la busqueda de homologos')
    log(ruta,'Ha ocurrido un error con la busqueda de homologos')"""
  
  """if proceso==12:
    #if os.path.isfile(rutal+'/Annotation/Prokka/'+lis+'_P.gbk') and os.path.isfile(rutal+'/Annotation/Prokka/'+lis+'_C.gbk') and os.path.isfile(rutal+'/Annotation/Prokka/'+lis+'_T.gbk') and os.path.getsize(rutal+'/Annotation/Prokka/'+lis+'_P.gbk')>0 and os.path.getsize(rutal+'/Annotation/Prokka/'+lis+'_C.gbk')>0 and os.path.getsize(rutal+'/Annotation/Prokka/'+lis+'_T.gbk')>0:
    if os.path.isfile(rutal+'/Annotation/Prokka/'+lis+'_T.gbk') and os.path.getsize(rutal+'/Annotation/Prokka/'+lis+'_T.gbk')>0:
      print('\n\tAnotacion, encontrada')
      log(ruta,'\n\tAnotacion, encontrada..')
      proceso=13
    else:
      print('\n\tIniciando proceso de Anotacion')
      if os.path.isdir(rutal+'/Annotation/Prokka/'):
        print('Carpeta: '+rutal+'/Annotation/Prokka/'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Annotation/Prokka/')
      log(ruta,'\n\tIniciando proceso de Anotacion:\n\t'+str(datetime.now()))
      #if os.path.isfile(rutal+'/Annotation/Prokka/'+lis+'_P.gbk') ==False:
      #  prokka(tablaP_dir, rutal,lis,lis+'_P')
      #if os.path.isfile(rutal+'/Annotation/Prokka/'+lis+'_C.gbk') ==False:
      #  prokka(tablaC_dir,rutal,lis,lis+'_C')
      if os.path.isfile(rutal+'/Annotation/Prokka/'+lis+'_T.gbk')==False:
        prokka(spades_filtrado,rutal,lis,lis+'_T')
      proceso=13
  else:
    print('\n\tHa ocurrido un error con la generacion de consolidados')
    log(ruta,'Ha ocurrido un error con la generacion de consolidados')"""
  time2=datetime.now()
  log(ruta,'Finalizado, el procesamiento de la muestra:\n\t'+str(lis))
  log(ruta,'Tiempo total:\n\t'+str(time2-time1)+'\n\t')
#---------------------------Enviar correo Electronico--------------------------#
def enviarnot(correo):
  msg = MIMEMultipart() 
  message = "\n\tSu proceso de creacion/actualizacion de la base de datos ha culminado con exito"
  message+='\n\tVerificar en la carpeta de destino el archivo info.log y la existencia de la base de datos blast'
  message+='\n\tGracias por implementar Create_UpdateBD.py \n Exito en sus investigaciones...'
  # setup the parameters of the message
  password = "Taler*1019"
  msg['From'] = "dtalero@unal.edu.co"
  msg['To'] = str(correo)
  msg['Subject'] = "Notificacion  Create_UpdateBD"
  msg.attach(MIMEText(message, 'plain'))
  #create server
  server = smtplib.SMTP('smtp.gmail.com: 587') 
  server.starttls()
  # Login Credentials for sending the mail
  server.login(msg['From'], password) 
  # send the message via the server.
  server.sendmail(msg['From'], msg['To'], msg.as_string())
  server.quit()
  print("\Correo enviado con exito %s:" % (msg['To']))
#------------------------Extraer correo de usuario-----------------------------#
def pedcorreo():
  correo='n'
  print('*_________________________--------*---------______________________________*')
  noti=str(input('\n\tDesea ser notificado por correo al terminar la operacion ? \n\n\t<1> SI \n\t<2> NO \n\tluego la tecla enter \n\t:'))
  if noti=='1':
    correo=str(input('\n\tIngrese su correo de gmail porfavor :'))
    return correo
  elif noti=='2':
    correo='n'
    return correo
  else: 
    print('\n\tEl dato ingresado no es valido ..')
    print('\t ----Adios----')
    sys.exit()
#-------------------------Extraer fastq del directorio-------------------------#
def fastqextract(ruta,muestra,indice,p1,p2,s1,s2):
  fastqs=['','','','']
  ensamble=0
  print(str(type(muestra)))
  if str(type(muestra))=='<class \'str\'>' or str(type(muestra))=='<type \'str\'>':
    direcc=ruta+'/'+muestra+'/'+indice
  else:
    direcc=ruta+'/'+str(muestra)[2:-2]+'/'+indice
  for archivo in os.listdir(direcc):
    if archivo[-len(p1):]==p1:#'fastq.gz'
      fastqs[0]=direcc+'/'+archivo
      ensamble=1
    elif archivo[-len(s1):]==s1:#'fastq.gz'
      fastqs[1]=direcc+'/'+archivo
      ensamble=1
    elif archivo[-len(p2):]==p2:#'fastq.gz'
      fastqs[2]=direcc+'/'+archivo
      ensamble=1
    elif archivo[-len(s2):]==s2:#'fastq.gz'
      fastqs[3]=direcc+'/'+archivo
      ensamble=1
  print('\tLecturas encontradas : '+fastqs[0])
  #print(fastqs)
  if ensamble ==1:
    return fastqs
  else:
    print('\t\nOcurrio un error encontrando las lecturas para el ensamble')
    sys.exit()
    #print(os.path.isfile(direcc+'/'+archivo))
  print(fastqs)
#------------------------------Funcion principal de pipeline-------------------#
def pipeline(ruta,lista,indice,kmer,BDc,BDp,p1,p2,s1,s2,repetir):
  str(type(lista))
  if str(type(lista))!='<class \'str\'>' and  str(type(lista))!='<type \'str\'>':#mientras abacas626,668
    for lis in lista:
      fastqs=fastqextract(ruta,lis,indice,p1,p2,s1,s2)
      if fastqs!=0:
        process_file(ruta,lis,fastqs,kmer,BDc,BDp,repetir)
      else:
        print('No se encontraron lecturas para la muestra '+str(lis)) 
  else:
    fastqs=fastqextract(ruta,lista,indice,p1,p2,s1,s2)
    if fastqs!=0:
      process_file(ruta,lista,fastqs, kmer,BDc,BDp,repetir)
    else:
      print('No se encontraron lecturas para la muestra '+str(lis))
  #print('Ejecutando calculo de reporte de profundidad de Ensamble ')#{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{
  #reportQC(ruta,lista)
#------------------------------Establecer kmers-.------------------------------#
def pedkmers():
  rep_ens=0
  print('*_________________________--------*---------______________________________*')
  print('\n\tEn referencia a los Ensambles de novo: ?')
  rep=str(input('\n\t<1> Desea usar Ensambles previos \n\t<2> Generar nuevos \n\tluego la tecla enter \n\t:')) 
  if rep=='1':
    rep_ens=0
    print('*_________________________--------*---------______________________________*')
    print('\n\tEn caso de que una o varias muestras no cuenten con un Ensamble de novo')
    print('\n\tDesea usar los kmers por defecto de la herramienta Spades ?\n')
    r=str(input('\n\t<1> SI \n\t<2> NO \n\tluego la tecla enter \n\t:'))
    if r=='1':
      return 0,rep_ens
    elif r=='2':
      print('*_________________________--------*---------______________________________*')
      print('\n\tPara reads con tamanios mayor a 250 pb es recomendable uzar 6 tamanos de kmers ')
      kmers=str(input('\n\tPorfavor ingerese los kmer separados por coma <#,#,#,#,#,#> luego la tecla enter \n\t:'))
      if ',' in kmers:
        return kmers,rep_ens
      else:
        print('\n\tEl dato ingresado no es valido ..')
        print('\t ----Adios----')
        sys.exit()
    else: 
      print('\n\tEl dato ingresado no es valido ..')
      print('\t ----Adios----')
      sys.exit()
  elif rep=='2':
    rep_ens=1
    print('*_________________________--------*---------______________________________*')
    print('\n\tDesea usar los kmers por defecto de la herramienta Spades ?\n')
    r=str(input('\n\t<1> SI \n\t<2> NO \n\tluego la tecla enter \n\t:'))
    if r=='1':
      return 0,rep_ens
    elif r=='2':
      print('*_________________________--------*---------______________________________*')
      print('\n\tPara reads con tamanios mayor a 250 pb es recomendable uzar 6 tamanos de kmers ')
      kmers=str(input('\n\tPorfavor ingerese los kmer separados por coma <#,#,#,#,#,#> luego la tecla enter \n\t:'))
      if ',' in kmers:
        return kmers,rep_ens
      else:
        print('\n\tEl dato ingresado no es valido ..')
        print('\t ----Adios----')
        sys.exit()
    else: 
      print('\n\tEl dato ingresado no es valido ..')
      print('\t ----Adios----')
      sys.exit()
  else: 
    print('\n\tEl dato ingresado no es valido ..')
    print('\t ----Adios----')
    sys.exit()
#--------------funcion para la administracion de formatos de lecturas----------#
def funcion_lecturas(ruta,lista,indice):
  print('*_________________________--------*---------______________________________*')
  print('\n\tEstablecer las lecturas para el ensamble de novo.')
  print('\n\tDesea usar los formatos por defecto:')
  print('\n\tR1_001.paired.fastq.gz,R1_001.unpaired.fastq.gz,R2_001.paired.fastq.gz,R2_001.unpaired.fastq.gz')
  r_def=str(input('\n\t<1> SI \n\t<2> NO \n\tluego la tecla enter \n\t:'))
  if r_def=='1':
    #p1,p2,s1,s2='R1_001.paired.fastq.gz','R2_001.paired.fastq.gz','R1_001.unpaired.fastq.gz','R2_001.unpaired.fastq.gz'
    p1,p2,s1,s2='1.p.fastq','2.p.fastq','1.u.fastq','2.u.fastq'
    return p1,p2,s1,s2
  elif r_def=='2':
    p1,p2,s1,s2=ped_lect(ruta,lista,indice)
    return p1,p2,s1,s2
  else:
    print('Seleccion incorrecta')
    sys.exit()
#------------------------------Establecer lecturas para ensamble---------------#
def ped_lect(ruta,listaa,indice):
  print('*_________________________--------*---------______________________________*')
  print('\n\tEjemplo de seleccion de formato:\n\tSi el formato de sus lecturas fuese el siguiente:\n\t<1>muestrafdt.r1.paired.fastq.gz/')
  print('\t<2>muestrafdt.r2.paired.fastq.gz/\n\t<3>muestrafdt.r1.unpaired.fastq.gz/\n\t<4>muestrafdt.r2.unpaired.fastq.gz/')
  print('\tUsted debe colocar en el formato de la ruta pareada derecha(Forward): r1.paired.fastq.gz')
  print('\ten el formato de la ruta pareada izquierda(Reverse): r2.paired.fastq.gz')
  print('\ty de tener lecturas no pareadas debe colocar tanto para lectura derecha como la izquierda: r#.unpaired.fastq.gz')
  lecturas={}
  lista=listaa
  print(type(lista))
  if str(type(lista))!='<class \'str\'>' and  str(type(lista))!='<type \'str\'>':
    for lis in lista[:1]:
      carp=os.listdir(ruta+'/'+lis+'/'+indice)
      for sub in carp:
        if sub.find('.') != -1:
         if lecturas.get(str(sub))==None:
           lecturas[str(sub)]=1
         else:
           lecturas[str(sub)]=lecturas[str(sub)]+1

  else:
    carp=os.listdir(ruta+'/'+lista+'/'+indice)
    for sub in carp:
      if sub.find('.') != -1:
        if lecturas.get(str(sub))==None:
          lecturas[str(sub)]=1
        else:
          lecturas[str(sub)]=lecturas[str(sub)]+1
  print('\n\tSe encontraron las siguientes lecturas:\n')
  conteo=1
  arreglo=[]
  for key in lecturas:
    print('\t<'+str(conteo)+'> '+key+'/')
    conteo+=1
  p1=(input('\n\tEscriba el formato de la ruta pareada derecha, Ejemplo:  "R1.fastq" \n\t:'))
  
  p2=(input('\n\tEscriba el formato de la ruta pareada izquierda, Ejemplo:  "R2.fastq" \n\t:'))
  
  ss1=(input('\n\tEscriba el formato de la ruta no pareada derecha, Ejemplo:  "R1.u.fastq" \n\tSi no tiene Una oprima n\n\t:'))
  if ss1=='n':
    s1=''
  else:
    s1=ss1
  ss2=(input('\n\tEscriba el formato de la ruta no pareada izquierda, Ejemplo:  "R2.u.fastq" \n\tSi no tiene Una oprima n\n\t:'))
  if ss2=='n':
    s2=''
  else:
    s2=ss2
  
  return p1,p2,s1,s2
#-------------------------------Establecer indice------------------------------#
def pedindice(ruta,listaa):
  print('*_________________________--------*---------______________________________*')
  print('\n\tEstablecer la ruta Base de los Reads para el ensamble \n\tEjemplo: <Reads_Clean>....')
  carpetas={}
  lista=listaa
  print(type(lista))
  if str(type(lista))!='<class \'str\'>' and  str(type(lista))!='<type \'str\'>':
    for lis in lista:
      subcar=os.listdir(ruta+'/'+lis)
      for sub in subcar:
        if carpetas.get(str(sub))==None:
          carpetas[str(sub)]=1
        else:
          carpetas[str(sub)]=carpetas[str(sub)]+1
  else:
    subcar=os.listdir(ruta+'/'+lista)
    for sub in subcar:
      if carpetas.get(str(sub))==None:
        carpetas[str(sub)]=1
      else:
        carpetas[str(sub)]=carpetas[str(sub)]+1
  if carpetas.get('Reads_Clean',None)!=None:
    print('\n\tCarpeta Reads_Clean/ encontrada... ')
    print('\n\tDesea usar Reads_Clean/ como folder de las lecturas de Ensamble de novo? \n\t:')
    r=str(input('\n\t<1> SI \n\t<2> NO \n\tluego la tecla enter \n\t:'))
    if r=='1':
      return 'Reads_Clean'
    elif r=='2':
      print('\n\tSe encontraron las siguientes carpetas recurrentes:\n')
      conteo=1
      arreglo=[]
      for key in carpetas:
        print('\t<'+str(conteo)+'> '+key+'/')
        arreglo.append(key)
        conteo+=1
      r=int(input('\n\tEscriba el codigo de la ruta que desea usar \n\t:'))
      try:
        rr=int(r)
        if rr>0 and rr<conteo:
          return str(arreglo[r-1])
        else:
          print('\n\tEl numero ingresado no esta dentro del rango')
          print('\t ----Adios----')
          sys.exit()
      except:
          print('\n\tEl dato ingresado no es un numero ..')
          print('\t ----Adios----')
          sys.exit()
  else:
    print('\n\tSe encontraron las siguientes carpetas recurrentes:\n')
    conteo=1
    arreglo=[]
    for key in carpetas:
      print('\t<'+str(conteo)+'> '+key+'/')
      arreglo.append(key)
      conteo+=1
    r=int(input('\n\tEscriba el codigo de la ruta que desea usar \n\t:'))
    try:
      rr=int(r)
      if rr>0 and rr<conteo:
        return str(arreglo[r-1])
      else:
        print('\n\tEl numero ingresado no esta dentro del rango')
        print('\t ----Adios----')
        sys.exit()
    except:
        print('\n\tEl dato ingresado no es un numero ..')
        print('\t ----Adios----')
        sys.exit()
#----------------------------Seleccion especifica de muestras-------------------#
def especific(ruta):
  listas=os.listdir(ruta)
  lista=[]
  for l in listas:
    if l.find('.')==-1:
      #if l!='test' and l!='SRR9042857'and l!='SRR3465532':#filtro para evaluaciion de muestras del experimento plasmidspa
      if l!='test' :#and l!='SRR9042857' and l!='SRR3465532' and l!='SRR8668707':#para prof
        lista.append(l)
  conteo=1
  print('*_________________________--------*---------______________________________*')
  print('\n\tMuestras Listadas en el directorio de entrada :\n')
  for lis in lista:
    print('\t<'+str(conteo)+'> '+lis+'/')
    conteo+=1
  print('\n\tEscriba la muestra, la lista o el rango de muestras que desea correr\n\tEjemplo:\n\t10\n\t1,5,10,23\n\t5-8')
  r=str(input('\n\tEscriba su respuesta y luego la tecla Enter\n\t:'))
  lism=[]
  if ',' in r:
    muestras=r.split(',')
    print('\n\tLa lista de muestras es:\n')
    for mu in muestras:
      print(mu)
      print('\t<'+mu+'>  '+str(lista[int(mu)-1]))
      lism.append(str(lista[int(mu)-1]))
    return lism
  elif '-' in r:
    muestras=r.split('-')
    lism=lista[int(muestras[0])-1:int(muestras[1])]
    print('\n\tLa lista de muestras seleccionadas es:\n')
    ini=int(muestras[0])
    for li in lism:
      print('\t<'+str(ini)+'>  '+str(li))
      ini+=1
    return lism
  else:
    try:
      uni=int(r)
      print('\tMuestra seleccionada :'+str(lista[uni-1]))
      return lista[uni-1]
    except:
      print('\n\tEl dato ingresado no es una opcion valida ..')
      print('\t ----Adios----')
      sys.exit()
#-------------------------------Establecer lista de muestras-------------------#
def pedmuestras(ruta):
  print('*_________________________--------*---------______________________________*')
  print('\n\tBuscando directorios referentes a muestras ....')
  print('\t'+ruta)
  listas=os.listdir(ruta)
  lista=[]
  for l in listas:
    if l.find('.')==-1:
      #if l!='SRR9042857'and l!='test'and l!='SRR3465532':#filtro para evaluaciion de muestras del experimento plasmidspades
      if l!='test' and l!='RGI' and l!='roary' and l!='deep' and l!='QC'and l!='SRR9042857' and l!='SRR3465532' and l!='SRR8668707':#evaluacion prof 
        lista.append(l)
      else:
        print('')
  print('\tEn el directorio seleccionado se encuentran '+str(len(lista))+' muestras\n\tDesea correr el Script en:\n')
  print('\t<1> Todas \n\t<2> Una muestra o una lista especifica')
  try:
    r=int(input('\n\tEscriba el codigo de la opcion \n\t:'))
    rr=r
    if rr>0 and rr<3:
      if rr==1:
        return lista
      elif rr==2:
        return especific(ruta)
    else:
      print('\n\tEl numero ingresado no esta dentro del rango')
      print('\t ----Adios----')
      sys.exit()
  except:
    print('\n\tEl dato ingresado no es un numero .u.')
    print('\t ----Adios----')
    sys.exit()
#--------------------------------pedir ruta------------------------------------#
def pedruta():
  print('\n\tEscriba la ruta de las lecturas ..\nPorfavor escriba la ruta absoluta /../../Ejemplo\n\t:')
  r=input()
  if os.path.isdir(r)==True:
    print('Ruta existente ...')
    return r
  else:
    print('\n\tLa ruta ingresada no existe .. \n')
    sys.exit()
#---------------------------------estableciendo ruta---------------------------#
def pruta():
  print('*_________________________--------*---------______________________________*')
  print('\n\tEstablecer la ruta principal de las Lecturas ....')
  ruta=raiz
  dr=''
  print('\tArchivo de directorios :'+str(ruta)+'/.params.log')
  if os.path.isfile(str(ruta)+'/.params.log')==True:
      rutalecs=str(parmext(str(raiz)+'/.params.log'))
      if rutalecs=='0':
        print('\n\t Nuevo ingreso seleccionado ..\n')
        r= pedruta()
        parmact(r)
        return r
      else:
        return rutalecs
  else:
    r=pedruta()
    a=open(ruta+'/.params.log','w')
    a.close()
    parmact(r)
    return r
#---------------------------------parametros.log actualizar--------------------#
def parmact(ruta):
  param=open(raiz+'/.params.log','r')
  lin = param.readlines()
  param.close()
  esta=0
  for li in lin:
    if li[:-1]==ruta:
      esta=1
      print('\n\tRuta ingresada ya residia en la lista de Bases de Datos Usadas ..')
  if esta ==0:
    param=open(raiz+'/.params.log','a')
    param.write(str(ruta)+'\n')
    param.close()
    print('\tArchivo de parametros de usuario actualizado ..\n')
  guar=[]
  if len(lin)>15:#------
    guar=lin[-15:]
    param=open(raiz+'/.params.log','w')
    param.writelines(guar)
    param.close()
#---------------------------------parametros.log extraer-----------------------#
def parmext(entrada):
  param=open(entrada,'r')
  lineas=param.readlines()
  param.close()
  can=1
  if len(lineas)>0:
    print('\t...... Ultimas usadas : \n\n')
    for l in lineas:
      print('\t<'+str(can)+'> '+l)
      can+=1
    print('\t<'+str(can)+'> NUEVO INGRESO\n')
    r=input('\n\tEscriba el codigo de la ruta que desea usar \n\t:')
    try:
      rr=int(r)
      if rr>0 and rr<can:
        return str(lineas[rr-1])[:-1]
      elif rr==can:
        return 0
      else:
        print('\n\tEl numero ingresado no esta dentro del rango')
        print('\t ----Adios----')
        sys.exit()
    except:
      print('\n\tEl dato ingresado no es un numero ..')
      print('\t ----Adios----')
      sys.exit()
  else:
    print('\tEl archivo esta vacio')
    return 0
#---crear archivo de direcciones para background-----------#
def crearlistam(ruta,lista):
  print('\n\tCreando archivo de lista de las muestras a procesar ..')
  dire=ruta+'/lista_muestras_back.txt'
  if str(type(lista))[-3:-2]=='t':
    with open(dire, 'w') as guardar:
      for li in lista:
        guardar.write(ruta+'/'+str(li)+'\n')
  else:
    with open(dire, 'w') as guardar:
      guardar.write(ruta+'/'+str(lista)+'\n')
  return dire
#----------------modo interactivo--------------------------#
def interactivo(raiz):
  ruta=pruta()#Establecer la ruta principal
  lista=pedmuestras(ruta)#Establecer muestras a implementar
  print(lista)
  indice=pedindice(ruta,lista)#Establecer la ruta Base de los Reads
  print(indice)
  p1,p2,s1,s2=funcion_lecturas(ruta,lista,indice)
  print(p1,p2,s1,s2)
  kmer=0
  kmer,repetir=pedkmers()#configurar kmers o repetir ensamble
  correo=pedcorreo()#Establecer correo electronico de notofocacion
  dirbds='/vault2/sgig/pipelinedbs/Cromosomas/Bacterias.txt'
  BDp,BDc=parmext_bd(dirbds)#Definir Base de Datos
  #tablanodos(ruta,lista)
  #experimento3(ruta,lista,indice)
  #explorprof(ruta,lista)
  #explorprofseq(ruta,lista)
  #recallprof(ruta,lista)
  #recallplasmid(ruta,lista)
  #recallhomo(ruta,lista)
  #recall(ruta,lista)
  #recallplasmid_puntas(ruta,lista)
  #recallhomocro(ruta,lista)
  #tabladiv(ruta,lista)
  #reportQC(ruta,lista)
  #sys.exit()
  print('*_________________________--------*---------______________________________*')
  print('\n\tDesea correr el script en background?')
  back=str(input('\n\t<1> SI \n\t<2> NO \n\tluego la tecla enter \n\t:'))
  if back=='1':
    dir_lis=crearlistam(ruta,lista)
    print('***********************************************************************************************************')
    print('*\tIniciando PipeLine de Ensamble, Identificacion de plasmidos y anotacion de genoma en backgorund..\t* ')
    print('***********************************************************************************************************')
    crom=BDc.split('/')
    cromo=str(crom[5])
    os.system('nohup python '+raiz+'/An_Bac_back.py -r '+dir_lis+' -c '+cromo+' -d '+indice+' -k '+str(kmer)+' -1 '+p1+' -2 '+p2+' -3 '+s1+' -4 '+s2+' &')
  elif back=='2':
    print('****************************************************************************************')
    print('\tIniciando PipeLine de Ensamble, Identificacion de plasmidos y anotacion de genoma .. ')
    print('****************************************************************************************')
    pipeline(ruta,lista,indice,kmer,BDc,BDp,p1,p2,s1,s2,repetir)
    #recall(ruta,lista)
  else:
    print('Seleccion incorrecta')
    sys.exit()
  if correo != 'n':
    enviarnot(correo)
  sys.exit()
#------------------------extraer lista de muestras----------#
def sacar_lista(archivo):
  lista=[]
  with open(archivo,'r') as arch:
    for linea in arch:
      #print(linea)
      lin=(linea.strip()).split('/')
      lista.append(str(lin[-1:])[2:-2])
      ruta=(linea)[:-(len(str(lin[-1:])[:-3]))]
    #print(lista,ruta)
  return ruta,lista
#-----------------------linea de comando--------------------#
def linea_comando(r,c,k,d,p1,p2,s1,s2,n):
  if len(r)!=0:
    ruta,lista=sacar_lista(r)
    #print(ruta,lista)
  print('\t\nCarpeta de lecturas ')
  print('\t\n'+d)
  if len(k)>1:
    kmers=k
    print('\t\nkmers seleccionados :')
    print(k)
  else:
    print('\t\nSe usaran los kmers por defecto para el ensamble de novo')
    kmers=0
  BDp='/vault2/sgig/pipelinedbs/Plasmidos/Blastdb/Plasmidos'
  BDc=('/vault2/sgig/pipelinedbs/Cromosomas/'+str(c)+'/Blastdb/Cromosomas')
  pipeline(ruta,lista,d,kmers,BDc,BDp,p1,p2,s1,s2)
  if n != 'n':
    enviarnot(n)
#-----------------------manual de usuario e instructivo-----# 
def help_message(raiz):
  print('\n\tBuscando manual de usuario .....:')
  print(str(raiz)+'/User_Manual.txt')
  ruta=raiz
  if os.path.isfile(str(ruta)+'/User_Manual.txt')==True:
    man=open(raiz+'/User_Manual.txt','r')
    lineas=man.readlines()
    man.close()
    for l in lineas[54:129]:
      print(l[:-1])
  else:
    print('\n\tNo se encuentra el Manual de usuario en la carpeta ')
def manual():
  print('\n\tBuscando manual de usuario .....:')
  print(raiz)
  ruta=raiz
  if os.path.isfile(str(ruta)+'/User_Manual.txt')==True:
    man=open(raiz+'/User_Manual.txt','r')
    lineas=man.readlines()
    man.close()
    for l in lineas:
      print(l[:-1])
  else:
    print('\n\tNo se encuentra el Manual de usuario en la carpeta ')
def main(a,raiz):# funcion si encuentra argumentos de usuario -t, -b ...
  r,k,c,d,n,h,m=0,0,0,0,0,0,0,
  correo='n'
  ruta=''
  bacteria=''
  indice_lec='Reads_Clean'
  p1='R1_001.paired.fastq.gz'
  s1='R1_001.unpaired.fastq.gz'
  p2='R2_001.paired.fastq.gz'
  s2='R2_001.unpaired.fastq.gz'
  b1,b2,b3,b4=0,0,0,0
  ss1=''
  ss2=''
  #p1=''
  #s1=''
  #p2=''
  #s2=''
  kmers=''
  try:
    opts, args = getopt.getopt(a, "hmr:c:k:d:n:1:2:3:4:")
    #print(opts)
  except getopt.GetoptError:
    help_message()
    print('\n\tError de argumentos de compilacion')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':    
      help_message(raiz)
      sys.exit()
    elif opt == '-m':    
     manual()
     sys.exit()
    elif opt == "-r":
      ruta= os.path.realpath(arg)#elimina los caracteres simbolicos de la direccion
      if os.path.isfile(ruta):
        print('\tArchivo de muestras a procesar:  '+str(ruta))
        r=1
      else:
        print('\tLa ruta del archivo no existe '+str(ruta)+'  no existe ......')
        sys.exit()
    elif opt=='-c':
      bacteria=str(arg)#.replace(' ','_')
      c=1
    elif opt == "-d": 
      indice_lec=str(arg)
      d=1
    elif opt=='-n':
      correo=str(arg)
      n=1
    elif opt=='-k':
      kmers=str(arg)
      k=1
    elif opt=='-1':
      p1=str(arg)
      b1=1
    elif opt=='-2':
      p2=str(arg)
      b2=1
    elif opt=='-3':
      ss1=str(arg)
      b3=1
    elif opt=='-4':
      ss2=str(arg)
      b4=1
    else:
      print(opt)
      print('\tError en el ingreso de parametros ..')
      help_message(raiz)
      sys.exit()
   
  if r==0 and c==0:
    print('\n\tDebe minimo ingresar un archivo con la lista de muestras a procesar parametro <-r> y el nombre de la bacteria parametro <-c>')
    sys.exit()
  else:
    if b3==0 and b4==0:
      s1=''
      s2=''
    elif b3==1 and b4==1:
      s1=ss1
      s2=ss2
    linea_comando(ruta,bacteria,kmers,indice_lec,p1,p2,s1,s2,correo)      
if __name__ == "__main__":
  print('Parametros :')
  print(sys.argv)
  raiz=str(sys.argv[0])[:-22]
  if len(raiz)==0:
      raiz='.'
  if len(sys.argv[1:]) != 0:
    main(sys.argv[1:],raiz)
  else:
    interactivo(raiz)
