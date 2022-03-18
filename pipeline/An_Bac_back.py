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
try:
    import joblib
except ImportError:
    print("libreria joblib no instalada")
try:
    import traceback
except ImportError:
    print("libreria traceback no instalada")
    
#--------------------------------unificacion de resultados---------------------
def multi_modelo(vector):
  modelo=joblib.load('/vault2/plasmid_pipeline/Scripts/Pipe_plasmids/clasificador5m.joblib')
  #entrada=
  print('entra a evaluar')
  arreglo=np.array([vector])
  salida=modelo.predict(arreglo)
  if salida[0]==1:
    return  True
  else:
    return  False
  #print(arreglo,salida)
#--------------reporte de profundidad en base a reads clean y ensamble    
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
#------------------------------filtrar y borrar residuo de alineamiento-------#
def filt_borr(arch,arch_n):
  with open(arch, 'r') as previo:
    with open(arch_n, 'w') as nuevo:
      contact=''
      contpas=''
      for lin in previo:
        linea=(lin.strip()).split(',')
        contact=linea[0]  
        if contact!=contpas:
          nuevo.write(lin)
        contpas=contact
  print('archivo filtrado creado, revisar y camellarle al borrado y corrida de TODO, va bien :*')
  os.remove(arch)
  print('Archivo de alineamiento filtrado y borrado')
#-----------------------------Prokka------------------------------------------#
def prokka(fasta_file, out_dir,sample, prefix):
  fasta_lines = []
  if prefix[-2:]=='_P':
    prok = ('source /vault2/soft/miniconda2/bin/activate prokka;prokka --outdir ' + out_dir + '/Annotation/Prokka/plasmids --force --cpus 0 --prefix ' + prefix +
            ' --addgenes --locustag IBUN_' + sample + ' --centre IBUN_INAS --genus Providencia ' +
            '--species rettgeri --strain ' + sample + ' --kingdom Bacteria --gcode 11 --usegenus ' +
            '--proteins /home/ebarreto/Downloads/CardDB.fasta --evalue 1e-9 --rfam ' + fasta_file)
    os.system(prok)
  elif prefix[-2:]=='_C':
    prok = ('source /vault2/soft/miniconda2/bin/activate prokka;prokka --outdir ' + out_dir + '/Annotation/Prokka/cromosome --force --cpus 0 --prefix ' + prefix +
            ' --addgenes --locustag IBUN_' + sample + ' --centre IBUN_INAS --genus Providencia ' +
            '--species rettgeri --strain ' + sample + ' --kingdom Bacteria --gcode 11 --usegenus ' +
            '--proteins /home/ebarreto/Downloads/CardDB.fasta --evalue 1e-9 --rfam ' + fasta_file)
    os.system(prok)
  elif prefix[-2:]=='_T':
    prok = ('source /vault2/soft/miniconda2/bin/activate prokka;prokka --outdir ' + out_dir + '/Annotation/Prokka/genome --force --cpus 0 --prefix ' + prefix +
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
def consolidar(spades_filtrado,rutal,result_nodes,etiq):
  dictado={}
  
  #print(result_nodes.keys())
  #dictado=result_nodes['plasm_filtrado_blast']
  #print(dictado.keys())
  #print(';P')
  #print(dictado)
  fasta_plas = rutal + '/consolidated_plasmids.fasta'
  fasta_crom = rutal + '/consolidated_chrom.fasta'
  with open(fasta_crom, 'w') as crom:
    with open(fasta_plas, 'w') as plasm:
      with open(spades_filtrado, 'r') as spades_fasta:
        for contig in SeqIO.parse(spades_fasta, 'fasta'):
          evaluar=[0,0,0,0,0]
          posi=0
          for i in etiq:#----------------utilizando todos
            dictado={}
            dictado=result_nodes[i]
            if dictado.get(contig.id, None):
                #print('\t' + contig.id)
                evaluar[posi]=1
                #
            posi+=1
          if multi_modelo(evaluar):
            line_fast='>'+str(contig.id)+'\n'+str(contig.seq)+'\n'
            plasm.write(line_fast)
            print('uno y  otro')
          else:
            line_fast='>'+str(contig.id)+'\n'+str(contig.seq)+'\n'
            crom.write(line_fast)
  return fasta_plas, fasta_crom
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
          cv = float(fields[11])
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
            if id <= 95 and cv <= 95 :#(cv <= 60 or cvs <= 60)#-----------------------------------------
              #print('verificando contigs que no pertenecen a cromosomas')
              if nodes.get(node, None) == None:
                nodes[node] = fields[1:]
                retorno[node]=1
                subjects.append(subject)
              salid.write(str(node)+'\n')
        contpas=contact
  #print(nodes, subjects)
  return retorno
#----------------------seleccion de contigs por homologia-------------------#
def blastn(fasta_file, out_dir,BDp,BDc,tipo):
    nodos={}
    subjects = []
    subjects1 = []
    if tipo =='c':
      #blast_csv = out_dir+'No_chromosomeDB_blastn3_olddb_out.csv'
      blast_csv = out_dir+'No_chromosomeDB_blastn.csv'
      #blast_csv = out_dir+'No_chromosomeDB_blastn_out.csv'
      #print(os.path.getsize(blast_csv))
      dir_nodos_filt=blast_csv[:-4]+'_filt.csv'
      if (os.path.isfile(blast_csv)==False or int(os.path.getsize(blast_csv))==0) and os.path.isfile(dir_nodos_filt)==False:
        print('Base de datos usada para homologia de secuencias en cromosoma: ')
        print(BDc)
        os.system('blastn -query '+fasta_file+' -db '+BDc+' -task blastn -max_target_seqs 1 -num_threads 20 -dust no -outfmt "10 qseqid sseqid pident length qstart qend sstart send qlen slen qcovus qcovs" -out ' + blast_csv)
        #---- filtrar tabla y borrar sobrantes
        filt_borr(blast_csv,dir_nodos_filt)
      else:
        print('\talineamiento blast prebio encontrado para cromosomas')
        if (os.path.isfile(dir_nodos_filt)==False or int(os.path.getsize(dir_nodos_filt))==0):
          filt_borr(blast_csv,dir_nodos_filt)
        else:
          print('\tAlineamiento Blastn hecho y filtrado')
      dir_nodos=blast_csv[:-10]+'nodos.txt'
      if os.path.isfile(dir_nodos)==False:
        nodos=process_blast_csv(dir_nodos_filt, subjects,'c',dir_nodos)
    elif tipo =='p':
      blast_csv =out_dir+'PlasmidsDB_blastn.csv'#prueba con base de datos antigua para comparar eficiencia y FP
      dir_nodos_filt=blast_csv[:-4]+'_filt.csv'
      if (os.path.isfile(blast_csv)==False or os.path.getsize(blast_csv)==0) and os.path.isfile(dir_nodos_filt)==False:
        print('Base de datos usada para homologia de secuencias en plasmido: ')
        print(BDp)
        os.system('blastn -query '+fasta_file+' -db '+BDp+' -task blastn -max_target_seqs 1 -num_threads 20 -dust no -outfmt "10 qseqid sseqid pident length qstart qend sstart send qlen slen qcovus qcovs" -out ' + blast_csv)
        filt_borr(blast_csv,dir_nodos_filt)
      else:
        print('\talineamiento blast prebio encontrado para cromosomas')
        if (os.path.isfile(dir_nodos_filt)==False or int(os.path.getsize(dir_nodos_filt))==0):
          filt_borr(blast_csv,dir_nodos_filt)
        else:
          print('\tAlineamiento Blastn hecho y filtrado')
      dir_nodos=blast_csv[:-10]+'nodos.txt'
      if os.path.isfile(dir_nodos)==False:
        nodos=process_blast_csv(dir_nodos_filt, subjects1,'p',dir_nodos)
    return nodos
#--------------------------Lectura de dictados desde dir--------------------#
def leer_dict(direc):
  nodos={}
  if os.path.isfile(direc):
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
  cov_file=salida+'/Profundidad/nodos_cov.txt'
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
  plasflow_txt= out_dir+ '/Plasflow/Plasflow_nodos.txt'
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
#-----------------------Buscar nodos de spades desde plasmidspades-----------#
def blastplasmid(plasm_filtrado,spades_filtrado,out_dir):
  blast_csv=plasm_filtrado.replace('fasta','blastn.csv')
  if os.path.isfile(blast_csv)==False:
    print('alineamineto blast no encontrado ..')
    os.system("makeblastdb -in %s -dbtype 'nucl' -out %s"%(plasm_filtrado,out_dir+'/db_spades_out'))
    os.system('blastn -query %s -db %s/db_spades_out -max_target_seqs 1 -max_hsps 3 -task blastn -dust no -outfmt "10 qseqid sseqid pident qcovus length slen qlen qstart qend sstart send" -out %s'%(spades_filtrado,out_dir, blast_csv))
  result_nodesbp = {}
  contact=''
  contpas=''
  nodos_file=out_dir+'/Plasmidspades_nodos.txt'
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
    print('\tCarpeta '+str(r)+'creada con exito ...\n')
  except Exception:
    traceback.print_exc()
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
  plasm_filtrado_blast={}
  bandera_plasmidsp=0
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
      log(ruta,'\tEnsamble previo encontrado..')
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
    print('\n\tFin proceso Spades')
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
      log(ruta,'\tHerramienta Quast para analisis de Ensamble ejecutado con exito')
      if os.path.isfile(rutal+'/Assembly/Spades/QC/report.txt'):
        proceso=3
    print('\n\tFin proceso Metrics')
    log(ruta,'\n\tFin proceso Metrics')
    #Warning plasmidspades
  if proceso==3:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/warnings.log') and os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.fasta')==False:
      print('\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, Ha culminado sin encontrar secuencias ..')
      log(ruta,'\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, Ha culminado sin encontrar secuencias ..')
      plasm_filtrado_blast=False#-----------
      proceso=7#111111111111111111111111111111111111Salto por warning plasmidspades
    else:
      proceso=4
    print('\n\tFin proceso warnings')
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
      print('\n\tInicio de ensamble de novo de plasmidos herramienta PlasmidSpades...')
      log(ruta,'\n\tInicio de ensamble de novo de plasmidos herramienta PlasmidSpades:\n\t'+str(datetime.now()))
      plasmispades_nodes=plasmidspades(fast1,fast2,fast1u,fast2u,kmer,rutal+'/Assembly/Plasmitest/Plasmidspades')
      #Warning plasmidspades
      if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/warnings.log')and os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/contigs.fasta')==False:
        print('\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, Ha culminado sin encontrar secuencias ..')
        log(ruta,'\n\tEnsamble de plasmidos previo herramienta PlasmidSpades, Ha culminado sin encontrar secuencias ..')
        plasm_filtrado_blast=False#-----------
        proceso=7
        print('\n\tFin proceso warnings')
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
    print('\n\tFin proceso profundidad_Plasmidspades')
    log(ruta,'\n\tFin proceso profundidad_Plasmidspades')
  if proceso ==6:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasmidspades/Plasmidspades_nodos.txt') and os.path.getsize(rutal+'/Assembly/Plasmitest/Plasmidspades/Plasmidspades_nodos.txt')>0 and repetir == 0:
      print('\n\tTabla de nodos Spades identificados por plasmidspades encontrado')
      log(ruta,'\n\tTabla de nodos Spades identificados por plasmidspades encontrado..')
      plasm_filtrado_blast_dir=rutal+'/Assembly/Plasmitest/Plasmidspades/Plasmidspades_nodos.txt'
      proceso=7
    else:
      print('\n\tIniciando busqueda de contigs ensamblados por PlasmidSpades en nodos Spades:\n\t'+str(datetime.now()))
      log(ruta,'\n\tIniciando busqueda de contigs ensamblados por PlasmidSpades en nodos Spades')
      plasm_filtrado_blast=blastplasmid(plasm_filtrado,spades_filtrado,rutal+'/Assembly/Plasmitest/Plasmidspades')
      print(plasm_filtrado_blast)
      proceso=7
    print('\n\tFin proceso nodos_plasmidspades')
    log(ruta,'\n\tFin proceso nodos_plasmidspades')
    log(ruta,'Ha ocurrido un error con el filtrado del archivo de plasmidSpades')
  if proceso==7:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Plasflow/Plasflow_nodos.txt')and os.path.getsize(rutal+'/Assembly/Plasmitest/Plasflow/Plasflow_nodos.txt')>0:
      print('\n\tTabla de contigs identificados por Plasflow encontrado')
      log(ruta,'\n\tTabla de nodos identificados por Plasflow encontrado..')
      plasflow_nodes_dir=rutal+'/Assembly/Plasmitest/Plasflow/Plasflow_nodos.txt'
      proceso=8
    else:
      print('\n\tIniciando busqueda de contigs de Plasmidos con herramienta Plasflow')
      log(ruta,'\n\tIniciando busqueda de contigs de Plasmidos con herramienta Plasflow:\n\t'+str(datetime.now()))
      plasflow_nodes = plasflow(spades_filtrado,rutal+'/Assembly/Plasmitest/')
      proceso=8
    print('\n\tFin proceso Plasflow')
    log(ruta,'\n\tFin proceso Plasflow'+str(proceso))
  if proceso==8:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Profundidad/nodos_cov.txt')and os.path.getsize(rutal+'/Assembly/Plasmitest/Profundidad/nodos_cov.txt')>0:
      print('\n\tTabla de nodos identificados por profundidad, encontrado')
      log(ruta,'\n\tArchivo de contigs identificados por profundidad, encontrado..')
      cov_nodes_dir=rutal+'/Assembly/Plasmitest/Profundidad/nodos_cov.txt'#**
      proceso=9
    else:
      print('\n\tIniciando proceso de extraccion de contigs por profundidad')
      log(ruta,'\n\tIniciando proceso de extraccion de contigs por profundidad:\n\t'+str(datetime.now()))
      cov_nodes = deep_cov(spades_filtrado,rutal+'/Assembly/Plasmitest/')
      proceso=9
    print('\n\tFin proceso profundidad')
    log(ruta,'\n\tFin proceso profundidad'+str(proceso))
  #proceso=11
  if proceso==9:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/No_chromosomeDB_nodos.txt') and os.path.getsize(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/No_chromosomeDB_nodos.txt')>0:
      print('\n\tTabla de nodos de cromosoma identificados por Homologia, encontrado')
      log(ruta,'\n\tTabla de nodos de cromosoma identificados por Homologia, encontrado..')
      no_chromosome_nodes_dir=rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/No_chromosomeDB_nodos.txt'
      proceso=10
    else:
      if os.path.isdir(rutal+'/Assembly/Plasmitest/Homologia/'):
        print('\tCarpeta de homologia cromosoma Encontrada')
      else:
        crearruta(rutal+'/Assembly/Plasmitest/Homologia/')
      if os.path.isdir(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/'):
        print('\tCarpeta de homologia cromosoma Encontrada')
      else:
        crearruta(rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/')
      print('\n\tIniciando proceso de extraccion de contigs por Homologia de cormosoma:')
      log(ruta,'\n\tIniciando proceso de extraccion de contigs por Homologia de cormosoma:\n\t'+str(datetime.now()))
      no_chromosome_nodes={}
      no_chromosome_nodes=blastn(spades_filtrado,rutal+'/Assembly/Plasmitest/Homologia/Cromosoma/',BDp,BDc,'c')
      proceso=10
    print('\n\tFinalizado proceso Homologia Cromosoma')
    log(ruta,'\n\t Finalizado proceso Homologia Cromosoma'+str(proceso))
  if proceso==10:
    if os.path.isfile(rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/PlasmidsDB_nodos.txt')and os.path.getsize(rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/PlasmidsDB_nodos.txt')>0:
      print('\n\tTabla de nodos identificados por Homologia Plasmidos, encontrado')
      log(ruta,'\n\tTabla de nodos identificados porHomologia, encontrado..')
      plasmid_nodes_dir=rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/PlasmidsDB_nodos.txt'
      proceso=11
    else:
      if os.path.isdir(rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/'):
        print('\tCarpeta de homologia plasmidos Encontrada')
      else:
        crearruta(rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/')
      print('\n\tIniciando proceso de extraccion de contigs por Homologia')
      log(ruta,'\n\tIniciando proceso de extraccion de contigs por Homologia:\n\t'+str(datetime.now()))
      plasmid_nodes={}
      plasmid_nodes=blastn(spades_filtrado,rutal+'/Assembly/Plasmitest/Homologia/Plasmidos/',BDp,BDc,'p')
      proceso=11
    print('\n\tFin proceso Homologia Plasmidos/')
    log(ruta,'\n\tFin proceso Homologia Plasmidos/')
  if plasm_filtrado_blast_dir!='None' :
    print('\tLeyendo nodos PlasmidSpades')
    plasm_filtrado_blast=leer_dict(plasm_filtrado_blast_dir)
  #print('nodos palsmidspades')
  #print(plasm_filtrado_blast)
  if plasflow_nodes_dir!='None':
    print('\tLeyendo nodos Plasflow')
    plasflow_nodes=leer_dict(plasflow_nodes_dir)
  #print('nodos plasflow')
  #print(plasflow_nodes)
  if cov_nodes_dir!='None':
    print('\tLeyendo nodos Cobertura')
    cov_nodes=leer_dict(cov_nodes_dir)
  #print('nodos covertura')
  #print(cov_nodes)
  if no_chromosome_nodes_dir!='None':
    print('\tLeyendo nodos BD de Cromosoma')
    no_chromosome_nodes=leer_dict(no_chromosome_nodes_dir)
  #print('nodos cromo homologos')
  #print(no_chromosome_nodes)
  if plasmid_nodes_dir!='None':
    print('\tLeyendo nodos BD de plasmidos')
    plasmid_nodes=leer_dict(plasmid_nodes_dir)  
  #print('nodos plasmidos homologos')
  #print(plasmid_nodes)
  #,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,Unificacion de resultados
  ingreso={}
  etiq=['plasm_filtrado_blast','plasflow_nodes','cov_nodes','plasmid_nodes','no_chromosome_nodes']
  conteoe=0
  for i in (plasm_filtrado_blast,plasflow_nodes,cov_nodes,plasmid_nodes,no_chromosome_nodes):
    if len(i.keys())>0:
      #print(str(i))
      ingreso[etiq[conteoe]]=i
    conteoe+=1
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
      tablaP_dir, tablaC_dir=consolidar(spades_filtrado,rutal,ingreso,etiq)
      proceso=12
  else:
    print('\n\tHa ocurrido un error con la busqueda de homologos')
    log(ruta,'Ha ocurrido un error con la busqueda de homologos')
  
  if proceso==12:
    if os.path.isfile(rutal+'/Annotation/Prokka/plasmids'+lis+'_P.gbk') and os.path.isfile(rutal+'/Annotation/Prokka/cromosome'+lis+'_C.gbk') and os.path.isfile(rutal+'/Annotation/Prokka/genome'+lis+'_T.gbk') and os.path.getsize(rutal+'/Annotation/Prokka/plasmids'+lis+'_P.gbk')>0 and os.path.getsize(rutal+'/Annotation/Prokka/cromosome'+lis+'_C.gbk')>0 and os.path.getsize(rutal+'/Annotation/Prokka/genome'+lis+'_T.gbk')>0:
    #if os.path.isfile(rutal+'/Annotation/Prokka/'+lis+'_T.gbk') and os.path.getsize(rutal+'/Annotation/Prokka/'+lis+'_T.gbk')>0:
      print('\n\tAnotacion, encontrada')
      log(ruta,'\n\tAnotacion, encontrada..')
      proceso=13
    else:
      print('\n\tIniciando proceso de Anotacion')
      if os.path.isdir(rutal+'/Annotation/'):
        print('Carpeta: '+rutal+'/Annotation/'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Annotation/')
      if os.path.isdir(rutal+'/Annotation/Prokka/'):
        print('Carpeta: '+rutal+'/Annotation/Prokka/'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Annotation/Prokka/')
      log(ruta,'\n\tIniciando proceso de Anotacion:\n\t'+str(datetime.now()))
      if os.path.isdir(rutal+'/Annotation/Prokka/plasmids'):
        print('Carpeta: '+rutal+'/Annotation/Prokka/plasmids'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Annotation/Prokka/plasmids')
      if os.path.isfile(rutal+'/Annotation/Prokka/plasmids/'+lis+'_P.gbk') ==False:
        prokka(tablaP_dir, rutal,lis,lis+'_P')
      if os.path.isdir(rutal+'/Annotation/Prokka/cromosome'):
        print('Carpeta: '+rutal+'/Annotation/Prokka/cromosome'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Annotation/Prokka/cromosome')
      if os.path.isfile(rutal+'/Annotation/Prokka/cromosome/'+lis+'_C.gbk') ==False:
        prokka(tablaC_dir,rutal,lis,lis+'_C')
      if os.path.isdir(rutal+'/Annotation/Prokka/genome'):
        print('Carpeta: '+rutal+'/Annotation/Prokka/genome'+'  Encontrado!!!!')
      else:
        crearruta(rutal+'/Annotation/Prokka/genome')
      if os.path.isfile(rutal+'/Annotation/Prokka/genome/'+lis+'_T.gbk')==False:
        prokka(spades_filtrado,rutal,lis,lis+'_T')
      proceso=13
  else:
    print('\n\tHa ocurrido un error con la generacion de consolidados')
    log(ruta,'Ha ocurrido un error con la generacion de consolidados')
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
  print('tamano lista de muestras = '+str(len(lista)))
  if str(type(lista))!='<class \'str\'>' and  str(type(lista))!='<type \'str\'>':#mientras abacas626,668 <class 'str'>
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
def linea_comando(r,c,k,d,p1,p2,s1,s2,n,repetir):
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
  pipeline(ruta,lista,d,kmers,BDc,BDp,p1,p2,s1,s2,repetir)
  if n != 'n':
    enviarnot(n)

def main(a,raiz):# funcion si encuentra argumentos de usuario -t, -b ...
  r,k,c,d,n=0,0,0,0,0
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
  kmers=''
  repetir=0
  try:
    opts, args = getopt.getopt(a, "er:c:k:d:n:1:2:3:4:")
    #print(opts)
  except getopt.GetoptError:
    print('\n\tError de argumentos de compilacion')
    sys.exit(2)
  for opt, arg in opts:
    if opt == "-r":
      ruta= os.path.realpath(arg)
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
    elif opt == "-e": 
      repetir=1#.............................
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
    elif opt=='-h':
      ss2=str(arg)
      b4=1
    else:
      print(opt)
      print('\tError en el ingreso de parametros ..')
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
    print(ruta,bacteria,kmers,indice_lec,p1,p2,s1,s2,correo)
    linea_comando(ruta,bacteria,kmers,indice_lec,p1,p2,s1,s2,correo,repetir)
    sys.exit()      
if __name__ == "__main__":
  print('Parametros :')
  print(sys.argv)
  raiz=str(sys.argv[0])[:-22]
  if len(raiz)==0:
      raiz='.'
  main(sys.argv[1:],raiz)

