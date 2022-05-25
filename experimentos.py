
#actualizacion 25/01/22
#--------------------------------------calculo de profundidad de nucleotido desde archivo generado por secuenciacion
def explorprofseq(ruta,lista):#explorar y hacer el histograma
  print(lista)
  #muestras='Porcentaje\t'
  muestras=''
  print(type(lista))
  #--------------buscar tamano promedio-
  lreads=ruta+'/reads1.txt'
  reads_len=300
  if os.path.isfile(lreads):
    print('\n\tArchivo de distribucion de tamanos de read encontrado')
    cantidad,tamanos,cantidad_n=[],[],[]
    suma=0
    with open(lreads, 'r') as lreads:
      for lr in lreads:
        lr=lr.split()
        cantidad.append(int(lr[0]))
        tamanos.append(int(lr[1]))
        suma+=int(lr[1])
    
    for can in cantidad:
      cantidad_n.append(can/suma)
    reads_len=0
    for i in range(0,len(cantidad)):
      reads_len+=cantidad_n[i]*tamanos[i]
    print('\t\tTamano calculado de la distribucion de reads '+str(reads_len))
  #-------------------
  if str(type(lista))[-3:-2]=='t':
    tablap=[]
    muestras='Muestra\tProfundidad_Ponderada\t'
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
            spade[str(fields[1])] = float(((float(fields[5])*reads_len)/(reads_len-127+1))*float(fields[3]))
            suma_tam+=int(fields[3])
      #print('dict'+str(spade))
      #print('suma genoma'+str(suma_tam))
      for contig in spade:
        prof_pon+=spade[contig]/suma_tam
      tablap.append(prof_pon)
  else:
    spade = {} 
    suma_tam=0
    prof_pon=0
    with open(spades_file, 'r') as input_fasta:
      for contig in SeqIO.parse(input_fasta, 'fasta'):
          line = contig.id.strip()
          fields = line.split('_')
          spade[str(fields[1])] = float(((float(fields[5])*reads_len)/(reads_len-127+1))*float(fields[3]))
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
  conteo=0
  for fila in tablap:
    linea=str(lista[conteo]+'\t'+str(fila))
    conteo+=1
    tab.write(linea+'\n')
  tab.close()
#----------------Recall variacion de distancia en las puntas----------------#
def recallplasmid_puntas(ruta,lista):
  print(lista)
  muestras='Distancia\t'
  if str(type(lista))[-3:-2]=='t':
    tablap=np.zeros((100,len(lista)+1))
    tablar=np.zeros((100,len(lista)+1))
    tablaf=np.zeros((100,len(lista)+1))
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
  fila=0
 
  for i in range(0,1000,10):
    tablap[fila][0]=i#indicar combinacion de id y cov 
    tablar[fila][0]=i
    tablaf[fila][0]=i
    neg_pr=np.zeros(len(lista))
    verdad=0
    conteo_m=0
    conteo_n=0
    #sys.exit()
    if str(type(lista))[-3:-2]=='t':
      vp=np.zeros(len(lista))
      tp=np.zeros(len(lista))
      for m in lista:
        blast_csv=ruta+'/'+m+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.blastn3_pos_csv'
        nodos={}
        contact=''
        contpas=''
        with open(blast_csv, mode='r') as blast_file:
          for line in blast_file:
            line = line.strip()
            fields = line.split(',')
            contact=fields[0]
            if contact!=contpas:
              qtam=float(fields[6])
              qal_ini=float(fields[7])
              qal_fin=float(fields[8])
              cvs= (float(fields[4])/float(fields[5]))*100
              cv = float(fields[3])
              ide = float(fields[2])
              if ide >= 90 and (cv >= 30 or cvs >= 30):
                if qal_ini < i or (qtam-qal_fin) < i :
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
        try:
          pre[conteo_m]=(vp[conteo_m]/me[conteo_m])
          #print(rec[ii])
          rec[conteo_m]=(vp[conteo_m]/(tp[conteo_m]))
          #print(pre[ii])
          fsc[conteo_m]=2*(((pre[conteo_m])*(rec[conteo_m]))/((pre[conteo_m])+(rec[conteo_m])))
          #print(fsc[ii])
        except:
          print('indeterminated (0>.<0)')
        fsc[np.isnan(fsc)]=0#eliminar los nan
        tablap[fila][conteo_m+1]=round(pre[conteo_m],4)#precision
        tablar[fila][conteo_m+1]=round(rec[conteo_m],4)#recall
        tablaf[fila][conteo_m+1]=round(fsc[conteo_m],4)#fscore
        conteo_m+=1
    else:
      vp=0
      tp=0
      blast_csv=ruta+'/'+lista+'/Assembly/Plasmitest/Plasmidspades/contigs.filter1000.blastn3_pos_csv'
      nodos={}
      contact=''
      contpas=''
      with open(blast_csv, mode='r') as blast_file:
        for line in blast_file:
          line = line.strip()
          fields = line.split(',')
          contact=str(fields[0])
          if contact != contpas:
            qtam=float(fields[6])
            qal_ini=float(fields[7])
            qal_fin=float(fields[8])
            ide = float(fields[2])
            cv = float(fields[3])
            cvs= (float(fields[4])/float(fields[5]))*100
            if ide >= 90 and (cv >= 30 or cvs >= 30):
              if qal_ini < i or (qtam-qal_fin) < i :
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
      #rec=(vp/me)
      #print(pre)
      #pre=(vp/(tp))
      #print(rec)
      print('nodos que pasan el filtro'+str(me))
      #print(vp)
      try:
        rec=(vp/me)
        print(pre)
        pre=(vp/(tp))
        print(rec)
        fsc=2*(((vp/me)*(vp/(tp)))/((vp/me)+(vp/(tp))))
      except:
        print('operacion indeterminada x/0')
      tablap[fila][1]=round(pre,4)#precision
      tablar[fila][1]=round(rec,4)#recall
      tablaf[fila][1]=round(fsc,4)#fscore
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
  tab=open(ruta+'/test/plasmidspades/tabla_fscore_bordes.txt','a')
  tab.write(muestras+'\n')
  for fila in tablaf:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
#----------------------------ABACAS
  """if proceso==9:
    print('Iniciando Abacas')#experimento
    if os.path.isfile(ruta+'/test/True/'+lis+'/Abacas/nodos.txt')==False:
      abacas(spades_filtrado,ruta+'/test/True/'+lis,lis)
    else:
      print('Archivo de etiquetado blast previo encontrado ..')
  else:
    log(ruta,'Ha ocurrido un error con la seleccion de contigs desde plasmidSpades a Spades')"""
  """if proceso==6:
    print('Iniciando extraccion de contigs verdaderos de variacion de selecion de contigs por homologia')#experimentotrulastn3_out_spe.csv'
    if os.path.isfile(ruta+'/test/True/'+lis+'/true_nodes_blastn3_out_spe.csv')==False:
      t_nodes=true_nodes(spades_filtrado,ruta+'/test/True/'+lis,lis)
    else:
      print('Archivo de etiquetado blast previo encontrado ..')
  else:
    log(ruta,'Ha ocurrido un error con la seleccion de contigs desde plasmidSpades a Spades')"""
#------------------------recall coverage------------------------------------#Experimento
def explorprof(ruta,lista):#explorar y hacer el histograma
  print(lista)
  #muestras='Porcentaje\t'
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
      muestras+=str(m)+'_ID\tKmer_Cov\tTipo_sec\t'
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
            spade[str(fields[1])] = float(fields[5])
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
    tab=open(ruta+'/test/profundidad/histograma.txt','a')
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
#-----------tabla de nodos sobre aliniamiento blastn de TRue nodes-----------#tambien para abacas
def tablanodos(ruta,lista):
  tabla=open('/vault2/homehpc/dtalero/Kpneumoniae/Reads/test/True/tabla_nodos_abacas.txt','a')#tabla_nodos3.txt'
  for m in lista:
    nodos=[]
    nodos_ord=[]
    agregar=m+'\t'
    #blast='/vault2/homehpc/dtalero/Kpneumoniae/Reads/test/True/'+m+'/true_nodes_blastn3_out_spe.csv'
    abacas='/vault2/homehpc/dtalero/Kpneumoniae/Reads/test/True/'+m+'/Abacas/'+m+'.MULTIFASTA.fa'
    contact=''
    contpas=''
    result_nodes={}
    """with open(blast, mode='r') as blast_file:
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
                  result_nodes[fields[0]] = 1
                  #print(fields[0])
            contpas=contact"""#funcion para tomar los nodos desde el aliniamiento de 3 MHPS
    with open(abacas, mode='r') as aba:
        for line in aba:
          if line[:1]=='>':
            datos=line.split('_')
            nodos.append(int(datos[1]))
    nodos_ord=np.sort(nodos)
    for nod in nodos_ord:
      agregar+=str(nod)+'\t'
    tabla.write(agregar+'\n')
  tabla.close()
#-------------- tabla de recall para todas las muestras---------------------#Experimento
def recallplasmid(ruta,lista):
  print(lista)
  if str(type(lista))[-3:-2]=='t':
    tablap=np.zeros((400,len(lista)+2))
    tablar=np.zeros((400,len(lista)+2))
    tablaf=np.zeros((400,len(lista)+2))
    vp=np.zeros(len(lista))
    me=np.zeros(len(lista))
    sc=np.zeros(len(lista))
    pc=np.zeros(len(lista))
    tp=np.zeros(len(lista))
    pre=np.zeros(len(lista))
    rec=np.zeros(len(lista))
    fsc=np.zeros(len(lista))
  else:
    tablap=np.zeros((400,3))
    tablar=np.zeros((400,3))
    tablaf=np.zeros((400,3))
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
  for i in range(90,100):
    jtam=0 
    for j in range(20,100,2):
      tablap[fila][0]=i
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
          spades=ruta+'/'+m+'/Ensamble/Spades/contigs.filter1000.fasta'
          renglones = os.popen("grep -c '>' "+spades).read()
          sc[conteo_m] = int(renglones.strip())#cabtidad de nodos spades
          plspades=ruta+'/'+m+'/Ensamble/Plasmidspades/contigs.filter1000.fasta'
          renglones = os.popen("grep -c '>' "+plspades).read()
          pc[conteo_m] = int(renglones.strip())#cabtidad de nodos Plasmidspades
          blast_csv=ruta+'/'+m+'/Ensamble/Plasmidspades/contigs.filter1000.blastn3_csv'
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
                    nodos[fields[0]]=1#son las etiquetas de los contigs de spades extraidas del alineamiento entre plasmid y spades
                  else:
                    print('repetido')
              contpas=contact           
          me[conteo_m]=conteo_n#Cantidad de nodos predichos por plasmid que pasan el filtro
          conteo_n=0
          #print('tamanio de nodos que pasaron el filtrado '+str(len(nodos)))
          true=ruta+'/test/True/'+m+'/verdaderos.txt'
          with open(true, mode='r') as verdaderos:
            for line in verdaderos:
              line = line.strip()
              fields = line.split(',')
              tp[conteo_m]+=1
              if nodos.get(str(fields[0]),None)!=None:
                vp[conteo_m]+=1          
          pre[conteo_m]=(vp[conteo_m]/me[conteo_m])
          #print(pre[ii])
          rec[conteo_m]=(vp[conteo_m]/(tp[conteo_m]))
          #print(rec[ii])
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
        spades=ruta+'/'+lista+'/Ensamble/Spades/contigs.filter1000.fasta'
        renglones = os.popen("grep -c '>' "+spades).read()
        sc = renglones.strip()#cabtidad de nodos spades
        plspades=ruta+'/'+lista+'/Ensamble/Plasmidspades/contigs.filter1000.fasta'
        renglones = os.popen("grep -c '>' "+plspades).read()
        pc = int(renglones.strip())#cabtidad de nodos Plasmidspades
        #print('nodos de plasmidspades'+str(pc))
        blast_csv=ruta+'/'+lista+'/Ensamble/Plasmidspades/contigs.filter1000.blastn3_csv'
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
                  nodos[fields[0]]=1#son las etiquetas de los contigs de spades extraidas del alineamiento entre plasmid y spades
                else:
                  print('repetido')
            contpas=contact
        me=conteo_n#Cantidad de nodos predichos por plasmid que pasan el filtro del alineamiento
        conteo_n=0
        
        true=ruta+'/test/True/'+lista+'/verdaderos.txt'#funcion de calificacion desde hits > 70 > 70
        with open(true, mode='r') as verdaderos:
          for line in verdaderos:
            line = line.strip()
            fields = line.split(',')
            #Los verdaderos contigs son seleccionados de un hit casi perfecto
            tp+=1
            if nodos.get(str(fields[0]),None)!=None:
              vp+=1
        pre=(vp/me)
        #print(pre)
        rec=(vp/(tp))
        #print(rec)
        print('nodos que pasan el filtro'+str(me))
        print(vp)
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
  tab=open(ruta+'/test/plasmidspades/tabla_precision.txt','a')
  for fila in tablap:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
  tab=open(ruta+'/test/plasmidspades/tabla_recall.txt','a')
  for fila in tablar:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
  tab=open(ruta+'/test/plasmidspades/tabla_fscore.txt','a')
  for fila in tablaf:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
  #por cada ID,COV sacar los contigs de cada muestra
  #comparar con los True
  #evaluar el conjunto en las metricas"""
#---------------------------------------------------------------------------------------
def experimento3(ruta,lista,indice):#funcion que obtiene la presicion de todo el set
#debe correrse desde interactivo sobre todas las muestras
  tabla_recall=[]
  dire=ruta.split('/')
  direbd=''
  for i in dire[:-1]:
    direbd+=i+'/'
  print(direbd)
  salida='/vault2/homehpc/dtalero/Kpneumoniae/Reads/test/plasmidspades/nodoplasmid_spads/'
  for i in range(90,100):
    for j in range(50,100):
      if os.path.isfile(salida+'/seq/ic_'+str(i)+'-'+str(j)+'.fasta')==False:
        fastasalida=open(salida+'/seq/ic_'+str(i)+'-'+str(j)+'.fasta','a')
        for lis in lista:
          blast_csv=ruta+'/'+lis+'/Ensamble/Plasmidspades/contigs.filter1000.blastn_csv'
          mfasta=ruta+'/'+lis+'/Ensamble/Spades/contigs.filter1000.fasta'
          nodos={}
          with open(blast_csv, mode='r') as blast_file:
            for line in blast_file:
              line = line.strip()
              fields = line.split(',')
              id = float(fields[2])
              cv = float(fields[3])
              if id > i and cv > j:
                if nodos.get(fields[1])==None:
                  nodos[fields[0]]=1
                else:
                  nodos[fields[0]]+=1
          for contig in SeqIO.parse(mfasta, 'fasta'):
            if nodos.get(str(contig.id)) != None:
              fastasalida.write('>'+str(contig.id)+'\n'+str(contig.seq)+'\n')
        fastasalida.close()
      else:
        print('archivo '+salida+'/seq/ic_'+str(i)+'-'+str(j)+'.fasta'+' ya existe')
  # blast
  secuencias=os.listdir(salida+'/seq')
  for sec in secuencias:
    if os.path.isfile(ruta+'/test/plasmidspades/nodoplasmid_spads/blast/'+sec):
      print('Archivo '+ruta+'/test/plasmidspades/nodoplasmid_spads/blast/'+sec+' encontrado ..')
    else:
      os.system('blastn -query %s -db %s -max_target_seqs 1 -max_hsps 1 -task blastn -dust no -num_threads 10 -outfmt "10 qseqid sseqid pident qcovus length" -out %s'%(salida+'seq/'+sec,direbd+'/BD_experimento/base_blast_exp', ruta+'/test/plasmidspades/nodoplasmid_spads/blast/'+sec))
  alineamientos=os.listdir(salida+'/blast')
  tabla=np.zeros((10,50))
  for alin in alineamientos:
    if alin[-5:]=='fasta':
      ypred=[]
      with open(salida+'/blast/'+alin, mode='r') as hits:
        for hit in hits:
          if (hit[int(hit.index('*'))+1]) == 'p':
            ypred.append(1)
          else:
            ypred.append(0)
      print(alin)
      pos_i=int(alin[int(alin.index('_'))+1]+alin[int(alin.index('_'))+2])-90
      pos_j=int(alin[int(alin.index('-'))+1]+alin[int(alin.index('-'))+2])-50
      yreal=np.ones(len(ypred))
      tabla[pos_i][pos_j]=(precision_score(ypred,yreal))     
  print(tabla)
  tab=open(salida+'/tabla/tabla_recall.txt','a')
  
  for fila in tabla:
    linea=''
    for col in fila:
      linea+=str(col)+'\t'
    tab.write(linea+'\n')
  tab.close()
#--------------Experimento contra base de datos local_plasmidspades__________#
def experimento2(plasm_filtrado,muestra,ruta):#funcion para verificar la presicion sobre los contigs de plasmidspades
  #correr blast.
  dire=ruta.split('/')
  direbd=''
  for i in dire[:-1]:
    direbd+=i+'/'
  print(direbd)
  blast_csv=plasm_filtrado.replace('fasta','blastn_csv')
  if os.path.isfile(ruta+'/test/'+muestra) == False:
    os.system('blastn -query %s -db %s -max_target_seqs 1 -max_hsps 1 -task blastn -dust no -outfmt "10 qseqid sseqid pident qcovus length" -out %s'%(plasm_filtrado,direbd+'/BD_experimento/base_blast_exp', ruta+'/test/'+muestra))  
  else:
    print('Archivo de alineamiento blast encontrado ')
  return True
#-------------Experimento de seleccion de contigs de plasmidspades-----------#
#funcion temporal para evaluar la seleccion de contigs de plasmidspades
def experimento1(plasm_filtrado,out_dir):# funcion que evalua la seleccion de contigs de blast a partir de conteo
  #saber cuantos contigs encuentra plasmidspades
  contigspl=0
  compon={}
  with open(plasm_filtrado) as multifasta:
    #print('el archivo tiene'+str(len(SeqIO.parse(multifasta, 'fasta'))))
    for contig in SeqIO.parse(multifasta, 'fasta'):
      contigspl+=1
      divi=str(contig.id).split('_')
      if compon.get(str(divi[7]))==None:
        compon[str(divi[7])]=int(divi[3])
      else:
        compon[str(divi[7])]+=int(divi[3])
      #print(divi[3],divi[7])
  print('muestra'+out_dir)
  print(compon)
  #print('cantidad de contigs = '+str(contigspl))
  blast_csv=plasm_filtrado.replace('fasta','blastn_csv')
  #result_nodes = {}
  tabla=[]
  iden=90
  cov=50
  cantidad=0
  linea=''
  abrir=open(out_dir+'/seleccionados.txt','w')
  abrir.close()
  nodos={}
  propor=[]
  with open(out_dir+'/seleccionados.txt','a') as salida:
  #salida.write('Cantidad de secuencias encontradas con plasmidspades = '+str(contigspl)+'\n')
    for i in range(90,100):
      for j in range(50,100):
        #print('rango identidad='+str(i)+' coverage='+str(j)+' cantidad='+str(cantidad))
        #cantidad=0
        nodos={}
        with open(blast_csv, mode='r') as blast_file:
          for line in blast_file:
            line = line.strip()
            fields = line.split(',')
            id = float(fields[2])
            cv = float(fields[3])
            if id > i and cv > j:
              plasm=str(fields[1]).split('_')
              spade=str(fields[1]).split('_')
              taman=str(fields[4])
              if nodos.get(str(plasm[7]))==None:
                nodos[str(plasm[7])]=int(taman)
              else:
                nodos[str(plasm[7])]+=int(taman)
        propor=[]
        for key in nodos:
          propor.append(round(nodos[key]/compon[key],5))
        if len(propor) > 0:
          linea+=str(round(statistics.mean(propor),4))+'\t'
        else:
          print(propor)
          linea+='0\t'
      linea+='\n'
      #print('escribiendo linea')
      salida.write(linea)
      linea=''
  print(nodos)
  return True