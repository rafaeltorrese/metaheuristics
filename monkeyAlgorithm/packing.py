import sys, shutil, os
import numpy as np
from datetime import date
from timeit import default_timer as timer
from gurobipy import *

def dMinFrontera(L, W, i):
    '''
    Calcula la distancia minima entre los puntos de la malla y la frontera
    :param L:int. Longitud horizontal del contenedor
    :param W:int. Longitud vertical del contenedor
    :i:dict[i]. Punto i seleccionado
:return:  
    '''
    x,y = i # unpack
    dminx = min(L-x,x ) # |----------(x)---| Calcula la menor distancia sobre eje x a una coordenada x
    dminy = min(W-y,y)  # |----------(y)---| Calcula la menor distancia sobre eje y a una coordenada y
    return min(dminx,dminy) # La menor distancia entre eje x e y



def Nik(i:dict,j:dict,Rk:dict,Rl:dict,figura:str)->bool:
    '''
    :i:dict[i]
    :j:dict[j]
    :R:dict.Lista de radios
    :figura:int. Tipo de figura {1,2,3,4,5,6,7}
    Definición de norma para calcular el conjunto Nik
    Para ver definición Nik ver artículo Litvinchev, Ozuna
    Sin telescopía
    '''
    if i != j: # If i different from j
        if figura == 'circulo':
            return normaCirculo(i,j) < Rk + Rl # Mejorar esta función
        elif figura == 'rombo':
            return normaRombo(i,j) < R[k] + R[l]
        elif figura == 'cuadrado':
            return normaCuadrado(i,j) < R[k] + R[l]
        elif figura == 'octagono':
            return normaOctagono(i,j) < R[k] + R[l]
        elif figura == 'elipse':
            return normaElipse(i,j) < R[k] + R[l]
        elif figura == 'hexagono':
            return normaHexagono(i,j) < R[k] + R[l]
        elif figura == 'rectangulo':
            return normaRectangulo(i,j) < R[k] + R[l]
        else:
            print("El tipo de objeto no existe")
#
#
def Oik(i,j,k,l,figura): # Para telescopia
    '''
    Esta función genera restricciones de traslape considerando telescopía 
    :i:dict[i]
    :j:dict[j]
    :rad[k]:Lista de radios
    :rad[l]:Lista de radios, esta es un alias de la anterior
    :figura:int. Tipo de figura {1,2,3,4,5,6,7}
    Definición de norma para calcular el conjunto Nik
    Para ver definición Nik ver artículo Litvinchev, Ozuna
    Sin telescopía
    '''
    if i != j: # If i different from j
        if figura == 1:
            return abs(rad[k] - rad[l]) < normaCirculo(i,j) < (rad[k] + rad[l])
        elif figura == 2:
            return abs(rad[k] - rad[l]) < normaRombo(i,j) < (rad[k] + rad[l])
        elif figura == 3:
            return abs(rad[k] - rad[l]) < normaCuadrado(i,j) < (rad[k] + rad[l])
        elif figura == 4:
            return abs(rad[k] - rad[l]) < normaOctagono(i,j) < (rad[k] + rad[l])
        elif figura == 5:
            return abs(rad[k] - rad[l]) < normaElipse(i,j) < (rad[k] + rad[l])
        elif figura == 6:
            return abs(rad[k] - rad[l]) < normaHexagono(i,j) < (rad[k] + rad[l])
        elif figura == 7:
            return abs(rad[k] - rad[l]) < normaRectangulo(i,j) < (rad[k] + rad[l])
        else:
            print("El tipo de objeto no existe")


def makeAndChangeDirectory(instName,anidamiento,fecha):
    '''Crea el directorio para guardar los resultados y se cambia al directorio creado
    :instName:str. Nombre de la instancia para generar directorio con mismo nombre
:fecha:str. Fecha de creación de directorio'''
    dirName , ext = os.path.splitext(instName)
    if anidamiento == 1:
        dirName = dirName.lower()+'_'+'anidamiento'+'_'+fecha
    else:
        dirName = dirName.lower()+'_'+fecha
    #pathname = os.path.dirname(dirName) # Expresa nombre directorio de acuerdo a OS
    #
    #
    if os.path.isdir(dirName): # Si Directorio Existe
        print('El directorio [{}] ya existe'.format(dirName))
        ans = 'l'
        while ans != 'y' or ans != 'n': # Obligar a escribir 'y/n'
            ans = input(print('¿Deseas eliminar sobreescribir directorio?(y/n): ')).lower()
        if ans=='y' or ans=='yes':
            shutil.rmtree(dirName) # Eliminar directorio
            print('Se eliminó directorio [{}]'.format(dirName))
            os.mkdir(dirName)
            print('Se creó el directorio: [{}]'.format(dirName))
            os.chdir(dirName) # Accesar al directorio creado
        else:
            print('========================================')
            print('Abandonando programa')
            print('========================================')
            sys.exit("Sugerencia: Guardar archivos en otra parte")
    else: # Intenta crear directorio
        os.mkdir(dirName)
        print('Se creó el directorio: [{}]'.format(dirName))
        os.chdir(dirName) # Accesar al directorio creado


def impresionResultados(nombreDeLaInstancia,tipoDeFigura,teles,circsK,pointSolutions):
    if m.isMIP == 0:
        print('El modelo NO ES entero')
        exit(1)
    if m.status == GRB.Status.UNBOUNDED:
        print("El modelo no puede resolverse porque es infactible o no está acotado")
        exit(0)
    elif m.status == GRB.Status.INFEASIBLE:
        print("El modelo no puede resolverse porque es infactible")
        exit(0)
    if m.status == GRB.Status.OPTIMAL or GRB.status.INTERRUPTED:
        areaContenedor = L * W # Si se modulariza esta sección, cuidar estos parámetros
        funcObj = m.objVal
        gap = m.MIPGap * 100
        #areaCirculos = np.pi * funcObj
        #ratioOcupacion = areaCirculos / areaContenedor
        #areaDesocupada = areaContenedor - areaCirculos
        tEjecucion = m.RunTime # Guarda el valor del tiempo que duró el programa
        #
        #
        #
        if m.status == GRB.STATUS.OPTIMAL:
            print("Optimización terminó con solución óptima {}")
            print("Función objetivo {0:3.3f}".format(m.objVal))
            print("Gap del {0:5.2}%".format(m.MIPGap*100))
        elif m.status == GRB.STATUS.INTERRUPTED::
            print("El propceso fue interrumpido por usauario")
        #
        #
        #
        solucion = list(zip(circsK,pointSolutions))
        with open('solCoord'+'_'+nombreDeLaInstancia,'wt') as fw:
            print(nombreDeLaInstancia,file=fw)
            print(solucion,file=fw)
            print("Se guardó solución en {}".format('solCoord'+'_'+nombreDeLaInstancia))
        #
        #
        print("La duración fue de {0:5.3f} segundos".format(tEjecucion))
        #
        #
        nombreArchivoResultados = 'resultados'+'_'+nombreDeLaInstancia
        print("\n==========================")
        print("Escribiendo resultados en '{}'".format(nombreArchivoResultados))
        #
        with open(nombreArchivoResultados,'wt') as fout:
            print(nombreDeLaInstancia,file=fout)
            print('Con_telescopia:,{}'.format(teles), file=fout)
            print('Figura,{}'.format(tipoDeFigura) , file=fout)
            #print("Contenedor",file=fout)
            #print("Largo,{}".format(Largo),file=fout)
            #print("Ancho,{}".format(Ancho) , file=fout)
            #print("AreaContenedor,{}".format(areaContenedor) , file=fout)
            #
            print("Info_de_la_Malla",file=fout)
            #print("dx,{}".format(deltaX),file=fout)
            #print("dy,{}".format(deltaY),file=fout)
            #
            print("Info_de_los_objetos",  file=fout)
            print("NumObjetos,{}".format(len(rad)) , file=fout )
            #
            print("InfoDeLaSolucion", file=fout)
            print("Gap,{}".format(gap) , file=fout)
            print("TiempoEjecucion,{}".format(tEjecucion) , file=fout)
            print("FuncionObjetivo,{}".format(funcObj) , file=fout)


def fijarTiempoEjecucion():
    opcion = input('El tiempo deseas fijarlo en segundos, minutos u horas (s,m,h): ')
    if opcion == '' or opcion.lower()=='s':
        t = int(input("Cuántos segundos ejecutarás el modelo: "))
        return t
    elif opcion.lower()=='m':
        t = int(input("Cuántos minutos ejecutarás el modelo: "))
        return t*60 # Para ponerlo en minutos
    elif opcion.lower()== 'h':
        t = float(input("Cuántas horas ejecutarás el modelo: "))
        return t*3600 # Tiempo en horas



def printElapsedTimes(tvar,tassign,tdomain,toverlap,tnesting,tfunction,toptimizer,tub,tlb,filename):
    '''
    This function print times for constraints generation
    :param tvar:float Variable Time Generation in seconds
    :param tassign: float. Assignment Constraint Time 
    :param tdomain: float. Domain Constraint Time Generation in seconds
    :param toverlap: float. Overlap Constraint  Time Generation in seconds
    :param tnesting: float. Nesting Constraint Time Generation in seconds
    :param tfunction: float. Objective Function Time Generation in seconds
    :param toptimizer: float. Optimization Process Time Generation in seconds
    :param tub: float.  Upper Bound Constraint Time Generation in seconds
    :param tlb: float. Lower Bound Time Generation in seconds
    :return:
    A file with information about generation times
    '''
    with open('tiemposRestricciones_'+filename+'.txt','wt') as fout:
        print("variables,cotainf,cotasup,asignacion,dominio,traslape,anidacion,objetivo,optimizacion",file=fout)
        print("{:4.2f},{:4.2f},{:4.2f},{},{:10.2f},{},{},{:10.2f},{:10.2f}".format(tvar,tlb,tub,tassign,tdomain,toverlap,tnesting,tfunction,toptimizer),file=fout)



        
def modelPacking(radii,points,tel,figura):
    '''
    :radii: conjunto de radios k
    :points: Conjunto de puntos I, que conforman la malla
    :tel:{0:'sin telescopia',1:'con telescopía'}
    :figura:int. Indica tipo de objeto

    '''
    K = range(len(radii))
    I = range(len(points))
    m = Model("packing") # Cambio en la función objetivo
    startTimeVariablesCreation = timer()
    centros = {(k,i):m.addVar(vtype=GRB.BINARY,name="({},{})".format(k,i)) for k in K for i in I}
    endTimeVariablesCreation = timer()
    elapsedTimeVariableCreation = endTimeVariablesCreation - startTimeVariablesCreation
    print("Se crearon {} variables en {:9.2f} segundos".format(m.numvars,elapsedTimeVariableCreation))
    #
    #
    m.update
    #
    #
    ## RESTRICCIONES DEL MODELO ##

    # Constraint Lower Bound
    stime_lbConst = timer()
    for k in K:
        m.addConstr(minK[k] - quicksum(centros[k,i] for i in I) <= 0, name="r2a_{}".format(k) )
    etime_lbConst = timer()
    elapsedTime_lbConst = etime_lbConst - stime_lbConst # tiempo total de generación restricción cota inferior
    print("Restricción de cota inferior. \t GENERADA en {:10.2f} segundos".format(elapsedTime_lbConst))

    #
    # Constraint Upper Bound
    stime_ubConst = timer() # tiempo de inicio para generar restricción de cota superior
    for k in K:
        m.addConstr( quicksum(centros[k,i] for i in I) - maxK[k]  <= 0, name="r2b_{}".format(k) )
    etime_ubConst = timer() # tiempo de finalización para generar restricción de cota superior
    elapsedTime_ubConst = etime_ubConst - stime_ubConst # tiempo total de generación restricción cota inferior
    print("Restricción de cota superior. \t GENERADA en {:14.2f} segundos".format(elapsedTime_ubConst))
    #
    #
    #
    # Generar restricción de dominio
    stime_domainConst = timer() # tiempo de inicio para generar restricción de dominio
    for i in I:
        for k in K:
            m.addConstr(rad[k]*centros[k,i] <= dMinFrontera(L=Largo,W=Ancho,points[i]))
    etime_domainConst = timer() # tiempo de inicio para generar restricción de dominio
    elapsedTime_domainConst = etime_domainConst - stime_domainConst
    print("Restricción de dominio.\t GENERADA en {:10.2f} segundos".format(elapsedTime_domainConst))
    #
    #
    #
    #
    #
    if tel == 0:# Restricción de traslape sin telescopia
        # Asignación
        stime_assignment = timer() # tiempo de inicio para generar restricción de cota superior
        for i in I:
            m.addConstr(quicksum(centros[k,i] for k in K) <= 1, name="r3_{}".format(i))
        etime_assignment = timer()
        elapsedTime_assignment = etime_assignment - stime_assignment # tiempo total de generación restricción cota inferior
        print("Restricción de asignación. \t GENERADA en {:14.2f} segundos".format(elapsedTime_assignment))
        #
        #
        
        #print("Generando restricción sin telescopía")
        stime_overlapConst = timer()
        for i in points:
            #for j in nPuntos[1:]:
            for j in points:
                for l in K:
                    for k in K:
                        if Nik(i,j,k,l,tipoDeObjeto): # precaución en esta función
                            m.addConstr(centros[k,i] + centros[l,j] <= 1)
        etime_overlapConst = timer()
        elapsedTime_overlapConst = etime_overlapConst - stime_overlapConst
        elapsedTime_nestConst = 'N/a' # No hay restricción de anidamiento
        print("Restricción SIN telescopía generada en {:10.2f} segundos".format(elapsedTime_overlapConst))
    #
    #
    #
    if tel == 1:# Restricción de traslape CON telescopia
        print("Generando restricción con telescopía")
        stime_nestConst = timer()
        for i in nPuntos:
            for j in nPuntos:
                for l in K:
                    for k in K:
                        if Oik(i,j,k,l,tipoDeObjeto): # precaución en esta función
                            m.addConstr(centros[k,i] + centros[l,j] <= 1)
        etime_nestConst = timer()
        elapsedTime_nestConst = etime_nestConst - stime_nestConst
        elapsedTime_overlapConst = 'N/a' # no hubo restriccion de traslape
        elapsedTime_assignment = 'N/a' # no hubo restriccion de asignación
        print("Restricción con TELESCOPÍA generada en {:10.2f} segundos".format(elapsedTime_nestConst))
        #
        #
        #
    # FUNCIÓN OBJETIVO #
    #m.setObjective(quicksum(rad[k]*rad[k] * centros[k,i] for i in nPuntos for k in K), sense=GRB.MAXIMIZE)
    print("Generando Función Objetivo")
    stime_fo = timer()
    #m.setObjective(quicksum(centros[k,i] for i in nPuntos for k in K), sense=GRB.MAXIMIZE)
    m.setObjective(quicksum(rad[k]*rad[k] * centros[k,i] for i in nPuntos for k in K), sense=GRB.MAXIMIZE)
    etime_fo = timer()
    elapsedTime_fo =  etime_fo - stime_fo
    print("Función objetivo.\t GENERADA en {:10.2f} segundos".format(elapsedTime_fo))
    print("\n=========")
    print("Comenzando con la optimización")
    #
    #
    stime_optimizationProcess = timer()
    #
    #
    #
    m.optimize() # Esta instrucción indica el comienzo de la optimización
    #
    #
    #
    print("\n============")
    etime_optimizationProcess = timer()
    elapsedTime_optimizationProcess = etime_optimizationProcess - stime_optimizationProcess
    print("¡¡Optimización Terminada!! en {:10.2f} segundos".format(elapsedTime_optimizationProcess))    
    #
    #
    #
    #

    if telescopia==1:
        nombreArchivoSolucion = InstanceName+"_"+"anidamiento"+".sol"
        nombreArchivoParams = InstanceName+"_"+"anidamiento"+".prm"
        m.write(nombreArchivoSolucion)
        m.write(nombreArchivoParams)
    else:
        nombreArchivoSolucion = InstanceName+".sol"
        nombreArchivoParams = InstanceName+".prm"
        m.write(nombreArchivoSolucion)
        m.write(nombreArchivoParams)

    if os.path.isfile(nombreArchivoSolucion):
        mostrarDibujo(objetosTipoK,listaSolucion,L,W,nombreInstancia,telescopia,tipoDeObjeto)
    
    #
    #
    #
    printElapsedTimes(tvar=elapsedTimeVariableCreation,tassign=elapsedTime_assignment,tdomain=elapsedTime_domainConst,toverlap=elapsedTime_overlapConst,tnesting=elapsedTime_nestConst,tfunction=elapsedTime_fo,toptimizer=elapsedTime_optimizationProcess,tub=elapsedTime_ubConst,tlb=elapsedTime_lbConst,filename=InstanceName)
        
    os.chdir(os.pardir) # regresar al directorio superior
