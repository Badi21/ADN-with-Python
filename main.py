import regex as re

PROTEINAS = {}
ADN = {}

CODONES = {'AGG': 'R', 'AGA': 'R', 'AGC': 'S', 'AGU': 'S',
           'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'AAU': 'N',
           'ACG': 'T', 'ACA': 'T', 'ACC': 'T', 'ACU': 'T',
           'AUG': 'M', 'AUA': 'I', 'AUC': 'I', 'AUU': 'I',

           'CGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGU': 'R',
           'CAG': 'Q', 'CAA': 'Q', 'CAC': 'H', 'CAU': 'H',
           'CCG': 'P', 'CCA': 'P', 'CCC': 'P', 'CCU': 'P',
           'CUG': 'L', 'CUA': 'L', 'CUC': 'L', 'CUU': 'L',

           'UGG': 'W',    'UGA': 'STOP', 'UGC': 'C', 'UGU': 'C',
           'UAG': 'STOP', 'UAA': 'STOP', 'UAC': 'Y', 'UAU': 'Y',
           'UCG': 'S',    'UCA': 'S',    'UCC': 'S', 'UCU': 'S',
           'UUG': 'L',    'UUA': 'L',    'UUC': 'F', 'UUU': 'F',

           'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GGU': 'G',
           'GAG': 'E', 'GAA': 'E', 'GAC': 'D', 'GAU': 'D',
           'GCG': 'A', 'GCA': 'A', 'GCC': 'A', 'GCU': 'A',
           'GUG': 'V', 'GUA': 'V', 'GUC': 'V', 'GUU': 'V',
           }


def Interprete(cadena):
    patronTituloGen = r'>.*\n' #ER que admite todos los titulos empezados por >C.
    patronTituloGenVer = r'\>.+\s{5}\d+.*' #ER que admite solo los que esten bien formulados
    erTituloGen = re.compile(patronTituloGen)
    erTituloGenVer = re.compile(patronTituloGenVer)
    patronFASTA = r'(((A|T|G|C){10}\s){5})|(((A|T|G|C){10}\s){0,4}(A|T|G|C){1,10}\n)' #Patron FASTA del cuerpo del ADN
    erFASTA = re.compile(patronFASTA)
    NameGen = ''    # Valor del nombre del Gen
    Gen = ''        # Valor del Gen
    bien = True     # Boolean que indica si esta correcta la expresion

    for linea in cadena:
        res = erTituloGen.fullmatch(linea) #Comprueba si es un titulo de Gen
        if res:
            ver = erTituloGenVer.fullmatch(NameGen.strip()) #Comprueba que el titulo este correctamente escrito
            if not ver:
                bien = False

            if bien: #Si se encuentra un titulo y esta bien, añade al diccionario ADN, la informacion anterior almacenada
                ADNaux = {NameGen.strip():Gen}
                ADN.update(ADNaux)
            if not bien:
                print('ERROR > La cadena de ADN: ', NameGen.strip(),' no se ha podido validar como modelo FASTA')
            Gen = ''
            NameGen = linea
            bien = True
        else:
            res2 = erFASTA.fullmatch(linea)
            if linea != '\n':
                if res2:
                    Gen = str(Gen + linea)
                else:
                    bien = False

    if bien:  #Introduce el ultimo gen si esta bien
        ADNaux = {NameGen.strip(): Gen}
        ADN.update(ADNaux)

def QuitarBlancos(archivo): #Pone en una linea y sin espacios la cadena que le entre
    adn = ''
    for linea in archivo:
        linea = linea.strip()
        linea = linea.replace(' ', '')
        adn = adn + linea
    return adn

def Cambiar_T_por_U(linea): #Cambia las Ts por Us
    return linea.replace('T','U')

def Generar_Proteina(clave, cadenaADN, n, l):
    i = 0
    cadena_proteina = ''
    terminar = False
    while terminar == False:
        cadena = cadenaADN[int(i):int(i+3)]
        proteina = str(CODONES.get(cadena))
        if i == n:
            terminar = True
        else:
            cadena_proteina = str(cadena_proteina) + str(proteina)
        i = i + 3

    comprobarFinalStop = str(cadena_proteina[len(cadena_proteina)-4:len(cadena_proteina)])
    cont = int(cadena_proteina.count('STOP'))
    if comprobarFinalStop == 'STOP' and cont == 1:
        if str(cadena_proteina[0:1]) == 'M':
            cadenaprefijo = clave[0:l]
            cadenasufijo1 = clave[l + len(str(n)):l + len(str(n)) + 3]
            cadenasufijo1 = cadenasufijo1.replace('nt', 'aa')
            cadenasufijo2 = clave[l + len(str(n)) + 3:]
            clave = cadenaprefijo + str(int(n / 3) - 1) + cadenasufijo1 + cadenasufijo2
            cadena_proteina = cadena_proteina.replace('STOP','')
            cadena_proteina = ArreglarCadena(cadena_proteina)

            PROTEINASaux = {clave:cadena_proteina}
            PROTEINAS.update(PROTEINASaux)
        else:
            print('ERROR > La cadena de ADN: ', clave, 'no comienza con el codon AUG')
    else:
        print('ERROR > La cadena de ADN: ', clave, 'no tiene el caracter de STOP en su sitio')

def ComprobarNumGen (clave, num, cadena):

    numero = int(num%3)

    if numero == 0:
        if len(cadena) == num:
            return True
        else:
            print('ERROR > La cade de ADN: ', clave, 'no posee la misma longitud como la que indica')
            return False
    else:
        print('ERROR > La cade de ADN: ', clave, 'no posee un numero multiplo de 3')
        return False

def LeerNombreArchivo():
    verificado = False
    NombreArchivo = r'(.*)(DNA.txt)'
    erNombreArchivo = re.compile(NombreArchivo)
    archivo = ''
    while not verificado:
        archivo = input('Nombre del Archivo > ')
        if (erNombreArchivo.fullmatch(archivo)):
            try:
                open(archivo)
                verificado = True

            except:
                print('ERROR > El Archivo no ha sido encontrado')
        else:
            print('ERROR > Nombre de Archivo incorrecto, vuelva a introdir uno nuevo del formato "...DNA.txt"')

    nom = erNombreArchivo.fullmatch(archivo)
    nom.group(1)
    return nom.group(1)

def ArreglarCadena(cadena):
    i = 0
    k = 0
    nuevacadena = ''
    for c in cadena:
        if i == 10:
            nuevacadena = nuevacadena + ' '
            k = k + 1
            i = 0
        if k == 5:
            k = 0
            nuevacadena = nuevacadena + '\n'
        i = i + 1
        nuevacadena = nuevacadena + c
    return nuevacadena.strip()

if __name__ == '__main__':
    nomArchivo = LeerNombreArchivo()
    nomProteninas = nomArchivo + 'Protein.txt'
    nomArchivo = nomArchivo + 'DNA.txt'
    archivo = open(nomArchivo, 'r')
    Interprete(archivo)
    proteinas = open(nomProteninas, 'w')
    patronNumero = r'(.*\s+)(\d+)(.*)'
    erNumeroADN = re.compile(patronNumero)
    patronParte2 = r'(.*\s+)(\d+)(.{3})(.*)'
    erParte2 = re.compile(patronParte2)

    for clave in ADN:
        result = erNumeroADN.fullmatch(clave)
        if result:
            leng = str(result.group(1))
            leng = len(leng)
            num = int(result.group(2))
            valor = str(ADN.get(clave))
            valor1 = str(Cambiar_T_por_U(valor))
            valor2 = str(QuitarBlancos(valor1))
            if ComprobarNumGen(clave, num, valor2):
                str(Generar_Proteina(clave, valor2, num, leng))

    for prot in PROTEINAS:
        proteinas.write(prot)
        proteinas.write('\n')
        proteinas.write(PROTEINAS.get(prot))
        proteinas.write('\n')
        proteinas.write('\n')

    proteinas.close()

    salir = False
    while not salir:
        print('Introduce una expresion regular que buscar dentro del Archivo. Enter en blanco para salir')
        Introducido = input(' ER >> ')

        if Introducido == '':
            salir = True
        else:
            print('===============================================================')
            Introducido = str('\>'+ str(Introducido))
            erIntroducido = re.compile(Introducido)#C\.Acc.*'
            coincidencias = 0
            for p in PROTEINAS:
                res = erIntroducido.fullmatch(str(p.strip()))
                if res:
                    res2 = erParte2.fullmatch(p)
                    nombreGen = res2.group(1)
                    coincidencias = coincidencias + 1
                    print('Nombre:', nombreGen[1:len(nombreGen)])
                    print()
                    nt =  int((int(res2.group(2)) + 1) * 3)
                    print('ADN:', nt, res2.group(3).replace('aa', 'nt'), res2.group(4))
                    print()
                    adn = res2.group(1) + str(nt).strip() + res2.group(3).replace('aa', 'nt') + res2.group(4)
                    print(ADN.get(adn))
                    print('Proteína:', res2.group(2), res2.group(3), res2.group(4))
                    print()
                    print(PROTEINAS.get(p))
                    print()
                    print('===============================================================')
                    print()
            print(coincidencias,'coincidencias' )
            print()
            print('===============================================================')
            print()

