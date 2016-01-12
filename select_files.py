import glob
import sys
import numpy as np

#####Informacion extra
if len(sys.argv) == 1:

	print
	print 'Como ejecutar:', sys.argv[0]
	print 'python', sys.argv[0], '<keywords> <output> [<rango de archivos> <intervalo de archivos>]'
	print
	print 'Genera listas de archivos'
	print     
	print 'Los parametros entre [] son opcionales, pero deben ir ambos si es que se usan'
	print 'El rango de archivos debe tener el formato "archivo_inicial-archivo_final" ej. "4-10" (sin comillas)'
	print 'El intervalo de archivos define cada cuantas epocas va a guardar'
	print 'Output: "4-10 2": guarda los archivos encontrados numero 4, 6, 8 y 10'
	print

	sys.exit(1)

keywords = sys.argv[1]
output   = sys.argv[2]

if len(sys.argv) > 3:
    try:
        rango = sys.argv[3]
        salto = sys.argv[4]
    except IndexError:
        print 'Numero de parametros invalidos!'
        sys.exit()

    try:
        salto = int(salto)
        ep1, ep2 = rango.split('-')
        ep1, ep2 = int(ep1), int(ep2)
    except ValueError:
        print 'Parametros invalidos!'
        sys.exit()

archivos = glob.glob(keywords)
print 'Filtrando %d archivos encontrados...' % len(archivos)
archivos = archivos[ep1-1:ep2+1:salto]

print 'Guardando lista con %d archivos...' % len(archivos)
np.savetxt(output,archivos,fmt='%s')
