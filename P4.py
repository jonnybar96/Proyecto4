# 
# ## Universidad de Costa Rica
# ### Escuela de Ingeniería Eléctrica
# #### IE0405 - Modelos Probabilísticos de Señales y Sistemas
# 
# Primer semestre del 2021
# 
# ---
# 
# * Estudiante: **Jonathan Barrantes Castillo**
# * Carné: **B50891**
# * Grupo: **2**

# 4.1. - Modulación 16-QAM
# (50%) Realice una simulación del sistema de comunicaciones como en la sección 3.2., pero utilizando una modulación 
#16-QAM en lugar de una modulación BPSK. Deben mostrarse las imágenes enviadas y recuperadas y las formas de onda.

# Librerías
import numpy as np
import matplotlib.pyplot as plt
import time
# Función moduladora
def modulador(bits, fc, mpp):
    '''Un método que simula el esquema de 
    modulación digital BPSK.

    :param bits: Vector unidimensional de bits
    :param fc: Frecuencia de la portadora en Hz
    :param mpp: Cantidad de muestras por periodo de onda portadora
    :return: Un vector con la señal modulada
    :return: Un valor con la potencia promedio [W]
    :return: La onda portadora c(t)
    :return: La onda cuadrada moduladora (información)
    '''
    # 1. Parámetros de la 'señal' de información (bits)
    N = len(bits) # Cantidad de bits

    # 2. Construyendo un periodo de la señal portadora c(t)
    Tc = 1 / fc  # periodo [s]
    t_periodo = np.linspace(0, Tc, mpp) # MPP Muestras por período
    portadora1 = np.sin(2*np.pi*fc*t_periodo) # Portadora 1
    portadora2 = np.cos(2*np.pi*fc*t_periodo) # Portadora 2

    # 3. Inicializar la señal modulada s(t)
    t_simulacion = np.linspace(0, N*Tc, N*mpp) 
    senal_Tx = np.zeros(t_simulacion.shape)
    moduladora = np.zeros(t_simulacion.shape)  # señal de información

    # 4. Asignar las formas de onda según los bits (QAM)
    x_1 = 1
    x_2 = 1
    for i, bit in enumerate(bits):
          
        if (bit == 0) and (x_1%2 != 0) and (x_2%2 != 0):
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1*-1
            moduladora[i*mpp : (i+1)*mpp] = 0
            x_1 += 1
            x_2 += 1
            
        if (bit == 1) and (x_1%2 != 0) and (x_2%2 != 0):
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1*1
            moduladora[i*mpp : (i+1)*mpp] = 1
            x_1 += 1
            x_2 += 1
         
        if (bit == 1) and (x_1%2 == 0) and (x_2%2 != 0):
            senal_Tx[i*mpp : (i+1)*mpp] = portadora2*1
            moduladora[i*mpp : (i+1)*mpp] = 1
            x_1 += 1
            x_2 += 1
            
        if (bit == 0) and (x_1%2 == 0) and (x_2%2 != 0):
            senal_Tx[i*mpp : (i+1)*mpp] = portadora2*-1
            moduladora[i*mpp : (i+1)*mpp] = 0
            x_1 += 1
            x_2 += 1
            
        if (bit == 0) and (x_1%2 != 0) and (x_2%2 == 0):
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1*-3
            moduladora[i*mpp : (i+1)*mpp] = 0
            x_1 += 1
            x_2 += 1
            
        if (bit == 1) and (x_1%2 != 0) and (x_2%2 == 0):
            senal_Tx[i*mpp : (i+1)*mpp] = portadora1*3
            moduladora[i*mpp : (i+1)*mpp] = 1
            x_1 += 1
            x_2 += 1
         
        if (bit == 1) and (x_1%2 == 0) and (x_2%2 == 0):
            senal_Tx[i*mpp : (i+1)*mpp] = portadora2*3
            moduladora[i*mpp : (i+1)*mpp] = 1
            x_1 += 1
            x_2 += 1
            
        if (bit == 0) and (x_1%2 == 0) and (x_2%2 == 0):
            senal_Tx[i*mpp : (i+1)*mpp] = portadora2*-3
            moduladora[i*mpp : (i+1)*mpp] = 0
            x_1 += 1
            x_2 += 1
         
            

    # 5. Calcular la potencia promedio de la señal modulada
    Pm = (1 / (N*Tc)) * np.trapz(pow(senal_Tx, 2), t_simulacion)
    
    LaPortadora = portadora1 + portadora2
    
    return senal_Tx, Pm, LaPortadora, moduladora  
    
    
def demodulador(senal_Rx, LaPortadora, mpp):
    '''Un método que simula un bloque demodulador
    de señales, bajo un esquema BPSK. El criterio
    de demodulación se basa en decodificación por 
    detección de energía.

    :param senal_Rx: La señal recibida del canal
    :param portadora: La onda portadora c(t)
    :param mpp: Número de muestras por periodo
    :return: Los bits de la señal demodulada
    '''
    # Cantidad de muestras en senal_Rx
    M = len(senal_Rx)

    # Cantidad de bits en transmisión
    N = int(M / mpp)

    # Vector para bits obtenidos por la demodulación
    bits_Rx = np.zeros(N)

    # Vector para la señal demodulada
    senal_demodulada = np.zeros(M)

    # Energía de un período de la portadora
    Es = np.sum(LaPortadora**2)

    # Demodulación
    for i in range(N):
        # Producto interno de dos funciones
        producto = senal_Rx[i*mpp : (i+1)*mpp] * LaPortadora
        senal_demodulada[i*mpp : (i+1)*mpp] = producto
        Ep = np.sum(producto) 

        # Criterio de decisión por detección de energía
        if Ep > Es*0:
            bits_Rx[i] = 1
        else:
            bits_Rx[i] = 0

    return bits_Rx.astype(int), senal_demodulada



# Parámetros
fc = 5000  # frecuencia de la portadora
mpp = 20   # muestras por periodo de la portadora
SNR = -5   # relación señal-a-ruido del canal

# Iniciar medición del tiempo de simulación
inicio = time.time()

# 1. Importar y convertir la imagen a trasmitir
imagen_Tx = fuente_info('arenal.jpg')
dimensiones = imagen_Tx.shape

# 2. Codificar los pixeles de la imagen
bits_Tx = rgb_a_bit(imagen_Tx)

# 3. Modular la cadena de bits usando el esquema BPSK
senal_Tx, Pm, LaPortadora, moduladora = modulador(bits_Tx, fc, mpp)

# 4. Se transmite la señal modulada, por un canal ruidoso
senal_Rx = canal_ruidoso(senal_Tx, Pm, SNR)

# 5. Se desmodula la señal recibida del canal
bits_Rx, senal_demodulada = demodulador(senal_Rx, LaPortadora, mpp)

# 6. Se visualiza la imagen recibida 
imagen_Rx = bits_a_rgb(bits_Rx, dimensiones)
Fig = plt.figure(figsize=(10,6))

# Cálculo del tiempo de simulación
print('Duración de la simulación: ', time.time() - inicio)

# 7. Calcular número de errores
errores = sum(abs(bits_Tx - bits_Rx))
BER = errores/len(bits_Tx)
print('{} errores, para un BER de {:0.4f}.'.format(errores, BER))

# Mostrar imagen transmitida
ax = Fig.add_subplot(1, 2, 1)
imgplot = plt.imshow(imagen_Tx)
ax.set_title('Transmitido')

# Mostrar imagen recuperada
ax = Fig.add_subplot(1, 2, 2)
imgplot = plt.imshow(imagen_Rx)
ax.set_title('Recuperado')
Fig.tight_layout()

plt.imshow(imagen_Rx)

# 4.2. - Estacionaridad y ergodicidad
# (30%) Realice pruebas de estacionaridad y ergodicidad a la señal modulada senal_Tx y obtenga conclusiones sobre estas.

# Librerías
import numpy as np
import matplotlib.pyplot as plt

# Tiempo para la muestra
t_x = np.linspace(0, 0.1,100)

# Posibles valores de A
A = [1,-1]

# Formas de onda
X_t = np.empty((4, len(t_x)))	   # 4 funciones del tiempo x(t) 

# Nueva figura 
plt.figure()

# Matriz con los valores de cada función posibles   
for i in A:
    x1 = i * np.cos(2*(np.pi)*fc*t_x) +  i*np.sin(2*(np.pi)*fc*t_x)
    x2 = -i * np.cos(2*(np.pi)*fc*t_x) +  i*np.sin(2*(np.pi)*fc*t_x) 
    X_t[i,:] = x1
    X_t[i+1,:] = x2
    plt.plot(t_x,x1, lw=2)
    plt.plot(t_x, x2, lw=2)       

# Promedio de las 4 realizaciones en cada instante 
P = [np.mean(X_t[:,i]) for i in range(len(t_x))]
plt.plot(t_x, P, lw=6,color='k',label='Promedio Realizaciones')

# 7. Graficar el resultado teórico del valor esperado
E = np.mean(senal_Tx)*t_x  # Valor esperado de la señal 
plt.plot(t_x, E, '-.', lw=3,color='c',label='Valor teórico')

# 8. Mostrar las realizaciones, y su promedio calculado y teórico
plt.title('Realizaciones del proceso aleatorio $X(t)$')
plt.xlabel('$t$')
plt.ylabel('$x_i(t)$')
plt.legend()
plt.show()


# 4.3. - Densidad espectral de potencia
# (20%) Determine y grafique la densidad espectral de potencia para la señal modulada senal_Tx.

#Realización del proceso = señal_Tx (quitamos el E y dejamos unicamente el (absXtw)2)

from scipy import fft

# Transformada de Fourier
senal_f = fft(senal_Tx)

# Muestras de la señal
Nm = len(senal_Tx)

# Número de símbolos (198 x 89 x 8 x 3)
Ns = Nm // mpp

# Tiempo del símbolo = periodo de la onda portadora
Tc = 1 / fc

# Tiempo entre muestras (período de muestreo)
Tm = Tc / mpp

# Tiempo de la simulación
T = Ns * Tc

# Espacio de frecuencias
f = np.linspace(0.0, 1.0/(2.0*Tm), Nm//2)

# Gráfica
plt.plot(f, 2.0/Nm * np.power(np.abs(senal_f[0:Nm//2]), 2))
plt.xlim(0, 20000)
plt.grid()
