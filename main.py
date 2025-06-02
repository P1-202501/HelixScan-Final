import logging
import os
from dotenv import load_dotenv
import sentry_sdk
import pandas as pd
# import matplotlib.pyplot as plt # Descomentar si se va a usar para visualizaciones
# import seaborn as sns # Descomentar si se va a usar para visualizaciones
from collections import Counter # Para el análisis de frecuencia de aminoácidos

# --- Carga de variables de entorno e inicialización de Sentry ---
load_dotenv() # Carga las variables de entorno desde el archivo .env

SENTRY_DSN = os.getenv("SENTRY_DSN")

if SENTRY_DSN:
    sentry_sdk.init(
        dsn=SENTRY_DSN,
        traces_sample_rate=1.0 # Tasa de muestreo para el monitoreo de rendimiento (ajustar según necesidad)
    )
    logging.info("Sentry inicializado correctamente.")
else:
    logging.warning("SENTRY_DSN no encontrado en .env. Sentry no estará activo.")


# --- Configuración del sistema de logging ---
# Crea un logger personalizado
logger_principal = logging.getLogger('sequence_processor')
logger_principal.setLevel(logging.DEBUG)

# Elimina cualquier handler existente para evitar que los mensajes vayan a la consola por defecto
for handler in logger_principal.handlers[:]:
    logger_principal.removeHandler(handler)
for handler in logging.root.handlers[:]: # También limpia el root logger
    logging.root.removeHandler(handler)

# Crea un FileHandler para escribir en el archivo
file_handler = logging.FileHandler('programa.log', encoding='utf-8')
file_handler.setLevel(logging.DEBUG) # Nivel mínimo para este handler

# Define el formato para el archivo
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(name)s - %(message)s')
file_handler.setFormatter(formatter)

# Añade el FileHandler al logger principal
logger_principal.addHandler(file_handler)

# --- Constantes para análisis de proteínas ---
MIN_LONGITUD_PROTEINA_FUNCIONAL = 5
MAX_LONGITUD_PROTEINA_FUNCIONAL = 100

# --- Tabla de codones ---
tabla_codones = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'STOP', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

































































































































































































































































































# --- Lógica principal del programa ---
def main():
    registrador = RegistradorSecuencias()
    gestor_adn = GestorADN(registrador)
    analizador_proteinas = AnalizadorProteinas(registrador)

    print("¡Bienvenido al Procesador de Secuencias Biológicas!")

    while True:
        print("\n--- Menú Principal ---")
        print("1. Procesar una nueva secuencia de ADN")
        print("2. Analizar datos de secuencias (desde 'datos_adn.csv')")
        print("3. Obtener cadena de ADN complementaria") # Nueva funcionalidad
        print("4. Analizar frecuencia de aminoácidos de una proteína") # Nueva funcionalidad
        print("5. Salir")
        
        opcion = input("Seleccione una opción (1, 2, 3, 4 o 5): ").strip()
        registrador.registrar_info(f"Usuario seleccionó la opción: {opcion}")

        if opcion == '1':
            try:
                adn_usuario = input("Ingrese la secuencia de ADN (ej. 5'-ATGC-3', o 'cancelar' para volver al menú): ").strip()
                if adn_usuario.lower() == 'cancelar':
                    print("Operación cancelada. Volviendo al menú principal.")
                    continue

                adn_limpio = gestor_adn.limpiar_adn(adn_usuario)
                
                if not gestor_adn.validar_adn(adn_limpio):
                    raise ValueError("La secuencia contiene caracteres inválidos. Solo se permiten A, T, C, G.")
                
                arn = gestor_adn.transcribir_adn(adn_limpio)
                proteina = analizador_proteinas.traducir_arn(arn)
                longitud_proteina = analizador_proteinas.analizar_longitud_proteina(proteina)

                analizador_proteinas.save_results_to_csv(adn_usuario, adn_limpio, arn, proteina)
                
                registrador.registrar_info("Transcripción y traducción completadas correctamente")
                
                print("\n--- Resultados del Procesamiento ---")
                print(f"ADN ingresado: {adn_usuario}")
                print(f"ADN limpio: {adn_limpio}")
                print(f"ARN: {arn}")
                print(f"Proteína: {proteina}")
                print(f"Longitud de Proteína: {longitud_proteina} aminoácidos")
                print("------------------------------------\n")

            except ValueError as e:
                registrador.registrar_error(f"Error de validación: {e}")
                print(f"Error de entrada: {e}\nPor favor, intente de nuevo o ingrese 'cancelar'.")
            except Exception as e:
                registrador.registrar_critico(f"Error inesperado: {e}", info_exc=True)
                print(f"Ocurrió un error inesperado: {e}\nPor favor, contacte al soporte técnico.")
                if SENTRY_DSN: sentry_sdk.capture_exception()
        
        elif opcion == '2':
            analizador_proteinas.cargar_y_analizar_datos()
        
        elif opcion == '3': # Nueva opción: Obtener cadena complementaria
            adn_usuario = input("Ingrese la secuencia de ADN para obtener su complementaria: ").strip()
            adn_limpio = gestor_adn.limpiar_adn(adn_usuario)
            if gestor_adn.validar_adn(adn_limpio):
                cadena_complementaria = gestor_adn.obtener_cadena_complementaria(adn_limpio)
                if cadena_complementaria: # Si no hubo problemas en obtenerla
                    print(f"\n--- Cadena Complementaria ---")
                    print(f"ADN original: {adn_usuario}")
                    print(f"ADN complementario: {cadena_complementaria}")
                    print("------------------------------\n")
            else:
                print("Error: La secuencia de ADN ingresada no es válida para obtener su complementaria. Solo se permiten A, T, C, G.")
                registrador.registrar_warning(f"Intento de obtener cadena complementaria con ADN inválido: '{adn_usuario}'")
        
        elif opcion == '4': # Nueva opción: Analizar frecuencia de aminoácidos
            proteina_usuario = input("Ingrese la secuencia de proteína para analizar su frecuencia (ej. FLSYWSTOP): ").strip().upper()
            analizador_proteinas.analizar_frecuencia_aminoacidos(proteina_usuario)

        elif opcion == '5':
            print("Gracias por usar el Procesador de Secuencias Biológicas. ¡Hasta pronto!")
            break
        
        else:
            print("Opción no válida. Por favor, seleccione una de las opciones del menú.")
            registrador.registrar_advertencia(f"Usuario ingresó una opción de menú inválida: {opcion}")

if __name__ == "main": 
    main()