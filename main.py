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