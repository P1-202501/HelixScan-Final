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
}# --- Clases ---

class RegistradorSecuencias:
    #Clase para encapsular las operaciones de logging para secuencias de ADN/ARN/Proteínas.
    
    def __init__(self, nombre_registrador='simulador_biologico'):
    
        
        self.registrador = logging.getLogger(nombre_registrador)

    def registrar_info(self, mensaje):
        """Registra un mensaje de información."""
        self.registrador.info(mensaje)

    def registrar_advertencia(self, mensaje, extra_data=None):
      
        self.registrador.warning(mensaje)
        if SENTRY_DSN:
            with sentry_sdk.push_scope() as scope:
                if extra_data:
                   
                    scope.set_extras(extra_data)
                sentry_sdk.capture_message(f"Advertencia/Anomalía: {mensaje}", level='warning')

    def registrar_error(self, mensaje, info_exc=False, extra_data=None):
        
        
        self.registrador.error(mensaje, exc_info=info_exc)
        if SENTRY_DSN:
            with sentry_sdk.push_scope() as scope:
                if extra_data:
                    # CORRECCIÓN CLAVE AQUÍ
                    scope.set_extras(extra_data)
                if info_exc:
                    sentry_sdk.capture_exception()
                else:
                    sentry_sdk.capture_message(f"Error: {mensaje}", level='error')

    def registrar_critico(self, mensaje, info_exc=False, extra_data=None):
        """
        Registra un mensaje crítico y, opcionalmente, envía la excepción a Sentry.
        Args:
            mensaje (str): Mensaje crítico.
            info_exc (bool): Si es True, incluye información de la excepción.
            extra_data (dict, optional): Datos adicionales para enviar a Sentry.
        """
        self.registrador.critical(mensaje, exc_info=info_exc)
        if SENTRY_DSN:
            with sentry_sdk.push_scope() as scope:
                if extra_data:
                    # CORRECCIÓN CLAVE AQUÍ
                    scope.set_extras(extra_data)
                if info_exc:
                    sentry_sdk.capture_exception()
                else:
                    sentry_sdk.capture_message(f"Crítico: {mensaje}", level='fatal')



class GestorADN:
    """
    Gestiona las operaciones de limpieza, validación, transcripción de ADN y obtención de complementarias.
    """
    def __init__(self, registrador: RegistradorSecuencias):  # <-- Cambia _init_ por __init__
        self.registrador = registrador
        self.complementos = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    def limpiar_adn(self, adn: str) -> str:
        """
        Elimina el formato 3' 5' y convierte la secuencia de ADN a mayúsculas.
        """
        adn_limpio = adn.replace("5'-","").replace("-3'","").upper()
        self.registrador.registrar_info(f"ADN limpiado: '{adn}' -> '{adn_limpio}'")
        return adn_limpio

    def transcribir_adn(self, adn: str) -> str:
        """
        Cambia las secuencias de ADN a ARN (reemplazando T por U).
        """
        arn = adn.replace('T','U')
        self.registrador.registrar_info(f"ADN transcrito a ARN: '{adn}' -> '{arn}'")
        return arn

    def validar_adn(self, adn: str) -> bool:
        """
        Verifica que la secuencia de ADN solo contenga caracteres válidos (A, T, C, G).
        """
        bases_invalidas = [base for base in adn if base not in "ATCG"]
        if bases_invalidas:
            self.registrador.registrar_advertencia(f"Validación de ADN fallida: '{adn}' contiene caracteres inválidos: {', '.join(set(bases_invalidas))}.")
            return False
        self.registrador.registrar_info(f"Validación de ADN exitosa para: '{adn}'")
        return True

    def obtener_cadena_complementaria(self, adn: str) -> str:
        """
        Genera la cadena de ADN complementaria (A-T, C-G).
        """
        if not self.validar_adn(adn):
            self.registrador.registrar_error(f"No se pudo obtener la cadena complementaria: ADN inválido '{adn}'.")
            return ""
        
        cadena_complementaria = "".join([self.complementos.get(base, '?') for base in adn])
        if '?' in cadena_complementaria:
            self.registrador.registrar_advertencia(f"Se encontraron bases desconocidas al generar la cadena complementaria para '{adn}'.")
            print(f"Advertencia: Se encontraron bases desconocidas al generar la cadena complementaria para '{adn}'.")

        self.registrador.registrar_info(f"Cadena complementaria de '{adn}' es '{cadena_complementaria}'")
        return cadena_complementaria

class AnalizadorProteinas:
    """
    Maneja la traducción de ARN a proteína y el análisis de la secuencia de proteínas.
    """
    def __init__(self, registrador: RegistradorSecuencias):  # <-- Cambia _init_ por __init__
        self.registrador = registrador

    def traducir_arn(self, arn: str) -> str:
        """
        Traduce la secuencia de ARN en una secuencia de proteínas usando la tabla de codones.
        También detecta y cuenta codones desconocidos, reportando anomalías.
        """
        proteina = ""
        conteo_codones_desconocidos = 0
        
        if len(arn) % 3 != 0:
            mensaje = f"La secuencia de ARN '{arn}' (longitud {len(arn)}) no es un múltiplo de 3. Se truncarán los últimos nucleótidos."
            self.registrador.registrar_advertencia(mensaje)
            print(f"Advertencia: {mensaje}")
            
        for i in range(0, len(arn) -2, 3):
            codon = arn[i:i+3]
            aminoacido = tabla_codones.get(codon, '?')
            if aminoacido == '?':
                conteo_codones_desconocidos += 1
            
            if aminoacido == 'STOP':
                self.registrador.registrar_info(f"Codón de parada 'STOP' encontrado en '{arn}' en posición {i}. Traducción terminada.")
                break
            proteina += aminoacido
        
        if conteo_codones_desconocidos > 0: 
            mensaje_usuario = f"¡Anomalía detectada! Se encontraron {conteo_codones_desconocidos} codones desconocidos en la secuencia de ARN. Esto podría indicar un error o una mutación."
            self.registrador.registrar_advertencia(f"Anomalía: {conteo_codones_desconocidos} codones desconocidos encontrados en la secuencia de ARN '{arn}'.")
            print(mensaje_usuario)
            if SENTRY_DSN:
                sentry_sdk.capture_message(
                    f"Anomalía Biológica: Alta cantidad de codones desconocidos ({conteo_codones_desconocidos})",
                    level='warning',
                    tags={"tipo_anomalia": "codones_desconocidos", "longitud_arn": len(arn), "conteo_desconocidos": conteo_codones_desconocidos},
                    extras={"secuencia_arn": arn, "proteina_traducida": proteina}
                )

        self.registrador.registrar_info(f"ARN '{arn}' traducido a proteína: '{proteina}'")
        return proteina

    def analizar_longitud_proteina(self, secuencia_proteina: str) -> int:
        """
        Calcula y registra la longitud de la secuencia de proteínas.
        También monitorea anomalías en la longitud de las proteínas (demasiado cortas/largas).
        """
        longitud = len(secuencia_proteina)
        self.registrador.registrar_info(f"Longitud de la proteína '{secuencia_proteina}': {longitud} aminoácidos.")
        
        if longitud < MIN_LONGITUD_PROTEINA_FUNCIONAL and longitud > 0:
            mensaje_usuario = f"¡Anomalía detectada! Proteína inusualmente corta ({longitud} aminoácidos). Esto podría afectar su función."
            self.registrador.registrar_advertencia(f"Anomalía: Proteína inusualmente corta detectada ({longitud} aminoácidos): '{secuencia_proteina}'.")
            print(mensaje_usuario)
            if SENTRY_DSN:
                sentry_sdk.capture_message(
                    f"Anomalía Biológica: Proteína inusualmente corta ({longitud} aa)",
                    level='warning',
                    tags={"tipo_anomalia": "proteina_corta", "longitud_proteina": longitud},
                    extras={"secuencia_proteina": secuencia_proteina}
                )
        elif longitud > MAX_LONGITUD_PROTEINA_FUNCIONAL:
            mensaje_usuario = f"¡Anomalía detectada! Proteína inusualmente larga ({longitud} aminoácidos). Esto podría ser un error de traducción."
            self.registrador.registrar_advertencia(f"Anomalía: Proteína inusualmente larga detectada ({longitud} aminoácidos): '{secuencia_proteina}'.")
            print(mensaje_usuario)
            if SENTRY_DSN:
                sentry_sdk.capture_message(
                    f"Anomalía Biológica: Proteína inusualmente larga ({longitud} aa)",
                    level='warning',
                    tags={"tipo_anomalia": "proteina_larga", "longitud_proteina": longitud},
                    extras={"secuencia_proteina": secuencia_proteina}
                )
        
        return longitud

    def analizar_frecuencia_aminoacidos(self, secuencia_proteina: str):
        """
        Calcula y muestra la frecuencia de cada aminoácido en una secuencia de proteína.
        """
        if not secuencia_proteina:
            print("No hay proteína para analizar la frecuencia de aminoácidos.")
            return

        frecuencia = Counter(secuencia_proteina)
        total_aminoacidos = len(secuencia_proteina)
        
        print("\n--- Frecuencia de Aminoácidos ---")
        for aminoacido, conteo in frecuencia.most_common():
            porcentaje = (conteo / total_aminoacidos) * 100
            print(f"  {aminoacido}: {conteo} ({porcentaje:.2f}%)")
        self.registrador.registrar_info(f"Análisis de frecuencia de aminoácidos completado para proteína de longitud {total_aminoacidos}.")

    def save_results_to_csv(self, adn_original: str, adn_limpio: str, secuencia_arn: str, proteina_secuencia: str, nombre_archivo='datos_adn.csv'):
        """
        Guarda los resultados del procesamiento en un archivo CSV.
        """
        datos = {
            'ADN_Original': [adn_original],
            'ADN_Limpio': [adn_limpio],
            'ARN_Secuencia': [secuencia_arn],
            'Proteina_Secuencia': [proteina_secuencia],
            'Longitud_Proteina': [len(proteina_secuencia)]
        }
        df = pd.DataFrame(datos)

        if os.path.exists(nombre_archivo):
            df.to_csv(nombre_archivo, mode='a', header=False, index=False)
            self.registrador.registrar_info(f"Resultados añadidos a '{nombre_archivo}'.")
        else:
            df.to_csv(nombre_archivo, mode='w', header=True, index=False)
            self.registrador.registrar_info(f"Archivo '{nombre_archivo}' creado y resultados guardados.")

    def cargar_y_analizar_datos(self, nombre_archivo='datos_adn.csv'):
        """
        Carga datos desde un CSV y realiza un análisis básico con pandas, presentando los resultados claramente.
        """
        if not os.path.exists(nombre_archivo):
            self.registrador.registrar_advertencia(f"El archivo '{nombre_archivo}' no existe para análisis.")
            print(f"\nNo se pudo encontrar el archivo '{nombre_archivo}' para análisis de datos. Por favor, procese algunas secuencias primero.")
            return

        try:
            df = pd.read_csv(nombre_archivo)
            self.registrador.registrar_info(f"Datos cargados desde '{nombre_archivo}'.")

            print("\n--- Análisis de Datos de Secuencias ---")
            print(f"Total de secuencias procesadas: {len(df)}")

            if not df.empty:
                print("\nEstadísticas de longitud de proteínas:")
                estadisticas_longitud_proteina = df['Longitud_Proteina'].describe()
                print(f"  Promedio: {estadisticas_longitud_proteina['mean']:.2f} aminoácidos")
                print(f"  Mínimo: {estadisticas_longitud_proteina['min']:.0f} aminoácidos")
                print(f"  Máximo: {estadisticas_longitud_proteina['max']:.0f} aminoácidos")
                print(f"  Desviación estándar: {estadisticas_longitud_proteina['std']:.2f} aminoácidos")

                # Visualización: Histograma de longitud de proteínas (requiere matplotlib/seaborn)
                # if 'matplotlib.pyplot' in globals() and 'seaborn' in globals():
                #     plt.figure(figsize=(10, 6))
                #     sns.histplot(df['Longitud_Proteina'], kde=True, bins=10)
                #     plt.title('Distribución de Longitudes de Proteínas')
                #     plt.xlabel('Longitud de Proteína (aminoácidos)')
                #     plt.ylabel('Frecuencia')
                #     plt.grid(axis='y', alpha=0.75)
                #     plt.show()
                # else:
                #     print("\nSugerencia: Instala matplotlib y seaborn para visualizar la distribución de longitudes de proteínas.")


                # Filtrar secuencias con proteinas de longitud > 10 (ejemplo de filtrado avanzado)
                print("\n--- Búsqueda de Proteínas Largas (Longitud > 10) ---")
                long_proteins = df[df['Longitud_Proteina'] > 10]
                if not long_proteins.empty:
                    print(f"Se encontraron {len(long_proteins)} secuencias con proteínas de longitud mayor a 10:")
                    for indice, fila in long_proteins.iterrows():
                        print(f"  - ADN Original: {fila['ADN_Original']} | Longitud Proteína: {fila['Longitud_Proteina']}")
                else:
                    print("No hay proteínas con longitud mayor a 10.")

                # Conteo de secuencias de ADN limpias únicas
                unique_adn_count = df['ADN_Limpio'].nunique()
                print(f"\nNúmero de secuencias de ADN limpias únicas: {unique_adn_count}")

                # Aquí podrías añadir más análisis con Pandas o incluso gráficos
                # Por ejemplo, un gráfico de barras de las bases más comunes si lo implementas

            else:
                print("\nEl archivo de datos está vacío. No hay estadísticas para mostrar.")

            print("--- Fin del Análisis ---")
            self.registrador.registrar_info("Análisis de datos con pandas completado.")

        except pd.errors.EmptyDataError:
            self.registrador.registrar_error(f"El archivo '{nombre_archivo}' está vacío.", info_exc=True)
            print(f"Error: El archivo '{nombre_archivo}' está vacío. No hay datos para analizar.")
        except Exception as e:
            self.registrador.registrar_error(f"Error al cargar o analizar el archivo '{nombre_archivo}': {e}", info_exc=True)
            print(f"Ocurrió un error inesperado al analizar los datos: {e}")
            if SENTRY_DSN: sentry_sdk.capture_exception()

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

# Al final del archivo, fuera de cualquier función:
if __name__ == "__main__":
    main()