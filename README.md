###################################
# **Proyecto glucuronidasas**
#Este proyecto tiene como objetivo presentar datos preliminares para la solicitud de plan nacional.
###################################

**Requisitos**
### Para replicar el análisis se necesitan los siguientes recursos:

- Descargarse las secuencias de NCBI estudio (https://doi.org/10.1002/1873-3468.13527)
- Obtener el catalogo de GUS (https://doi.org/10.1016/j.str.2017.05.003
- Modules: FastQC, Trimmomatic, Bowtie2, DIAMOND, R

### Estructura del Pipline
1. Descargar archivos brutos. 
2. Limpieza de las secuencias. (01_QC_host_removal_SE.sbatch|02_ResumenQC.sbatch)
3. Preparación del catálogo.(03_DIAMOND_sb_GUS.sbatch
4. Mapear las lecturas frente al catalogo GUS. (04_GUS_blastx_dm.sbatch)
5. Añadir anotaciones de clase. (05_add_GUS_annot.sbatch) 
6. Filtrar los mapeos débiles y seleccionar un uno hit por lectura.(06_filter_all_GUS.sbatch)
7. Evaluar diferencias de abundancia por grupos.(07_count_GUS.sbatch | 08_GUSclass_abudance.sbatch | 09_GUS_abundance.sbatch)
