# post_expresion_diferenicial_sleuth
# compara los resultados de la expresión diferencial con Sleuth con la tabla uniprot recortada de Annocript

qVal <-  0.01 # threshold de q valor 

###### CARGA LA TABLA MAESTRA DE BLAST VS UNIPROT
uprot_anno_file <-read.csv2("anno_uniprot_exp_dif_lncRNA.csv2", stringsAsFactors = FALSE)
 
###### EXAMINA FICHERO, CARGA TODAS LAS TABLAS DE SLEUTH COMO DATA FRAMES
base_dir <- paste0(getwd(),"/tablas_sleuth") # directorio de los output de sleuth
sample_id <-dir(file.path(base_dir)) # recupera los nombres de los archivos en ese directorio

# Carga todas las tablas de sleuth contenidas en el directiorio "tablas_sleuth y les hace una variable con 
# el nombre del archivo" 
for (i in 1:length(sample_id)){
  assign(paste(sample_id[i]), # Crea variable con el nombre del archivo
         read.csv(paste0(base_dir,"/", sample_id[i]),stringsAsFactors = FALSE, header = TRUE))
  # (eval(parse(text=sample_id[1]))) esto llama a la variable desde la lista!!!!
  assign(paste(sample_id[i]), subset((eval(parse(text=sample_id[i]))), qval <= qVal)) # filtra por qValue
  uniprot_rows<-which(uprot_anno_file$ID_mx_corto %in% paste0(">",eval(parse(text=sample_id[i]))$target_id))
  x <- uprot_anno_file[uniprot_rows,]
  assign(paste(sample_id[i]), with(eval(parse(text=sample_id[i])),eval(parse(text=sample_id[i]))[order(target_id),]))
  x <- with(x, x[order(ID_mx_corto),]) # esta linea y la de arriba ordenan por nombre de transcrito
  assign(paste(sample_id[i]), cbind.data.frame(eval(parse(text=sample_id[i])) , x)) # une las dos tablas en una sola
  assign(paste(sample_id[i]), with(eval(parse(text=sample_id[i])),eval(parse(text=sample_id[i]))[order(qval),]))
  write.csv(eval(parse(text=sample_id[i])),paste0("uref_",sample_id[i],".csv"))
}